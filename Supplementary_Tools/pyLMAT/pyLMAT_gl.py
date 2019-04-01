#!/usr/bin/env python3
#
# Python 3.X site-specific launcher for LMAT run_gl.sh (gene labeling step)
#
# NOTE: the LMAT module should be loaded! Tested module configuration:
# module purge; module load gcc/gcc-4.9.0 python/3.4.0 LMAT/1.2.4

import argparse
import os
import sys
import errno
import subprocess

# pyLMAT release information
__version__ = '0.0.4'
_verdata = 'Mar 2015'
_devflag = True

# Predefined internal constants
_PATH = './kwashiorkor/'
_FLST = 'rl_output.flst'
_GLDB = '/disk/disksom2/ramdisk/lmat.genes.7-14.db'
_MTSC = 'auto'

# Argument Parser Configuration
parser = argparse.ArgumentParser(
    description='pythonic LMAT launcher',
    epilog='pyLMAT_gl - DLS team - by J.Mn.Marti - ' + _verdata
)
parser.add_argument(
    '-V', '--version',
    action='version',
    version='pyLMAT_gl.py release ' + __version__ + ' - ' + _verdata
)
parser.add_argument(
    '-p', '--path',
    action='store',
    default=_PATH,
    help=('relative path of the data files (if omitted, \'' +
          _PATH + '\' will be tried)')
)
parser.add_argument(
    '-d', '--db',
    action='store',
    default=_GLDB,
    help=('absolute path of the gene database (if omitted, \'' +
          _GLDB + '\' will be tried)')
)
parser.add_argument(
    '-m', '--min_tax_score',
    action='store',
    default=_MTSC,
    help=('minimum score for matched tax id, where \'auto\' will try to ' +
          ' guess it from the fastsummary filename (if omitted, \'' +
          _MTSC + '\' will be tried)')
)

# Parse arguments (except min_tax_score)
args = parser.parse_args()
path = args.path
db = args.db

# Program Header
print('\n=-= pyLMAT_gl =-= v' + __version__ + ' =-= ' +
      _verdata + ' =-= by DLS team =-=')
if(_devflag):
    print('\n>>> WARNING! THIS IS JUST A DEVELOPMENT SUBRELEASE.' +
          ' USE AT YOUR OWN RISK!')

# LMAT script and options
_pgrm = 'run_gl.sh'
_dbfl = '--db_file='
_ilst = '--ilst='
_fsum = '--filesum='
_odir = '--odir='
_mtsc = '--min_tax_score='
_ovw = '--overwrite'

(root, dirs, files) = next(os.walk(path))
# Main loop: every subdir in idir
for d in dirs:
    (root_, dirs_, files_) = next(os.walk(path + '/' + d))
    # 2ndary loop: launch run_cs in processed subdirs
    filelst = []
    numfilelst = 0
    flst = root_ + '/' + _FLST
    mtsc = args.min_tax_score
    for f in files_:
        if f.endswith('.out') and not f.startswith(_FLST):
            filelst.append(root_ + '/' + f + '\n')  # EOL for 'writelines'
            numfilelst += 1
        if f.endswith('.fastsummary') and '.out.' not in f:
            sumf = root_ + '/' + f
            if mtsc == 'auto':
                try:
                    mtsc = f.split('.')[-3]
                    dumbvar = int(mtsc)  # Exception if mtsc is not a number
                except:
                    mtsc = '0'
                    print('Failed guessing the min_tax_score! Setting it to ' +
                          mtsc)
    if len(filelst):
        print(str(numfilelst) + ' .out files added to ' + flst)
        with open(flst, 'w') as lstfile:
            lstfile.writelines(filelst)
        argums = [_pgrm, _dbfl + db, _ilst + flst, _fsum + sumf, _odir + root_,
                  _mtsc + mtsc, _ovw]
        print('\tLaunching subprocess: ', argums)
        sys.stdout.flush()
        try:
            p = subprocess.call(argums)
        except (PermissionError):
            print('\nERROR! Unable to launch ' + repr(_pgrm))
            print('TIP: Check if the *right* LMAT module is loaded!\n')
            raise
