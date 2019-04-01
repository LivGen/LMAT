#!/usr/bin/env python3
#
# Python 3.X site-specific launcher for LMAT run_cs.sh (after run_rl.sh step)
#
# NOTE: the LMAT module should be loaded! Tested module configuration:
# module purge; module load gcc/gcc-4.9.0 python/3.4.0 LMAT/1.2.4

import argparse
import os
import sys
import errno
import subprocess

# pyLMAT release information
__version__ = '0.0.3'
_verdata = 'Jan 2015'
_devflag = True

# > MAIN

# Predefined internal constants
_PATH = './kwashiorkor/'
_FLST = 'rl_output.flst'

# Argument Parser Configuration
parser = argparse.ArgumentParser(
    description='pythonic LMAT launcher',
    epilog='pyLMAT_cs - DLS team - by J.Mn.Marti - ' + _verdata
)
parser.add_argument(
    '-v', '--version',
    action='version',
    version='pyLMAT.py release ' + __version__ + ' - ' + _verdata
)
parser.add_argument(
    '-p', '--path',
    action='store',
    default=_PATH,
    help=('relative path of the data files (if omitted, \'' +
          _PATH + '\' will be tried)')
)

# Parse arguments
args = parser.parse_args()
path = args.path

# Program Header
print('\n=-= pyLMAT_cs =-= v' + __version__ + ' =-= ' +
      _verdata + ' =-= by DLS team =-=')
if(_devflag):
    print('\n>>> WARNING! THIS IS JUST A DEVELOPMENT SUBRELEASE.' +
          ' USE AT YOUR OWN RISK!')

# LMAT script and options
_pgrm = 'run_cs.sh'
_ilst = '--ilst='
_fsum = '--filesum='
_odir = '--odir='
_ovw = '--overwrite'

(root, dirs, files) = next(os.walk(path))
# Main loop: every subdir in idir
for d in dirs:
    (root_, dirs_, files_) = next(os.walk(path + '/' + d))
    # 2ndary loop: launch run_cs in processed subdirs
    filelst = []
    flst = root_ + '/rl_output.flst'
    sumf = root_ + '/*.fastsummary'
    for f in files_:
        if f.endswith('.out'):
            filelst.append(root_ + '/' + f + '\n')  # EOL for 'writelines'
    if len(filelst):
        with open(flst, 'w') as lstfile:
            lstfile.writelines(filelst)

    argums = [_pgrm, _ilst + flst, _fsum + sumf, _odir + root_, _ovw]
    print('\tLaunching subprocess: ', argums)
    sys.stdout.flush()
    try:
        p = subprocess.call(argums)
    except (PermissionError):
        print('\nERROR! Unable to launch ' + repr(_pgrm))
        print('TIP: Check if the *right* LMAT module is loaded!\n')
        raise
