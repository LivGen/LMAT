#!/usr/bin/env python3
#
# Python 3.X site-specific launcher for LMAT run_rl.sh (1st step)
#
# NOTE: the LMAT module should be loaded! Tested module configuration:
# module purge; module load gcc/gcc-4.8.2 python/3.3.4rc1 LMAT/1.2.1

import argparse
import os
import sys
import errno
import subprocess

# pyLMAT release information
__version__ = '0.0.10'
_verdata = 'Mar 2015'
_devflag = True

# Predefined internal constants
_PATH = './kwashiorkor/'
_BDIR = '/disk/disksom2/martijm/'
_FULL = '/disk/disksom2/ramdisk/'
# _BDIR = '/home/martijm/lmat/'
# _FULL = '/ramdisk/'

# Argument Parser Configuration
parser = argparse.ArgumentParser(
    description='pythonic LMAT launcher',
    epilog='pyLMAT_rl - DLS team - by J.Mn.Marti - ' + _verdata
)
parser.add_argument(
    '-v', '--version',
    action='version',
    version='pyLMAT_rl.py release ' + __version__ + ' - ' + _verdata
)
parser.add_argument(
    '-p', '--path',
    action='store',
    default=_PATH,
    help=('relative path of the data files (if omitted, \'' +
          _PATH + '\' will be tried)')
)
parser.add_argument(
    '-b', '--bdir',
    action='store',
    metavar='PATH',
    default=_BDIR,
    help=('base directory (if omitted, \'' +
          _BDIR + '\' will be tried)')
)
parser.add_argument(
    '-f', '--fulldbdir',
    action='store',
    metavar='PATH',
    default=_FULL,
    help=('path of the full LMAT database lmat-4-14.20mer.db (if omitted, \'' +
          _FULL + '\' will be tried)')
)
parser.add_argument(
    '-t', '--threads',
    action='store',
    default='64',
    help=('number of OpenMP threads to use (64 by default)')
)
parser.add_argument(
    '-m', '--minscore',
    action='store',
    default='0',
    help=('minimum score assigned to read for it to be included in binning' +
          ' (0 by default)')
)
# Parser group: sequencing
groupSeq = parser.add_mutually_exclusive_group(required=True)
groupSeq.add_argument(
    '-w', '--wgs',
    action='store_true',
    default=True,
    help=('suppose Whole Genome Shotgun sequencing')
)
groupSeq.add_argument(
    '-s', '--s16',
    action='store_true',
    default=False,
    help=('suppose 16S sequencing instead of WGS')
)

# Parse arguments
args = parser.parse_args()
path = args.path
bdir = args.bdir
s16 = '_16S' if args.s16 else ''

# Program Header
print('\n=-= pyLMAT_rl =-= v' + __version__ + ' =-= ' +
      _verdata + ' =-= by DLS team =-=')
if(_devflag):
    print('\n>>> WARNING! THIS IS JUST A DEVELOPMENT SUBRELEASE.' +
          ' USE AT YOUR OWN RISK!')

# LMAT script and options
_pgrm = 'run_rl.sh'
_db = '--db_file=' + args.fulldbdir + '/lmat-4-14.20mer.db'
_fst = '--query_file='
_odir = '--odir='
_nthr = '--threads=' + args.threads
_ovw = '--overwrite'
_mins = '--min_score=' + args.minscore
_null = ('--nullm=/software/LMAT-1.2.4/runtime_inputs/' +
         'lmat-4-14.20mer.16bit.g200.grand.with-adaptors.db.null_lst.txt')

# Main loop
for root, dirs, files in os.walk(path):
    for file in files:
        argums = [_pgrm, _db, _mins, _nthr, _null]
        fileNOext = os.path.splitext(file)[0]
        outdir = bdir + path + '/' + fileNOext + s16
        print('\nStoring results in the dir: ', outdir)
        try:
            os.makedirs(outdir)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        argums.append(_fst + path + '/' + file)
        argums.append(_odir + outdir)
        print('\tLaunching subprocess: ', argums)
        sys.stdout.flush()
        try:
            p = subprocess.call(argums)
        except (PermissionError):
            print('\nERROR! Unable to launch ' + repr(_pgrm))
            print('TIP: Check if the *right* LMAT module is loaded!\n')
            raise
        print('Copying krona HTML from the dir ' + outdir + ' to ~/krona')
        sys.stdout.flush()
        try:
            c1 = subprocess.call('cp -v ' + outdir + '/' + fileNOext +
                                 '.*.lineage.html ' + '~/krona/' + fileNOext +
                                 s16 + '.lineage.html', shell=True)
# print('Copying summary HTML from the dir ' + outdir + ' to ~/krona')
# c2 = subprocess.call('cp -v ' + outdir + '/' + fileNOext +
#                      '.*.fastsummary.html ' + '~/krona/' + fileNOext + s16 +
#                      '.fastsummary.html', shell=True)
        except:
            print('WARNING! The copy of the HTML files failed for: ', file)
