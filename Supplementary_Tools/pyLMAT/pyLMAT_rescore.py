#!/usr/bin/env python3
#
# Python 3.X site-specific launcher for LMAT losummary_fast_mc.sh and further
#
# NOTE: the LMAT module should be loaded! Tested module configuration:
# module purge; module load gcc/gcc-4.9.0 python/3.4.0 LMAT/1.2.4

import argparse
import os
import sys
import subprocess as sp

# pyLMAT release information
__version__ = '0.0.8'
_scriptname = 'pyLMAT_rescore'
_verdata = 'Mar 2015'

# Predefined internal constants
_PATH = './kwashiorkor/'
_TAXID = '9606'


def printcall(arglst, stdout=None, stderr=sp.STDOUT, shell=False, minerr=0):
    """Aux method: Print and manage an external call"""
    if type(arglst) is not list:
        arglst = [arglst]
    print('\033[90mRUNNING\033[0m:', ' '.join(arglst), end='')
    sys.stdout.flush()
    try:
        retcode = sp.call(arglst, stdout=stdout, stderr=stderr, shell=shell)
    except (PermissionError, FileNotFoundError):
        print('\n\033[91m ERROR! \033[0m Unable to launch ' + arglst[0])
        print('\033[34m TIP\033[0m: Check if the *right* ' +
              'LMAT module is loaded!\n')
        raise
    except:
        print('\n\033[91m ERROR! \033[0m Unable to launch ' + arglst[0])
        raise
    else:
        if retcode > minerr:
            print('\n\033[91m ERROR! \033[0m Failed running of ' + arglst[0])
        else:
            print(' ...\033[92m OK! \033[0m')

# Argument Parser Configuration
parser = argparse.ArgumentParser(
    description='pythonic LMAT losummary_fast_mc.sh&Co launcher',
    epilog=_scriptname + ' - DLS team - by J.Mn.Marti - ' + _verdata
)
parser.add_argument(
    '-v', '--version',
    action='version',
    version=_scriptname + '.py release ' + __version__ + ' - ' + _verdata
)
parser.add_argument(
    '-p', '--path',
    default=_PATH,
    help=('relative path of the sequence file (if omitted, \'' +
          _PATH + '\' will be tried)')
)
parser.add_argument(
    '-f', '--file',
    required=True,
    help=('target filename')
)
parser.add_argument(
    '-t', '--threads',
    default='64',
    help=('number of OpenMP threads to use (64 by default)')
)
parser.add_argument(
    '-m', '--minscore',
    nargs='+',
    default='0',
    help=('minimum score assigned to read for it to be included in binning' +
          ' (0 by default, but example list: -1 0 1)')
)
parser.add_argument(
    '-q', '--freq',
    default='10',
    help=('minimum frequency for taxa to be included in Krona plot' +
          ' (10 by default)')
)
parser.add_argument(
    '-x', '--taxid',
    default=_TAXID,
    help=('tax id for pulling from reads (if omitted, \'' +
          _TAXID + '\' will be tried)')
)

# Parse arguments
args = parser.parse_args()
path = args.path
mfrq = args.freq
file = args.file
txid = args.taxid
lstmsco = [msco.strip() for msco in args.minscore]

# Program Header
print('\n=-= ' + _scriptname + ' =-= v' + __version__ + ' =-= ' +
      _verdata + ' =-= by DLS team =-=\n')

# LMAT scripts and options
_stdo = _scriptname + '.' + file + '.log'
_p_lo = 'losummary_fast_mc.sh'
_p_fs = 'fsreport.py'
_p_to = 'tolineage.py'
_p_kt = 'ktImportText'
_p_pu = 'pull_reads_mc.sh'
_flst = '--file_lst=rl_output.lst'
_nthr = '--threads=' + args.threads
_taxl = 'plasmid,species,genus'
_ifil = '--idfile=tid.txt'
_ncbi = ('/software/LMAT-1.2.4/runtime_inputs/' +
         'ncbi_taxonomy_rank.segment.pruned.txt')

# Main
os.chdir(path)
with open(_stdo, 'at') as stdo:
    print('\033[34m NOTE\033[0m: Details logged in', path + _stdo)
    print('\n\033[36m>>>\033[0m Initialization steps:')
    arg = 'ls -1 *output*.out > rl_output.lst'
    printcall(arg, stdout=stdo, shell=True)
    ftar = file + '.tar'
    # tar doesn't remove old versions, so delete entire file
    printcall(['rm', '-vf', ftar], stdout=stdo)
    for msco in lstmsco:
        print('\n\033[36m>>>\033[0m Obtaining results for minimum \033[36m' +
              'score', msco, '\033[0m:')
        _mins = '--min_score=' + msco
        printcall([_p_lo, _flst, _nthr, _mins], stdout=stdo)
        fastsum = file + '.lmat-4-14.20mer.db.lo.rl_output.' + \
            msco + '.30.fastsummary'
        printcall([_p_fs, fastsum, _taxl, '.'], stdout=stdo)
        flineage = fastsum + '.lineage'
        printcall([_p_to, _ncbi, fastsum, flineage, mfrq, msco], stdout=stdo)
        fhtml = flineage + '.html'
        printcall([_p_kt, flineage, '-o', fhtml], stdout=stdo)
        printcall('echo ' + txid + ' > tid.txt', stdout=stdo, shell=True)
        printcall([_p_pu, _ifil, _flst, _mins], stdout=stdo)
        fsco = file + '.lmat-4-14.20mer.db.tid.txt.' + txid + '.fna'
        fsto = file + '.lmat-4-14.20mer.db.minsco' + msco + '.' + txid + '.fna'
        printcall(['mv', '-vf', fsco, fsto], stdout=stdo)
        tarsum = fastsum + '*'
        # tar exit code 1 means no error but 'Some files differ'
        printcall('tar uvf ' + ftar + ' ' + tarsum + ' ' + fsto,
                  stdout=stdo, shell=True, minerr=1)
    print('\n\033[36m>>>\033[0m Finalization steps:')
    printcall(['gzip', '--best', '--force', ftar], stdout=stdo)
