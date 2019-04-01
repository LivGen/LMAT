#!/usr/bin/env python3.4
"""Python 3.X converter from LMAT-1.2.X gene labeling results to cmplxcruncher input."""

import argparse
import os
import subprocess
import pandas as pd

# pyLMAT release information
__version__ = '0.0.1'
_verdata = 'Dec 2014'
_devflag = True

# Predefined internal constants
_PATH = './LMAT/'


class LMATgl2CC(object):

    """Convert multiple-times LMAT files to multiple-sheet Excel file."""

    def __init__(self, args):
        """Class initialization."""
        # Attributes
        self.datapath = args.path
        self.join = args.join
        try:
            os.chdir(self.datapath)
        except (OSError):
            raise Exception(
                'ERROR: Cannot change to directory \'' + self.datapath + '\'')
        self.xlsx_gene = pd.ExcelWriter('LMATgene2cmplx.xlsx')
        self.processedcounter = 0          
        # Main loop
        for root, dirs, files in os.walk('.'):
            print('Looking for target files in: ',root)
            genefiles=[]
            for file in files:
                if file.endswith('.genesummary'):
                    genefiles.append(file)
            if len(genefiles):
                genefiles.sort()
                print('   Parsing GENESUMMARY files...')
                self.summary2cmplx('gene', root, genefiles)
        # Finish
        if (self.processedcounter):
            self.xlsx_gene.save()
            print('END: LMATgl2cmplx successfully converted ' +
                str(self.processedcounter) + ' datasets under the dir \'' + 
                self.datapath+'\'.')
        else:
            print('END: LMATgl2cmplx was not able to parse any file under \'' + 
                self.datapath+'\'')
            
    def summary2cmplx(self, target, root, summaryfiles):
        """Convert multiple-times LMAT-genesummary files to Excel sheet."""
        DFlist = []
        for filename in summaryfiles:
            if (filename[0] == '_'):
                print('\tNOT processing file:', filename)
            else:
                print('\tParsing file: ', filename)
                fullfilename = os.path.join(root, filename)
                fileNOext = os.path.splitext(filename)[0]
                filesplit = fileNOext.split('.')
                dataset = filesplit[0]
                time = filesplit[1]
                DF = pd.read_table(
                    fullfilename, 
                    usecols=[4,1],
                    names=[time, 'gene'], 
                    header=None, 
                    index_col=1
                    )
                #DF=DF.sort_index(by=time,ascending=False)[:99]
                DF=DF.groupby(level=0).sum()
                #print(DF)
                DFlist.append(DF)
        DFall = pd.concat(
            DFlist, axis=1,
            verify_integrity=True, join=self.join
            )
        DFall.fillna(value=0, method=None, inplace=True)
        DFall['SUM']=DFall.sum(axis=1)
        DFall=DFall.sort_index(by='SUM',ascending=False)[:100]
        DFall.drop('SUM',axis=1,inplace=True)
        if (target == 'gene'):
            DFall.to_excel(self.xlsx_gene, sheet_name=dataset)
        else:
            raise Exception('ERROR: Unknown target!')
        self.processedcounter += 1


#########################################################################> MAIN
#_GENUS_COLS = ['AvrReadScore', 'WeightedReadCount',
#        'AssignedReads', 'NCBItaxNumber']
#_DF_COLS = ['WeightedReadCount']
#_DEL_COLS = _GENUS_COLS.copy()
#for column in _DF_COLS: _DEL_COLS.remove(column) 

if __name__ == '__main__':
    # Argument Parser Configuration
    parser = argparse.ArgumentParser(
        description='LMAT-1.2.x genesummary results to cmplxcruncher input',
        epilog='LMATgl2cmplx - DLS team - by J.Mn.Marti - ' + _verdata
        )
    parser.add_argument(
        '-V', '--version',
        action='version',
        version='LMATgl2cmplx.py release ' + __version__ + ' - ' + _verdata
        )
    parser.add_argument(
        '-p', '--path', 
        action='store',
        default=_PATH,
        help=('path of the LMAT genesummary results (if omitted, \'' + 
              _PATH + '\' will be tried)')
        )
    parser.add_argument(
        '-j', '--join', 
        action='store',
        choices={'outer', 'inner'},
        default='outer',
        help=('how to handle indexes on the elements concatenation ' +
            'of time series: outer for union and inner for intersection' +
            ' (outer is the default)')
        )
    # Parse arguments
    args = parser.parse_args()

    # Program Header
    print('\n=-= LMATgl2cmplx =-= v' + __version__ + ' =-= ' +
          _verdata + ' =-= by DLS team =-=')
    if(_devflag): 
        print('\n>>> WARNING! THIS IS JUST A DEVELOPMENT SUBRELEASE.' + 
              ' USE AT YOUR OWN RISK!\n')  

    LMATgl2CC(args)

