#!/usr/bin/env python3
"""Python3 converter from LMAT-1.2.X results to CC input."""

import argparse
import os
import subprocess
import pandas as pd

# pyLMAT release information
__version__ = '0.1.2'
_verdata = 'Dec 2016'
_devflag = True

# Predefined internal constants
_PATH = './LMAT/'


class LMAT2CC(object):

    """Convert multiple-times LMAT files to multiple-sheet Excel file."""

    def __init__(self, args):
        """Class initialization."""
        # Attributes
        self.datapath = args.path
        self.join = args.join
        self.verb = (args.verbose if args.verbose is not None else 0)
        try:
            os.chdir(self.datapath)
        except (OSError):
            raise Exception(
                'ERROR: Cannot change to directory \'' + self.datapath + '\'')

    def parse(self, step):
        """Parse the fastsummary files resulting from a LMAT step."""
        if (step == '1'):
            self.xlsx_gen = pd.ExcelWriter('LMAT_RL_gen2cmplx.xlsx')
            self.xlsx_spc = pd.ExcelWriter('LMAT_RL_spc2cmplx.xlsx')
            self.file_gen = 'fastsummary.genus'
            self.file_spc = 'fastsummary.species'
        elif (step == '2'):
            self.xlsx_gen = pd.ExcelWriter('LMAT_CS_gen2cmplx.xlsx')
            self.xlsx_spc = pd.ExcelWriter('LMAT_CS_spc2cmplx.xlsx')
            self.file_gen = 'fastsummary.ordered.genus'
            self.file_spc = 'fastsummary.ordered.species'
        elif (step == '3'):
            raise Exception('ERROR: LMAT step 3 is not yet implemented here!')
        else:
            raise Exception('ERROR: Unknown LMAT step!')
        self.processedcounter = 0
        print('Processing step', step, 'results under', self.datapath, '...')
        # Main loop
        for root, dirs, files in os.walk('.'):
            if self.verb > 0:
                print('Looking for target files in: ', root)
            genusfiles = []
            speciesfiles = []
            for file in files:
                if file.endswith(self.file_gen):
                    genusfiles.append(file)
                elif file.endswith(self.file_spc):
                    speciesfiles.append(file)
            if len(genusfiles):
                genusfiles.sort()
                if self.verb > 1:
                    print('   Parsing GENUS files...')
                self.summary2cmplx(step, 'genus', root, genusfiles)
            if len(speciesfiles):
                speciesfiles.sort()
                if self.verb > 1:
                    print('   Parsing SPECIES files...')
                self.summary2cmplx(step, 'species', root, speciesfiles)
        # Finish
        if (self.processedcounter):
            print('Saving final spreadsheets under', self.datapath, '...')
            self.xlsx_gen.save()
            self.xlsx_spc.save()
            print('-> lmat2cmplx successfully converted ' +
                  str(int(self.processedcounter / 2)) + ' step-' + step +
                  ' datasets under the dir "' + self.datapath + '".')
        else:
            print('-> lmat2cmplx was not able to parse any step-' +
                  step + ' file under "' + self.datapath + '"')

    def summary2cmplx(self, step, level, root, summaryfiles):
        """Convert multiple-times LMAT-fastsummary files to Excel sheet."""
        DFlist = []
        for filename in summaryfiles:
            if (filename[0] == '_'):
                if self.verb > 2:
                    print('\tNOT processing file:', filename)
                else:
                    pass
            else:
                if self.verb > 2:
                    print('\tParsing file: ', filename)
                fullfilename = os.path.join(root, filename)
                fileNOext = os.path.splitext(filename)[0]
                filesplit = fileNOext.split('.')
                dataset = filesplit[0]
                time = filesplit[1]
                if (step == '1'):
                    DF = pd.read_table(fullfilename,
                                       engine='c',
                                       names=['Average Read Score', time,
                                              'Read Count', 'TaxID', level,
                                              'Strain ID', 'Strain Name'],
                                       header=None,
                                       index_col=4
                                       )
                    DF.drop(['Average Read Score', 'Read Count', 'TaxID',
                             'Strain ID', 'Strain Name'],
                            axis='columns',
                            inplace=True)
                    try:
                        DF.drop('Name', axis='index', inplace=True)
                    except:
                        print('WARNING: Title row was not deleted!')
                elif (step == '2'):
                    DF = pd.read_table(fullfilename,
                                       engine='c',
                                       names=['% of Reads', 'Avg Read Score',
                                              time, 'Read Count (RC)',
                                              'Original WRC', 'Original RC',
                                              level, 'Taxid', 'X1', 'X2', 'X3',
                                              'X4', 'X5', 'X6', 'X7', 'X8'],
                                       header=None,
                                       comment='%',
                                       index_col=6
                                       )
                    DF.drop(['% of Reads', 'Avg Read Score', 'Read Count (RC)',
                             'Original WRC', 'Original RC', 'Taxid',
                             'X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8'],
                            axis='columns',
                            inplace=True)
                    # try:
                    #    DF.drop('% of Reads', axis='index', inplace=True)
                    # except:
                    #    print('WARNING: Title row was not deleted!')
                else:
                    raise Exception('summary2cmplx ERROR: Unknown LMAT step!')
                if level == 'species':
                    try:
                        DF.drop('synthetic construct',
                                axis='index', inplace=True)
                    except:
                        pass
                DFlist.append(DF)
        if (DFlist):
            DFall = pd.concat(
                DFlist, axis=1,
                verify_integrity=True, join=self.join
            )
            DFall.fillna(value=0, method=None, inplace=True)
            if (level == 'genus'):
                DFall.to_excel(self.xlsx_gen, sheet_name=dataset)
            elif (level == 'species'):
                DFall.to_excel(self.xlsx_spc, sheet_name=dataset)
            else:
                raise Exception('ERROR: Unknown fastsummary taxonomic-level!')
            self.processedcounter += 1
        else:
            pass


# > MAIN

if __name__ == '__main__':
    # Argument Parser Configuration
    parser = argparse.ArgumentParser(
        description='LMAT-1.2.4 fastsummary results to CC input',
        epilog='lmat2cmplx - DLS team - by J.Mn.Marti - ' + _verdata
    )
    parser.add_argument(
        '-V', '--version',
        action='version',
        version='lmat2cmplx release ' + __version__ + ' - ' + _verdata
    )
    parser.add_argument(
        '-v', '--verbose',
        action='count',
        help='increase output verbosity (-v, -vv or even -vvv)'
    )
    parser.add_argument(
        '-p', '--path',
        action='store',
        default=_PATH,
        help=('path of the LMAT fastsummary results (if omitted, \'' +
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
    print('\n=-= lmat2cmplx =-= v' + __version__ + ' =-= ' +
          _verdata + ' =-= by DLS team =-=')
    if(_devflag):
        print('\n>>> WARNING! THIS IS JUST A DEVELOPMENT SUBRELEASE.' +
              ' USE AT YOUR OWN RISK!\n')

    conv = LMAT2CC(args)
    conv.parse('1')
    conv.parse('2')
