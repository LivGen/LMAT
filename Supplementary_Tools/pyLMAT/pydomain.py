#!/usr/bin/env python3
"""
pythonic LMAT multi-domain assignment analyzer

"""

import argparse
from collections import Counter
import importlib.util
import multiprocessing as mp
import os
import shutil
import sys
import time

#
# Biopython initialization
#
from Bio import SeqIO
# import lmat_io from Recentrifuge specific location
spec = importlib.util.spec_from_file_location(
    'lmat_io', '/home_local/martijm/recentrifuge/recentrifuge/lmat_io.py')
lmat_io = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lmat_io)
# Addition of new custom formats to SeqIO
# pylint: disable=protected-access
SeqIO._FormatToIterator["lmat"] = lmat_io.lmat_out_iterator
SeqIO._FormatToWriter["lmat"] = lmat_io.LmatOutWriter
# pylint: enable=protected-access

#
# Matplotlib initialization
#
import matplotlib as mpl
# Set right MPL back-end before any other MPL import
if os.environ.get('DISPLAY', '') == '':
    print('Warning! No display found. Defaulting to Agg backend')
    mpl.use('Agg')
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt

# Script release information
__version__ = '0.2.0'
__author__ = 'Jose Manuel Martí'
__date__ = 'Nov 2018'

# Constants
DOMAIN_TPL = ('___V', 'abeV', 'abEv', 'abEV', 'aBev', 'aBeV', 'aBEv', 'aBEV',
              'Abev', 'AbeV', 'AbEv', 'AbEV', 'ABev', 'ABeV', 'ABEv', 'ABEV',
              '____U', 'abevU', 'over', 'else')
DIFTYPE_TPL = ('canVfinA', 'canVfinB', 'canVfinE')


def procLmatOut(*args, **kwargs):
    """Process a LMAT .out file (to be usually called in parallel!)."""
    lmatoutfile = args[0]
    lmatoutnum = (os.path.splitext(lmatoutfile)[0]).split('output')[-1]
    seqs_lst = []
    # Initialize Domain presence lists
    seqs_dic = dict()
    for domtype in DOMAIN_TPL:
        seqs_dic[domtype] = []
    # Initialize Domain differences lists
    canVfinA_seqs = []
    canVfinB_seqs = []
    canVfinE_seqs = []
    with open(lmatoutfile, 'rU') as lmo_handle:
        seqs_lst = list(SeqIO.parse(lmo_handle, "lmat"))
    for seq in seqs_lst: # Main loop
        finalcallid = {seq.annotations['final_taxid']}
        candidset = set(seq.annotations['candidict'].keys())
        tags = seq.annotations['tags']
        # Tags by candidate list
        if arch_set & candidset:
            tags |= {'candi_arch', 'arch'}
        if bact_set & candidset:
            tags |= {'candi_bact', 'bact'}
        if virs_set & candidset:
            tags |= {'candi_virs', 'virs'}
        if euka_set & candidset:
            tags |= {'candi_euka', 'euka'}
        if unkn_set & candidset:
            tags |= {'candi_unkn', 'unkn'}
        if over_set & candidset:
            tags |= {'candi_over', 'over'}
        # Tags by final call
        if arch_set & finalcallid:
            tags |= {'final_arch', 'arch'}
        elif bact_set & finalcallid:
            tags |= {'final_bact', 'bact'}
        elif virs_set & finalcallid:
            tags |= {'final_virs', 'virs'}
        elif euka_set & finalcallid:
            tags |= {'final_euka', 'euka'}
        elif unkn_set & finalcallid:
            tags |= {'final_unkn', 'unkn'}
        elif over_set & finalcallid:
            tags |= {'final_over', 'over'}
        else:
            tags |= {'final_else', 'else'}
        # Populate seq lists by Domain presence
        if tags >= {'arch', 'bact',' euka', 'virs'}:
            seqs_dic['ABEV'].append(seq)
        #
        if tags >= {'arch', 'bact', 'euka'} and 'virs' not in tags:
            seqs_dic['ABEv'].append(seq)
        if tags >= {'arch', 'bact', 'virs'} and 'euka' not in tags:
            seqs_dic['ABeV'].append(seq)
        if tags >= {'arch', 'euka', 'virs'} and 'bact' not in tags:
            seqs_dic['AbEV'].append(seq)
        if tags >= {'bact', 'euka', 'virs'} and 'arch' not in tags:
            seqs_dic['aBEV'].append(seq)
        #
        if tags >= {'euka', 'virs'} and not tags & {'arch', 'bact'}:
            seqs_dic['abEV'].append(seq)
        if tags >= {'bact', 'virs'} and not tags & {'arch', 'euka'}:
            seqs_dic['aBeV'].append(seq)
        if tags >= {'arch', 'virs'} and not tags & {'bact', 'euka'}:
            seqs_dic['AbeV'].append(seq)
        if tags >= {'bact', 'euka'} and not tags & {'arch', 'virs'}:
            seqs_dic['aBEv'].append(seq)
        if tags >= {'arch', 'bact'} and not tags & {'euka', 'virs'}:
            seqs_dic['ABev'].append(seq)
        if tags >= {'arch', 'euka'} and not tags & {'bact', 'virs'}:
            seqs_dic['AbEv'].append(seq)
        #
        if tags & {'arch', 'bact', 'euka', 'virs'} == {'virs'}:
            seqs_dic['abeV'].append(seq)
        if tags & {'arch', 'bact', 'euka', 'virs'} == {'euka'}:
            seqs_dic['abEv'].append(seq)
        if tags & {'arch', 'bact', 'euka', 'virs'} == {'bact'}:
            seqs_dic['aBev'].append(seq)
        if tags & {'arch', 'bact', 'euka', 'virs'} == {'arch'}:
            seqs_dic['Abev'].append(seq)
            seqs_dic['aBev'].append(seq)
        if tags & {'arch', 'bact', 'euka', 'virs', 'unkn'} == {'unkn'}:
            seqs_dic['abevU'].append(seq)
        if tags & {'arch', 'bact', 'euka', 'virs', 'unkn', 'over'} == {'over'}:
            seqs_dic['over'].append(seq)
        if tags & {'arch', 'bact', 'euka', 'virs', 'unkn', 'over', 'else'} == {
            'else'}:
                seqs_dic['else'].append(seq)
        #
        if 'virs' in tags:
            seqs_dic['___V'].append(seq)
        if 'uncla' in tags:
            seqs_dic['____U'].append(seq)
        # Populate seq lists by Domain differences between candidates and final
        if tags >= {'candi_virs','final_arkea'}:
            canVfinA_seqs.append(seq)
        if tags >= {'candi_virs','final_bacte'}:
            canVfinB_seqs.append(seq)
        if tags >= {'candi_virs','final_eukar'}:
            canVfinE_seqs.append(seq)

    # Get partial domain results
    dom_cnt = Counter()
    for domtype in DOMAIN_TPL:
        dom_cnt[domtype] = len(seqs_dic[domtype])

    # Save partial domain difference results
    tmpfile_dic = dict()
    for diftype in DIFTYPE_TPL:
        diftype_seqs = eval(diftype + '_seqs')
        if diftype_seqs:
            filename = fileprefix + '.' + diftype + lmatoutnum + '.tmp'
            count = 0
            with open(filename, 'w') as out_handle:
                writer = SeqIO.LmatIO.LmatOutWriter(out_handle)
                count = writer.write_file(diftype_seqs)
            tmpfile_dic[diftype] = (filename, count) # Tuple per each dic entry
    if verb:
        print('\t Parsed LMAT file ' + lmatoutfile)
    else:
        print('\033[90m.\033[0m', end='')
    # Return num of sequences and a dictionary with the temporary files
    #  created (if any)
    return len(seqs_lst), dom_cnt, tmpfile_dic


def plotVenn(dom_cnt, totnum, name, format):
    """Plot Venn diagram for the 4 Domains"""

    def txt(x, y, count, fs=12):
        """Aux for plotting text labels"""
        plt.text(x, y, '%.3f%%' % (count*100.0/totnum),
                 fontsize=fs, ha='center', va='center')
    
    # Create the figure
    fig = plt.figure()
    ax = fig.gca()
    #ax.grid()
    plt.axis('scaled')
    ax.set_xlim(left=-10, right=+10)
    ax.set_ylim(bottom=0, top=+15)
    ax.get_xaxis().set_ticks([])
    plt.setp(ax.get_xticklines(), visible=False)
    ax.get_yaxis().set_ticks([])
    plt.setp(ax.get_yticklines(), visible=False)
    ax.set_title("LMAT candidates multidomain analysis")
    ax.add_artist(Ellipse((0, 9), 14, 8, angle=+45, color='b', alpha=0.2))
    ax.add_artist(Ellipse((-3, 6), 14, 8, angle=-45, color='g', alpha=0.25))
    ax.add_artist(Ellipse((0, 9), 14, 8, angle=-45, color='y', alpha=0.3))
    ax.add_artist(Ellipse((+3, 6), 14, 8, angle=+45, color='r', alpha=0.25))
    plt.text(-7.5, 3, 'Archaea', fontsize=14, color='g',
             ha='center', va='center')
    plt.text(-7, 14, 'Bacteria', fontsize=14, color='y',
             ha='center', va='center')
    plt.text(7, 14, 'Eukarya', fontsize=14, color='b',
             ha='center', va='center')
    plt.text(7.5, 2.5, 'Viruses\n&viroids', fontsize=14, color='r',
             ha='center', va='center')
    plt.text(-9.8, 0.6, 'Number of tagged reads: %i' % totnum, fontsize=6)
    plt.text(-9.8, 0.2, 'Reads over domain level: %i (%.2f%%)' %
        (dom_cnt['over'], dom_cnt['over']*100.0/totnum), fontsize=6)

    # Data
    txt(-7, 9, dom_cnt['Abev'])
    txt(-3.5, 13, dom_cnt['aBev'])
    txt(3.5, 13, dom_cnt['abEv'])
    txt(7, 9, dom_cnt['abeV'])
    #
    txt(4, 5.5, dom_cnt['aBeV'])
    txt(-4, 5.5, dom_cnt['AbEv'])
    txt(0, 11, dom_cnt['aBEv'])
    txt(0, 2, dom_cnt['AbeV'])
    txt(-4.6, 11, dom_cnt['ABev'], fs=10)
    txt(4.6, 11, dom_cnt['abEV'], fs=10)
    #
    txt(2.5, 8.5, dom_cnt['aBEV'])
    txt(-2.5, 8.5, dom_cnt['ABEv'])
    txt(-1.5, 4, dom_cnt['AbEV'], fs=10)
    txt(1.5, 4, dom_cnt['ABeV'], fs=10)
    #
    txt(0, 6, dom_cnt['ABEV'])
    
    # Plot or save
    if (format == ''):
        plt.show()
    else:
        name = name + '.venn'
        if format == 'png':  # Raster image
            pngname = name + '.png'
            fig.savefig(pngname, dpi=200)
        else:  # Vectorial image
            imagename = name + '.' + format
            fig.savefig(imagename)
    # Close the plot to release memory
    plt.close()


def grow_tree(parent, path, subtree):
    """Recursive function to build the taxonomy tree"""
    if parent not in path:
        if parent in children:
            for child in children[parent]:
                subtree[child] = {}
                grow_tree(child, path + [parent], subtree[child])
    elif verb:
        print('[\033[93mWARNING\033[0m: Avoided loop for repeated taxid',
              parent, end=']')


def print_tree(subtree):
    """Recursive function to print the taxonomy tree"""
    for tid in subtree:
        if subtree[tid]:
            print(tid, end='->(')
            print_tree(subtree[tid])
        else:
            print(tid, end=',')
    print(')', end='')


def get_tree(subtree, node_set, maxdepth=0):
    """
    Recursive function to get the taxonomy subtree until some depth level
    
    maxdepth null (default) does not stop the search at any depth level.
    node_set is input/output set: at entry it should contain the root of the
    subtree, while at output it will contain all the nodes from the root to
    the maxdepth level (if any).
    """
    for tid in subtree:
        if subtree[tid]:
            if tid in node_set:
                maxdepth -= 1
                if not maxdepth:
                    break
                for child in subtree[tid]:
                    node_set.add(child)
            get_tree(subtree[tid], node_set, maxdepth)


# > MAIN
if __name__ == '__main__':
    # Argument Parser Configuration
    parser = argparse.ArgumentParser(
        description='pythonic LMAT multi-domain assignment analyzer',
        epilog=f'%(prog)s  - Release {__version__} - {__date__}'
    )
    parser.add_argument(
        '-V', '--version',
        action='version',
        version=f'%(prog)s version {__version__} released in {__date__}'
    )
    parser.add_argument(
        '-f', '--fileprefix',
        action='store',
        metavar='PREFIX',
        required=True,
        help=('file prefix of LMAT .out files')
    )
    groupMode = parser.add_mutually_exclusive_group(required=True)
    groupMode.add_argument(
        '-s', '--sequential',
        action='store_true',
        default=True,
        help=('Sequential parse the LMAT .out files')
    )
    groupMode.add_argument(
        '-p', '--parallel',
        action='store_true',
        help=('Parallel parse of the LMAT .out files')
    )
    parser.add_argument(
        '-v', '--verbose',
        action='count',
        help='increase output verbosity'
    )
    parser.add_argument(
        '-n', '--nodesfile',
        action='store',
        metavar='PATH/FILE',
        default='./nodes.dmp',
        help=('nodes file (nodes.dmp from NCBI)')
    )
    # Parser group: Mode
    groupMode = parser.add_mutually_exclusive_group(required=True)
    groupMode.add_argument(
        '-i', '--interactive',
        action='store_true',
        help=('The plots are interactive')
    )
    groupMode.add_argument(
        '-a', '--automatic',
        action='store',
        choices={'png', 'pdf', 'svg', 'ps', 'eps'},
        help=('The plots are saved with the given format')
    )

    # Parse arguments
    args = parser.parse_args()
    fileprefix = args.fileprefix
    verb = args.verbose
    nodesfile = args.nodesfile
    format = ''
    if not args.interactive:
        format = args.automatic

    # timing initialization
    start_time: float = time.time()
    # Program header
    print(f'\n=-= {sys.argv[0]} =-= v{__version__} - {__date__}'
          f' =-= by {__author__} =-=\n')
    sys.stdout.flush()

    print('\033[90mLoading NCBI nodes...\033[0m', end='')
    sys.stdout.flush()
    # build dict parents: parent for a given taxid (key)
    parents = {}
    with open(nodesfile, 'r') as f:
        for line in f:
            taxid, parent, *_ = line.split('\t|\t')
            parents[taxid] = parent

    # build dict children: dict of children for a given parent taxid (key)
    children = {}
    for tid in parents :
        if parents[tid] not in children:
            children[parents[tid]]= {}
        children[parents[tid]][tid] = 0
    print('\033[92m OK! \033[0m')

    print('\033[90mBuilding taxonomy tree and domain sets...\033[0m', end='')
    sys.stdout.flush()
    root = '1'
    tree = {root:{}}
    grow_tree(root, [], tree[root])
    #print_tree(tree)

    viroids_set = {'12884'}
    get_tree(tree, viroids_set)
    virus_set = {'10239'}
    get_tree(tree, virus_set)
    virs_set = viroids_set | virus_set

    arch_set = {'2157'}
    get_tree(tree, arch_set)

    bact_set = {'2'}
    get_tree(tree, bact_set)

    euka_set = {'2759'}
    get_tree(tree, euka_set)

    unclass_set = {'12908'}
    get_tree(tree, unclass_set)
    other_set = {'28384'}
    get_tree(tree, other_set)
    unkn_set = unclass_set | other_set

    over_set = {root, '131567'} # Over domain: root and cellular organisms
    print('\033[92m OK! \033[0m')

    lmatoutfiles = []
    #all_seqs_lst = []
    all_numseq = 0
    all_dom_cnt = Counter()
    tmpfile_dic_lst = []
    for file in os.listdir('.'):
        if (file.startswith(fileprefix) and file.endswith('.out') and
                not 'canVfinA' in file and not 'canVfinB' in file and
                not 'canVfinE' in file and not 'pyLCA' in file):
            lmatoutfiles.append(file)
    numoutfiles = len(lmatoutfiles)
    if numoutfiles:
        # Read LMAT output from the file
        print('\033[90mLoading and tagging LMAT output seqs:\033[0m ', end='')
        if verb:
            print('')
        sys.stdout.flush()
        if args.parallel:
            with mp.Pool(processes=numoutfiles) as pool:  # start worker prcs
                async_results = [pool.apply_async(
                    procLmatOut,
                    args=[lmatoutfiles[i]]
                ) for i in range(numoutfiles)]
                pool.close()
                map(mp.pool.ApplyResult.wait, async_results)
                for lmatoutfile, (numseq, dom_cnt, tmpfile_dic) in zip(
                        lmatoutfiles, [r.get() for r in async_results]):
                    all_numseq += numseq
                    all_dom_cnt += dom_cnt
                    tmpfile_dic_lst.append(tmpfile_dic)

        else:
            for lmatoutfile in lmatoutfiles:
                numseq, tmpfile_dic = procLmatOut(lmatoutfile)
                tmpfile_dic_lst.append(tmpfile_dic)
                all_numseq += numseq
        #all_seqs_lst = list(itertools.chain.from_iterable(seqs_lst_lst))
        if verb:
            print('\t', all_numseq, 'LMAT output seqs read and tagged')
        else:
            print('\033[92m OK! \033[0m')
    else:
        raise Exception('\033[91m ERROR! \033[0m No LMAT .out files found!')
    sys.stdout.flush()

    # Merge and print domain difference results
    print('\033[90mSaving extra LMAT out files: \033[0m', end='')
    if verb:
        print('')
    for diftype in DIFTYPE_TPL:
        diftype_tmp_lst = []
        count = 0
        for dic in tmpfile_dic_lst:
            if diftype in dic:
                diftype_tmp_lst.append(dic[diftype][0])
                count += dic[diftype][1]
        if diftype_tmp_lst:
            outfilename = fileprefix + '.' + diftype + '.out'
            with open(outfilename, 'wb') as outfile:
                for tmpfile in diftype_tmp_lst:
                    with open(tmpfile, 'rb') as readtmpfile:
                        shutil.copyfileobj(readtmpfile, outfile)
                    os.remove(tmpfile)
            #count = 0
            #with open(filename, 'w') as out_handle:
            #    writer = SeqIO.LmatIO.LmatOutWriter(out_handle)
            #    count = writer.write_file(diftype_seqs)
            if verb:
                print("Number of seqs labelled as V but with", diftype[-1],
                      "candidates: \033[94m%i\033[0m" % count)
                print('\t', count, 'LMAT seqs written to', outfilename)
            else:
                print('\033[90m.\033[0m', end='')
    if not verb:
        print('\033[92m OK! \033[0m')

    # Print domain results
    print("\n\033[1mResults\033[0m (KEY: upper-case = present, lower-case = excluded)")
    for domtype in DOMAIN_TPL:
        count = all_dom_cnt[domtype]
        if count or verb:
            print("Number of sequences of type %s:\033[94m %i\033[0m"
                  % (domtype, count))
    print("Total of seqs with viral taxid: %i (%.4f%%)" %
        (all_dom_cnt['___V'], all_dom_cnt['___V']/all_numseq*100))

    plotVenn(all_dom_cnt, all_numseq, fileprefix, format)

    # Timing results
    print('Total elapsed time:', time.strftime(
        "%H:%M:%S", time.gmtime(time.time() - start_time)))


# ################ OLD
"""
    if viruscall:
        print(seq.id, "labelled as virus", seq.annotations['finalcallid'])
        print(" ... but with eukariotic candidates:",
              ", ".join(str(e) for e in eukarcand))
    elif eukarcall:
        print(seq.id, "labelled as eukaria", seq.annotations['finalcallid'])
        print(" ... but with virus candidates:",
              ", ".join(str(v) for v in viruscand))
    else:
        print(seq.id, "NOT labelled as virus nor eukaria (" +
              str(seq.annotations['finalcallid']) +
              "), but with virus candidates:",
              ", ".join(str(v) for v in viruscand) +
              ", and with eukaria candidates:",
              ", ".join(str(e) for e in eukarcand))
    if viruscall:
        print(seq.id, "annotated as virus", seq.annotations['finalcallid'])
        virusnofinal = viruscand - {seq.annotations['finalcallid']}
        if virusnofinal:
            print(" ... but with other virus candidates:",
                  ", ".join(str(v) for v in virusnofinal))
    else:
        if viruscand:
            print(seq.id, "NOT annotated as virus (" +
                  str(seq.annotations['finalcallid']) +
                  "), but with virus candidates:",
                  ", ".join(str(v) for v in viruscand))

#print("".join(str(seq) for seq in vir_euk_seqs))

    print('\033[90mLoading NCBI categories...\033[0m', end='')
    sys.stdout.flush()
    with open(catfile, 'r') as f:
        for line in f:
            taxid = int(line.strip().split('\t')[2])
            if line.startswith('A'):
                arch_set.add(taxid)
            elif line.startswith('B'):
                bact_set.add(taxid)
            elif line.startswith('V'):
                virs_set.add(taxid)
            elif line.startswith('E'):
                euka_set.add(taxid)
            else:
                unkn_set.add(taxid)
    print('\033[92m OK! \033[0m')

"""

