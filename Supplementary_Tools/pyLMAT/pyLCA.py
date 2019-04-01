#!/usr/bin/env python
"""
pythonic LMAT LCA tool

"""

# Python2 compability
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os
import sys
import argparse
import shutil
import multiprocessing as mp
from collections import Counter
import numpy as np
from Bio import SeqIO, Phylo
import matplotlib as mpl
mpl.use('Agg')  # Force mpl backend to avoid troubles with unset $DISPLAY    
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
from Bio import SeqIO, Phylo

# pyLMAT release information
__version__ = '0.2.2'
_scriptname = 'pyLCA'
_verdata = 'May 2015'

# Predefined internal constants
_PATH = '.'
_NODESFILE = 'nodes.dmp'
_PLASMIDFILE = 'plasmid.names.txt'


def grow_tree(parent, path, subtree):
    """Recursive function to build the taxonomy tree"""
    if parent not in path:
        if parent in children:
            for child in children[parent]:
                subtree[child] = {}
                grow_tree(child, path + [parent], subtree[child])
    elif verb > 2:
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


def trace_node(subtree, target_lst, node_lst):
    """Recursive function to build a list of nodes from root to target taxid"""
    for tid in subtree:
        if subtree[tid]:
            node_lst.append(tid)
            if target_lst[0] in subtree[tid]:
                target = target_lst.pop()
                if target != root:  # Avoid to append twice root node if target
                    node_lst.append(target)
                break
            trace_node(subtree[tid], target_lst, node_lst)
            if not target_lst:
                break
            else:
                node_lst.pop()
    #print(node_lst)


def procLmatOut(*args, **kwargs):
    """Process a LMAT .out file (to be usually called in parallel!)."""
    output = ''
    lmatoutfile = args[0]
    lmatoutnum = (os.path.splitext(lmatoutfile)[0]).split('output')[-1]
    with open(lmatoutfile, 'rU') as lmo_handle:
        seqs_lst = list(SeqIO.parse(lmo_handle, "lmat-out"))
    targetseqs = [seq for seq in seqs_lst if (
        seq.annotations['finalcallid'] == targetaxid)]
    interest_path_set = set()
    for seq in targetseqs:
        candidict = seq.annotations['candidict']
        # Minimum score for LCA the final call score minus wstd*stdev
        try:
            targetscore = candidict[targetaxid]
        except KeyError:
            break  # List of candidates type {-1:-1}, abort further processing!
        minscore = targetscore - wstd * seq.annotations['statsdict']['stdev']
        if minscore < 1:
            minscore = 1  # Jonathan: avoid candidates with score under 1
        if verb > 2:
            output += ('Minscore = t - w*s = ' + str(targetscore) +
                ' - ' + str(wstd) + ' * ' +
                str(seq.annotations['statsdict']['stdev']) + (' = %.4f' %
                minscore) + '\n')
        tid_lst = [id for id, score in candidict.items() if score >= minscore]
        # Build dict with taxid as keys, whose values are the list of nodes
        #   in the tree in the path from root to such a taxid
        tid_node_dict = {}
        for tid in tid_lst:
            if tid in parents:
                node_lst = []
                target_lst = [tid]
                trace_node(tree, target_lst, node_lst)
                tid_node_dict[tid] = node_lst
                if verb > 2:
                    output += (str(node_lst) + ' (score=' +
                               str(candidict[tid]) + ')\n')
            elif verb:
                output += ('[\033[93mWARNING\033[0m: Discarded unknown tid ' +
                           str(tid) + ']\n')

        # Traverse walk the tid_node_dict searching for the LCA
        lca = '1'
        done = False
        level = 1
        while not done:
            try:
                trav_lst = [tid_node_dict[tid][level] for tid in tid_node_dict]
            except IndexError:
                done = True
            else:
                if len(set(trav_lst)) == 1:
                    lca = trav_lst[0]
                    level += 1
                else:
                    done = True
        # 1st search for interesting candidates (not in the final call path)
        inter_path_set = set()  # set of tuples got from path (list of nodes)
        remov_node_set = set()  # set of nodes whose path will be removed
        for tid in tid_node_dict:
            tid_node_set = set(tid_node_dict[tid])
            inter_path_set |= {tuple(tid_node_dict[tid])} #[level:]
            remov_node_set |= tid_node_set - {tid}
        if verb > 2:
            output += ('Remove node set: ' + str(remov_node_set) + '\n')
            output += ('Interesting lineages found: \n')
        for path in inter_path_set:
            if set(path) - remov_node_set:
                interest_path_set |= {path}
                if verb > 2:
                    if path[-1] == targetaxid:
                        output += ('  \033[1m' + str(path) + '\033[0m\n')
                    else:
                        output += ('  ' + str(path) + '\n')
        if verb:
            if verb > 1:
                output += ("LCA at level " + str(level) + " of " +
                    str(len(tid_lst)) + " taxids in seq " + str(seq.name) +
                    ": " + str(lca) + " (call id was: " +
                    str(seq.annotations['finalcallid']) + ")\n")
        # Update the final call info
        seq.annotations['finalcallid'] == lca
        finalcall = seq.annotations['finalcall'].split()
        finalcall[0] = lca
        if lca in candidict:
            finalcall[1] = str(candidict[lca])
            finalcall[2] = "LCA"
        else:
            finalcall[1] = str(minscore)
            finalcall[2] = "LCA-HIGH"
        seq.annotations['finalcall'] = ' '.join(finalcall)

    # Save partial changed seqs
    tmpfile = []
    count = 0
    if targetseqs:
        tmpfile = fileprefix + '.pyLCA.' + lmatoutnum + '.tmp'
        with open(tmpfile, 'w') as out_handle:
            writer = SeqIO.LmatIO.LmatOutWriter(out_handle)
            count = writer.write_file(targetseqs)
    if verb:
        output += ('\t Parsed LMAT file ' + lmatoutfile + '\n')
    else:
        output += '\033[90m.\033[0m'
    print(output, end='')
    sys.stdout.flush()
    # Return num of sequences, the temporary file and the number of seqs in it
    return len(seqs_lst), tmpfile, count, interest_path_set


def lineagelst2tree(lineage_lst, targetaxid):
        """From list of lineages to a Biopython.Phylo tree"""
        # Traverse walk the lineage_lst searching for the overall LCA
        _TREECOLOR = "gray"
        lca = '1'
        done = False
        level = 1
        while not done:
            try:
                trav_lst = [lineage[level] for lineage in lineage_lst]
            except IndexError:
                done = True
            else:
                if len(set(trav_lst)) == 1:
                    lca = trav_lst[0]
                    level += 1
                else:
                    done = True
        # Build the Phylo.Tree from the LCA
        cld_lca = Phylo.BaseTree.Clade(name=lca, clades=[])
        tree = Phylo.BaseTree.Tree(rooted=True, root=cld_lca)
        tree.root.color = _TREECOLOR
        done = False
        for lineage in lineage_lst:
            cld_last = cld_lca
            for nivel in range(level, len(lineage)):
                cld_tid_lst = [cld.name for cld in cld_last]
                current_tid = lineage[nivel]
                if current_tid not in cld_tid_lst:
                    cld_new = Phylo.BaseTree.Clade(name=lineage[nivel], clades=[])
                    if lineage[-1] == targetaxid:
                        cld_new.color = "blue"
                    else:
                        cld_new.color = _TREECOLOR
                    cld_last.clades.append(cld_new)
                    cld_last = cld_new
                else:
                    for cld in cld_last:
                        if cld.name == current_tid:
                            cld_last = cld
                            break
        return(tree)

# > MAIN
if __name__ == '__main__':
    # Argument Parser Configuration
    parser = argparse.ArgumentParser(
        description='pythonic phague tool',
        epilog=_scriptname + ' - LLNL - May 2015 - ' + _verdata
    )
    parser.add_argument(
        '-V', '--version',
        action='version',
        version=_scriptname + ' release ' + __version__ + ' - ' + _verdata
    )
    parser.add_argument(
        '-f', '--fileprefix',
        action='store',
        metavar='PREFIX',
        required=True,
        help=('file prefix of LMAT .out files')
    )
    parser.add_argument(
        '-v', '--verbose',
        action='count',
        default=0,
        help='increase output verbosity'
    )
    parser.add_argument(
        '-n', '--nodespath',
        action='store',
        metavar='PATH',
        default='./',
        help=('path for the nodes information files (nodes.dmp from NCBI ' +
              'and, optional, the plasmid.names.txt LMAT file')
    )
    parser.add_argument(
        '-t', '--targetaxid',
        action='store',
        required=True,
        help=('target taxid')
    )
    parser.add_argument(
        '-w', '--wstd',
        action='store',
        default=0.5,
        help=('weight to stdev filter (default is 0.5)')
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
    fileprefix = args.fileprefix
    verb = args.verbose
    nodesfile = os.path.join(args.nodespath, _NODESFILE)
    plasmidfile = os.path.join(args.nodespath, _PLASMIDFILE)
    targetaxid = args.targetaxid
    path = args.path
    wstd = float(args.wstd)

    # Program Header and chdir
    print('\n=-= ' + _scriptname + ' =-= v' + __version__ + ' =-= ' +
          _verdata + ' =-= LLNL =-=\n')

    # Load NCBI and plasmid nodes
    print('\033[90mLoading NCBI and plasmid nodes...\033[0m', end='')
    sys.stdout.flush()
    # build dict parents: parent for a given taxid (key)
    parents = {}
    try:
        with open(nodesfile, 'r') as f:
            for line in f:
                taxid, parent, *_ = line.split('\t|\t')
                parents[taxid] = parent
    except:
        raise Exception('\n\033[91mERROR!\033[0m Cannot read "' +
                        nodesfile + '"')
    # Add also plasmid taxid (optional, if file found)
    try:
        with open(plasmidfile, 'r') as f:
            for line in f:
                taxid, parent, *_ = line.split('\t')
                parents[taxid] = parent
    except:
        print('\033[93mWARNING\033[0m: Cannot read "' +
                    plasmidfile + '". Plasmid taxid not loaded!')
    # Build dict children: dict of children for a given parent taxid (key)
    children = {}
    for tid in parents :
        if parents[tid] not in children:
            children[parents[tid]]= {}
        children[parents[tid]][tid] = 0
    print('\033[92m OK! \033[0m')
    # Change directory to LMAT output files path
    try:
        os.chdir(path)
    except (OSError):
        raise Exception('\033[91mERROR!\033[0m Cannot change to directory "' +
                        path + '"')
    else:
        if verb and path is not _PATH:
            print('\033[90mNOTE: Base working directory changed to "' + path +
                  '"\033[0m')

    # Build taxonomy tree
    print('\033[90mBuilding taxonomy tree...\033[0m', end='')
    sys.stdout.flush()
    root = '1'
    tree = {root:{}}
    grow_tree(root, [], tree[root])
    ## print_tree(tree)
    print('\033[92m OK! \033[0m')

    lmatoutfiles = []
    all_numseq = 0
    all_target = 0
    tmpfile_lst = []
    interest_path_set = set()
    for file in os.listdir('.'):
        if (file.startswith(fileprefix) and file.endswith('.out') and
                not 'canVfinA' in file and not 'canVfinB' in file and
                not 'canVfinE' in file and not 'pyLCA' in file):
            lmatoutfiles.append(file)
    numoutfiles = len(lmatoutfiles)
    if numoutfiles:
        # Read LMAT output from the file
        print('\033[90mProcessing LMAT output seqs:\033[0m ', end='')
        if verb:
            print('')
        sys.stdout.flush()
        with mp.Pool(processes=numoutfiles) as pool:  # start worker prcs
            async_results = [pool.apply_async(
                procLmatOut,
                args=[lmatoutfiles[i]]
            ) for i in range(numoutfiles)]
            pool.close()
            map(mp.pool.ApplyResult.wait, async_results)
            for lmatoutfile, (numseq, tmpfile, count, inter_path_set) in zip(
                    lmatoutfiles, [r.get() for r in async_results]):
                all_numseq += numseq
                interest_path_set |= inter_path_set
                if count:
                    all_target += count
                    tmpfile_lst.append(tmpfile)
        if verb:
            print('\t', all_numseq, 'LMAT output seqs read')
        else:
            print('\033[92m OK! \033[0m')
    else:
        raise Exception('\033[91m ERROR! \033[0m No LMAT .out files found!')
    sys.stdout.flush()

    # Final search for interesting candidates (not in the final call path)
    alternative_path_set = set();
    remov_node_set = set()  # set of nodes whose path will be removed
    for path in interest_path_set:
        remov_node_set |= set(path[:-1])
    for path in interest_path_set:
        if set(path) - remov_node_set:
            alternative_path_set |= {path}
    lineage_lst = list(alternative_path_set)
    lineage_lst.sort()
    if len(lineage_lst) > 1:
        print('\033[90mLineage of', len(lineage_lst),
              'interesting candidates:\033[0m')
        for path in lineage_lst:
            if path[-1] == targetaxid:
                print('  \033[94m\033[1m' + str(path) + '\033[0m')
            else:
                print('  ' + str(path))
    else:
        print('\033[1mSorry, no alternative candidates found!\033[0m')

    tree = lineagelst2tree(lineage_lst, targetaxid)
    mpl.rcParams.update({'font.size': 3})
    Phylo.draw(tree, do_show=False)
    ax = plt.gca()
    ax.set_title(fileprefix + ": LMAT output LCA analysis for taxid " +
                 targetaxid)
    imagename = fileprefix + '.lineage.pyLCA.' + args.wstd + '.pdf'
    plt.savefig(imagename)
    plt.close()
    try:
        Phylo.draw_ascii(tree)
    except:
        pass

    # Merge result pyLCA files
    if tmpfile_lst:
        print('\033[90mSaving .pyLCA.' + args.wstd +
              '.out file: \033[0m', end='')
        if verb:
            print('')
        sys.stdout.flush()
        outfilename = fileprefix + '.pyLCA.' + args.wstd + '.out'
        with open(outfilename, 'wb') as outfile:
            for tmpfile in tmpfile_lst:
                with open(tmpfile, 'rb') as readtmpfile:
                    shutil.copyfileobj(readtmpfile, outfile)
                os.remove(tmpfile)
        if verb:
            print('\t', all_target, 'LMAT seqs written to', outfilename)
        else:
            print('\033[92m OK! \033[0m')
