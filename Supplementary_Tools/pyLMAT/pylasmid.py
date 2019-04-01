#!/usr/bin/env python
"""
pythonic LMAT plasmid tool

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

# pyLMAT release information
__version__ = '0.4.5'
_scriptname = 'pylasmid'
_verdata = 'Ago 2015'

# Predefined internal constants
_PATH = '.'
_NODESFILE = 'nodes.dmp'
_NAMESFILE = 'names.dmp'
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


def procLmatOut(*args, **kwargs):
    """Process a LMAT .out file (to be usually called in parallel!)."""
    output = ''
    lmatoutfile = args[0]
    lmatoutnum = (os.path.splitext(lmatoutfile)[0]).split('output')[-1]
    with open(lmatoutfile, 'rt') as lmo_handle:
        seqs_lst = list(SeqIO.parse(lmo_handle, "lmat-out"))
    targetseqs = [seq for seq in seqs_lst if (
        seq.annotations['finalcallid'] in plasmids)]
    lineage_cnt_dic = {}
    parsedseqs = []
    for seq in targetseqs:
        targetaxid = seq.annotations['finalcallid']
        candidict = seq.annotations['candidict']
        # Minimum score for LCA the final call score minus wstd*stdev
        try:
            targetscore = candidict[targetaxid]
        except KeyError:
            continue  # Candidates list type {-1:-1}, abort further processing!
        minscore = (targetscore -
                    float(wstd) * seq.annotations['statsdict']['stdev'])
        if minscore < setminscore:
            minscore = setminscore  # Avoid candidates with score under min
        if verb > 2:
            output += ('Minscore = t - w*s = ' + str(targetscore) +
                       ' - ' + wstd + ' * ' +
                       str(seq.annotations['statsdict']['stdev']) +
                       (' = %.4f' % minscore) + '\n')
        tid_lst = [id for id, score in candidict.items() if score >= minscore]
        # Build dict with taxid as keys, whose values are the list of nodes
        #   in the tree in the path from root to such a taxid
        tid_node_dic = {}
        for tid in tid_lst:
            if tid in parents:
                node_lst = []
                target_lst = [tid]
                trace_node(tree, target_lst, node_lst)
                if node_lst:
                    tid_node_dic[tid] = node_lst
                elif verb:
                    output += ('[\033[93mWARNING\033[0m: Failed trace_node' +
                               ' of tid ' + str(tid) + ']\n')
                if verb > 2:
                    output += (str(node_lst) + ' (score=' +
                               str(candidict[tid]) + ')\n')
            elif verb:
                output += ('[\033[93mWARNING\033[0m: Discarded unknown tid ' +
                           str(tid) + ']\n')

        # Traverse walk the tid_node_dic searching for the LCA
        lca = '1'
        done = False
        level = 1
        while not done:
            try:
                trav_lst = [tid_node_dic[tid][level] for tid in tid_node_dic]
            except IndexError:
                done = True
            else:
                if len(set(trav_lst)) == 1:
                    lca = trav_lst[0]
                    level += 1
                else:
                    done = True
        # 1st search for relevant alternative candidates
        if targetaxid not in lineage_cnt_dic:  # Initialize Counter of lineages
            lineage_cnt_dic[targetaxid] = Counter()
        inter_path_set = set()  # set of tuples got from path (list of nodes)
        remov_node_set = set()  # set of nodes whose path will be removed
        for tid in tid_node_dic:
            tid_node_set = set(tid_node_dic[tid])
            inter_path_set |= {tuple(tid_node_dic[tid])}  # [level:]
            remov_node_set |= (tid_node_set - {tid})
        if verb > 2:
            output += ('Remove node set: ' + str(remov_node_set) + '\n')
            output += ('Relevant lineages found: \n')
        for path in inter_path_set:
            if set(path) - remov_node_set:
                lineage_cnt_dic[targetaxid].update([path])
                if verb > 2:
                    if path[-1] == targetaxid:
                        output += ('  \033[1m' + str(path) + '\033[0m\n')
                    else:
                        output += ('  ' + str(path) + '\n')
        if verb:
            if verb > 1:
                output += ("LCA at level " + str(level) + " of " +
                           str(len(tid_lst)) + " taxids in seq " +
                           str(seq.name) + ": " + str(lca) +
                           " (call id was: " +
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
        parsedseqs.append(seq)
        if verb > 1:
            print(output, end='')
            output = ''
            sys.stdout.flush()
        elif verb == 1:
            print('\033[90m.\033[0m', end='')
            sys.stdout.flush()

    # Save partial changed seqs
    tmpfile = []
    count = 0
    if parsedseqs:
        tmpfile = (fileprefix + '.pylasmid.' + wstd + '.' + lmatoutnum +
                   '.tmp')
        with open(tmpfile, 'w') as out_handle:
            writer = SeqIO.LmatIO.LmatOutWriter(out_handle)
            count = writer.write_file(parsedseqs)
    if verb:
        output += ('\n\t Parsed ' + str(len(parsedseqs)) + ' seqs of ' +
                   str(len(seqs_lst)) + ' in ' + lmatoutfile + '\n')
    else:
        output += '\033[90m.\033[0m'
    print(output, end='')
    sys.stdout.flush()
    # Return num of sequences, the temporary file, the number of seqs in it
    #   and a dictionary of plasmid:(set of interesting lineages)
    return len(seqs_lst), tmpfile, count, lineage_cnt_dic


def lineagelst2tree(lineage_lst, targetaxid, maxoccurs):
        """From list of lineages to a Biopython.Phylo tree"""
        # Traverse walk the lineage_lst searching for the overall LCA
        _TREECOLOR = "#DDDDDD"
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
        leaves = 0
        for lineage in lineage_lst:
            # Filtering: avoid elements under half the max frequency
            #  (the targetaxid occurrence), or very low frequency,
            #  and keep total number reasonable
            try:
                counts_last_id = int(lineage[-1].split('(')[0])
            except:
                pass
            else:
                if counts_last_id < maxoccurs/2:
                    continue
                elif counts_last_id < 2:
                    continue
                elif counts_last_id < maxoccurs:
                    if leaves > 19:
                        break                    
            # Algorithm to include lineages in the tree
            cld_last = cld_lca
            for nivel in range(level, len(lineage)):
                cld_tid_lst = [cld.name for cld in cld_last]
                current_tid = lineage[nivel]
                if current_tid not in cld_tid_lst:
                    cld_new = Phylo.BaseTree.Clade(name=lineage[nivel],
                                                   clades=[])
                    try:
                        lineage_last_id = (lineage[-1].split(')')[0]
                                           ).split('(')[-1]
                        if lineage_last_id == targetaxid:
                            cld_new.color = "blue"
                        else:
                            cld_new.color = _TREECOLOR
                    except:
                        cld_new.color = _TREECOLOR
                    cld_last.clades.append(cld_new)
                    cld_last = cld_new
                    leaves += 1
                else:
                    for cld in cld_last:
                        if cld.name == current_tid:
                            cld_last = cld
                            break
        if leaves == 0:
            tree = False
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
        help='increase output verbosity (from -v to -vvv)'
    )
    parser.add_argument(
        '-n', '--nodespath',
        action='store',
        metavar='PATH',
        default='./',
        help=('path for the nodes information files (nodes.dmp and names.dmp' +
              ' from NCBI and the plasmid.names.txt LMAT file')
    )
    parser.add_argument(
        '-t', '--targetaxid',
        action='store',
        required=False,
        default=False,
        help=('target taxid if looking for a single plasmid id')
    )
    parser.add_argument(
        '-w', '--wstd',
        action='store',
        default='0.5',
        help=('weight to stdev filter (default is 0.5)')
    )
    parser.add_argument(
        '-p', '--path',
        action='store',
        default=_PATH,
        help=('relative path of the data files (if omitted, \'' +
              _PATH + '\' will be tried)')
    )
    parser.add_argument(
        '-m', '--minscore',
        action='store',
        default='1',
        help=('minimum filtering score (default is 1)')
    )

    # Parse arguments
    args = parser.parse_args()
    fileprefix = args.fileprefix
    verb = args.verbose
    nodesfile = os.path.join(args.nodespath, _NODESFILE)
    namesfile = os.path.join(args.nodespath, _NAMESFILE)
    plasmidfile = os.path.join(args.nodespath, _PLASMIDFILE)
    targetaxid = args.targetaxid
    path = args.path
    wstd = args.wstd
    setminscore = float(args.minscore)

    # Program header and chdir
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
    # build dict names: name for a given taxid (key)
    names = {}
    try:
        with open(namesfile, 'r') as f:
            for line in f:
                if 'scientific name' in line:
                    taxid, name, *_ = line.split('\t|\t')
                    names[taxid] = name
    except:
        raise Exception('\n\033[91mERROR!\033[0m Cannot read "' +
                        namesfile + '"')
    # Add also plasmid taxid and name
    plasmids = set()
    try:
        with open(plasmidfile, 'r') as f:
            count = 0
            for line in f:
                count += 1
                taxid, parent, _, _, longid, *_ = line.split('\t')
                if taxid not in parents:
                    parents[taxid] = parent
                    plasmids.add(taxid)
                    try:
                        _longid = (longid.split('|')[-1]).split('[')[0]
                        names[taxid] = (_longid).split(',')[0]
                    except:
                        print('\033[93mWARNING\033[0m: Failed parse of' +
                              'plasmid ' + longid)
                elif verb > 3:
                    print('[Discarding bad id plasmid ' + taxid + ']', end='')
    except:
        raise Exception('\n\033[91mERROR!\033[0m Cannot read "' +
                        plasmidfile + '"')
    if targetaxid:
        plasmids = {targetaxid}  # Only targeting entered plasmid

    # Build dict children: dict of children for a given parent taxid (key)
    children = {}
    for tid in parents:
        if parents[tid] not in children:
            children[parents[tid]] = {}
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
    tree = {root: {}}
    grow_tree(root, [], tree[root])
    print('\033[92m OK! \033[0m')

    # Initialize and search for valid .out files
    lmatoutfiles = []
    all_numseq = 0
    all_target = 0
    tmpfile_lst = []
    lineage_cnt_dic = {}
    for file in os.listdir('.'):
        if (file.startswith(fileprefix) and file.endswith('.out') and
                all([eval("x not in file") for x in ('canVfinA', 'canVfinB',
                    'canVfinE', 'pyLCA', 'pylasmid', 'pyhage')])):
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
            for lmatoutfile, (numseq, tmpfile, count, line_cnt_dic) in zip(
                    lmatoutfiles, [r.get() for r in async_results]):
                all_numseq += numseq
                for targetaxid, interest_path_cnt in line_cnt_dic.items():
                    if targetaxid not in lineage_cnt_dic:
                        lineage_cnt_dic[targetaxid] = Counter()
                    lineage_cnt_dic[targetaxid].update(interest_path_cnt)
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

    # Final search for relevant alternative candidates
    for targetaxid, interest_path_cnt in lineage_cnt_dic.items():
        alternative_path_cnt = Counter()
        remov_node_set = set()  # set of nodes whose path will be removed
        for path in interest_path_cnt:
            remov_node_set |= set(path[:-1])
        for path, count in interest_path_cnt.items():
            if set(path) - remov_node_set:
                alternative_path_cnt.update({path: count})
        mostcomm_lst = alternative_path_cnt.most_common()
        print('\033[90mLineage of relevant among ' + str(len(mostcomm_lst)) +
              ' candidates for plasmid ' + targetaxid +
              ':\033[0m')
        if len(mostcomm_lst) > 1:
            # Convert list of (lineage, occurrences) to list of canditates and
            #   substitute candidate taxid by occurrences:taxid:name compound
            lineage_lst = []
            for (lineage, occurs) in mostcomm_lst:
                tmp_lineage = list(lineage)
                tmp_lineage[-1] = (str(occurs) + '(' +
                                   lineage[-1] + ')' + names[lineage[-1]])
                lineage_lst.append(tmp_lineage)
            Tree = lineagelst2tree(lineage_lst, targetaxid,
                mostcomm_lst[0][1])  # Maximum absolute frequency 
            if Tree:
                mpl.rcParams.update({'font.size': 4})
                Phylo.draw(Tree, do_show=False)
                ax = plt.gca()
                ax.set_title(fileprefix + ': LMAT out LCA analysis for plasmid ' +
                             targetaxid + ' (' + names[targetaxid] + ')')
                imagename = (fileprefix + '.' + targetaxid + '.lineage.pylasmid.' +
                             wstd + '.pdf')
                plt.savefig(imagename)
                plt.close()
                try:
                    Phylo.draw_ascii(Tree)
                except:
                    pass
            else:
                print('\033[1mSorry, all the candidates found were filtered ' +
                      'for plasmid ' + targetaxid + '.\033[0m\n')
        else:
            print('\033[1mSorry, no relevant alternative candidates found ' +
                  'for plasmid ' + targetaxid + '.\033[0m\n')
        sys.stdout.flush()

    # Merge result file
    if tmpfile_lst:
        print('\033[90mSaving .pylasmid.' + wstd + '.out file: \033[0m',
              end='')
        if verb:
            print('')
        sys.stdout.flush()
        outfilename = fileprefix + '.pylasmid.' + wstd + '.out'
        with open(outfilename, 'wb') as outfile:
            for tmpfile in tmpfile_lst:
                with open(tmpfile, 'rb') as readtmpfile:
                    shutil.copyfileobj(readtmpfile, outfile)
                os.remove(tmpfile)
        if verb:
            print('\t', all_target, 'LMAT seqs written to', outfilename)
        else:
            print('\033[92m OK! \033[0m')
