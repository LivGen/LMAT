#!/usr/bin/env python

from sys import *

usage = '''
usage: %s input_dir
where: input_dir contains the unzipped files:  gi_taxid_nucl.dmp.gz,  gi_taxid_nucl.dmp.gz
''' % argv[0]

if len(argv) != 2 :
  print usage
  exit(1)

d =argv[1]

#get map: tax_id -> parent
a = open(d + '/nodes.dmp')
parents = {}
ranks = {}
for line in a :
  t = line.split('|')
  tid = t[0].strip()
  parent = t[1].strip()
  parents[tid] = parent
  rank = t[2].strip()
  ranks[tid] = rank

#get map: tax_id -> name
a = open(d + '/names.dmp')
names = {}
children = {}
for line in a :
  if line.find('scientific name') != -1 :
    t = line.split('|')
    tid = t[0].strip()
    name = t[1].strip()
    names[tid] = name
    children[tid] = {}

#build map: tax_id -> {children}
for tid in parents :
  children[parents[tid]][tid] = 0

#write taxonomy file in our preferred format
#'''
out = open('ncbi_taxonomy.dat', 'w')
out.write('#format, line 1: tid num_children list_of_children parent\n')
out.write('#format, line 2: name; third line in file is the number of entries\n')
out.write(str(len(parents)) + '\n')
for tid in parents.keys() :
  out.write(tid + ' ' + str(len(children[tid])) + ' ')
  for child in children[tid].keys() :
    out.write(child + ' ')
  out.write(parents[tid] + '\n')
  out.write(names[tid] + '\n')
out.close()
#'''

#build depth map
depths = {}
for tid in parents.keys() :
  t = tid
  depth = 0
  while parents[t] != t :
    depth += 1
    t = parents[t]
  depths[tid] = depth 

#next, build ncbi_taxonomy_rank.txt and depth file
out = open('ncbi_taxonomy_rank.txt', 'w')
out_depth = open('depth_for_ncbi_taxonomy.dat', 'w')
for tid in parents.keys() :
  out.write('depth=' + str(depths[tid]) + ',taxid=' + tid + ',ktaxid=' + tid + ',entries=-1\t')
  t = tid
  lineage = []
  while True :
    lineage.append(ranks[t]+','+names[t])
    if parents[t] == t : break
    t = parents[t]
  for j in range(len(lineage)-1, -1, -1) :
    if lineage[j].find('no rank,root') == -1 :
      out.write(lineage[j]+'\t')
  out.write('\n')
  out_depth.write(tid + ' ' + str(depths[tid])+'\n')

out.close()
out_depth.close()

