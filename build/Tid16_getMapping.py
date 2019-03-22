#!/usr/bin/python

from sys import *

usage = '''
usage: %s tid_fn taxonomy_fn output_fn

where: 'tid_fn' is the output from TID16_getMapping.py

note: 'taxonomy_fn' should be the "complete" taxonomy

''' % (argv[0])

if len(argv) != 4 :
  print usage
  exit(9)

tid_32 = {}
a = open(argv[1])
for line in a :
  assert(not tid_32.has_key(line[:-1]))
  tid_32[line[:-1]] = 0
print 'tid leaf count=', len(tid_32)


#read in the taxonomy
a = open(argv[2])
#discard header lines
a.readline()
a.readline()
a.readline()

parents = {}
names = {}
children = {}

while True :
  line = a.readline()
  if len(line) < 2 : break
  name = a.readline()
  assert(len(name))
  t = line.split()
  tid = t[0]
  parent = t[-1]
  parents[tid] = parent
  names[tid] = name
  children[tid] = {}
  for child in t[2:-1] :
    children[tid][child] = 0

#collect all 32-bit tax ids; this gets the nodes that
#are between the leaves (the IDs from the fasta file)
#and the root.
tid_32_needed = {}
for tid in tid_32.keys() :
  if not tid in parents:
    print(str(tid)+"\n")
    continue
  parent = tid
  while True :
    if (parent > 1):
      tid_32_needed[parent] = 0

    next = parents[parent]

    #exit condition: we've reached the root
    if next == parent : 
      assert(next == "1")
      break

    parent = next
print 'number of IDs in taxonomy subtree:', len(tid_32_needed)

#form 32bit -> 16bit mapping
mp = {}
id = 2

mp["1"] = 1

for tid in tid_32_needed.keys() :

  if not tid == "1":
    mp[tid] = id
    id += 1

#write mapping file
out = open(argv[3],'w')
for id in mp.keys() :
  out.write(str(id) + ' ' + str(mp[id]) + '\n')
out.close()
