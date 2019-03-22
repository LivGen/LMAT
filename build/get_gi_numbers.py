#!/usr/bin/python

from sys import *
import gzip

usage = '''
usage: %s fasta.fa output_directory
outputs: fasta.fa.gi, which contains one gi number per line
''' % argv[0]

if len(argv) != 3 :
  print usage
  exit(1)


#get gi numbers for all sequences
b = argv[1].rfind('/')
name = argv[1] + '.gi'
if b != -1 :
  name = argv[1][b+1:] + '.gi'
out = open(argv[2] + '/' + name, 'w')
j = 0
a = open(argv[1])
for line in a :
  j += 1
  if line[0] == '>' :
    line = line[:-1] 
    t = line.split('|')
    gi = -1
    #attempt to find the gi
    for k in range(len(t)) :
      if t[k] == 'gi' :
        gi = t[k+1]
        break

    if gi == -1 :
      x = line.find('>gi|')
      if x != -1 :
        j1 = line.find('|', x+5)
        assert(j1 != -1)
        gi = line[x+4:j1]
        j = gi.find(':')
        if j != -1 :
          gi = gi[:j]

    #warn if we failed to find the gi
    if gi == -1 :
      was_found = False
      print 'discarding sequence since gi was not found in this header:'
      print line
      continue

    out.write(gi + '\n')
