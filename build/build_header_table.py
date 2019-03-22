#!/usr/bin/python

from sys import *
import gzip

usage = '''
usage: %s  fasta.fa gi_to_tid_map output_directory

  transforms fasta file containing reference genome into format with a
  single id for header; the header is the NCBI taxonomy node ID associated
  with the sequences gi number. If the input file is 'fasta.fa,'
  the output files will be:

    fasta.fa.int        - same as fasta.fa, but with transformed headers
    fasta.fa.gi.table   - each two-line entry contains (1) gi number
                          (2) header from the input file
    fasta.tid.table     - each two-line entry contains (1) tax node ID;
                          (2) header from the input file

''' % argv[0]

if len(argv) != 4 :
  print usage
  exit(1)

#read gi to tid mapping
print 'reading:', argv[2]
gi_to_tid = {}
a = open(argv[2])
for line in a :
  line=line.rstrip()
  t = line.split('\t')
  
  gi_to_tid[t[4]] = t[0]

fasta = argv[1]
output_dir = argv[3]
j = fasta.rfind('/')
j += 1
base = argv[1][j:]


#open output files
print 'opening for write:', output_dir + '/' + base + '.int'
out_seq = open(output_dir + '/' + base + '.int', 'w')

print 'opening for write:', output_dir + '/' + base + '.gi.table'
out_gi = open(output_dir + '/' + base + '.gi.table', 'w')

print 'opening for write:', output_dir + '/' + base + '.tid.table'
out_tid = open(output_dir + '/' + base + '.tid.table', 'w')

print 'reading:', argv[1]
a = open(argv[1])
ctr = 0
discard = False
first = True
for line in a :
  ctr += 1
  #parse the header line
  if line[0] == '>' :
    if not first and not discard :
      out_seq.write('\n')
    discard = False
    line = line[:-1] 
    first = False
    #t = line.split('|')
    gi = line
    #warn if we failed to find the NCBI gi
    if gi == -1 :
      was_found = False
      print 'discarding sequence since gi was not found in this header:'
      discard = True
      print line
      print 'at line',ctr
      continue

    #map the gi to tid
    if not gi_to_tid.has_key(gi)  or gi_to_tid[gi] == '-1' :
      print 'discarding sequence; failed to find mapping in gi_to_tid for gi:', gi
      print 'at line',ctr
      discard = True
      continue

    tid = gi_to_tid[gi];
    si = tid

    #warn if taxonomy doesn't contain the id
    '''
    if not tids.has_key(tid) :
      print 'discarding sequence that maps to tid', tid, 'since tid was not found in the taxonomy'
      print line
      print 'at line',ctr
      discard = True
      continue
    '''

    #write the header with the fake id
    out_seq.write('>' +si + '\n')
    out_gi.write(gi + '\n')
    out_gi.write(line + '\n')
    out_tid.write(tid + '\n')
    out_tid.write(line + '\n')

  #write the next part of the sequence
  else :  
    if not discard :
      line = line.strip()
      out_seq.write(line)

out_seq.write('\n')

out_seq.close()
out_gi.close()
out_tid.close()
