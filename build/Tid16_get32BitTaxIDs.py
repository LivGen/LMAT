#!/usr/bin/python

from sys import *

usage = '''
usage: %s input_fasta_filename.fa output_fn

assumes: each sequence in 'fasta_input.fa'is contained on a single line

output file will contain all tax IDs from the headers in input_fasta_filename.fa,
one tax ID per line
''' % (argv[0])

if len(argv) != 3 :
  print usage
  exit(9)

#collect tax IDs from the fasta file, or read tax IDs from intermediate file
tid_32 = {}
a = open(argv[1])
while True :
    header = a.readline()
    if len(header) < 2 : break
    assert(header[0] == '>') #weak attempt to detect 
                             #if previous sequence spanned multiple lines
    seq = a.readline()
    assert(len(seq))
    tid_32[header[1:-1]] = 0

print 'tid count:', len(tid_32)
print 'writing tids to file:' + argv[2]
out = open(argv[2], 'w')
for tid in tid_32.keys() :
  out.write(tid + '\n')
out.close()

