#!/usr/bin/env python

from sys import *

usage = '''
usage: %s path_to_ncbi_taxonomy_rank.txt input_fn output_fn num include_only
where: num          - drop entries if the occur less frequently than num
       include_only - drop entries that do not have 'include_only' in their lineage. Example,
                      use 'Bacteria' if you just want to see bacteria. Use 'all' to disable
       
''' % argv[0]

if len(argv) != 6 :
  print usage
  exit(1)

out = open(argv[3], 'w')
num = int(argv[4])
min_avg=float(argv[5])
include_only='all'
#ranks = 'superkingdom','family','genus','species',
#order subgenus subspecies

a = open(argv[1])
tax = {}
for line in a :
  t = line.split(',')
  t = t[2]
  t = t.split('=')
  ktaxid = t[1]
  tax[ktaxid] = line

a = open(argv[2])
for line in a :
  t = line.split()
  count = t[1]
  avg=float(t[0])/float(t[1])
  ktaxid = t[2]
  descrip=''
  if not tax.has_key(ktaxid) :
    ktaxid=-1 
    print 'error: failed to find ktaxid', ktaxid,'for entry:'
    print line
    continue

  if int(ktaxid) == 1 :
    e2 = ['Root,Root\n']
  elif ktaxid != -1 :
    e = tax[ktaxid]
    j = e.find('\t')
    if j == -1 :
       print "set to root?",line,ktaxid
       e2 = ['Root,Root\n']
    else :
       e = e[j+1:]
       e2 = e.split('\t')
  else :
    e2 = descrip  
     #e2 = ['unknown,unkown\n']
  if include_only == 'all' or e.find(include_only) != -1 :
     if int(count) > num and avg >= min_avg :
       out.write(count + '\t')
       for x in e2[:-1] :
         if x.find('no rank') == -1 :
           x2 = x.split(',')
           out.write(x2[1] + '\t')
         
       j = e2[-1].find(',')
       out.write(e2[-1][j+1:])
