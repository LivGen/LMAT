#!/usr/bin/env python

import sys
import pprint

newstartid=100000000
#read in remapping data
mapfile=sys.argv[1]

#a = open('/usr/mic/post1/metagenomics/ref_sets/fasta/20130206update/microbe3.20130207.headers.mapping.ncbionly.multi')
a = open(mapfile,'r')

### note that the first time around the altered tax ids will be identified on this line 
###      print 'failed to find name entry for tid=', tid
### after this, you have to go to the web site to see what the new updated tax ids are
manual_conv={}
if len(sys.argv) >= 3 :
   u = open(sys.argv[2])
   for lin in u :
      lin=lin.rstrip()
      vals=lin.split()
      manual_conv.setdefault(vals[0],vals[1])

##
print("Manually identified the following recent changes in NCBI taxonomy",manual_conv)

#build map: old_tid -> [list of new (working) tids]
ntids={}
tid_old_to_new = {}
save_hdr={}
counter=0
for line in a :
 # counter+=1
  #print(str(counter)+"\n")
  t = line.rstrip().split('\t')
  tid_old =t[1]
  if tid_old in manual_conv.keys() :
     tid_old = manual_conv[t[1]]
  tid_new = t[0]
  ntids.setdefault(tid_new,1)
  if not tid_old in tid_old_to_new.keys() :
     tid_old_to_new[tid_old] = []
  tid_old_to_new[tid_old].append(tid_new)
  save_hdr.setdefault(tid_new,t[4])

#remap: kpath_taxonomy.dat
a = open('taxonomy/ncbi_taxonomy.dat')
a.readline()
a.readline()
a.readline()

out = open('taxonomy/ncbi_taxonomy.segment.dat', 'w')
out.write('#format: ID num_children child_1 ... parent\n')
out.write('#best_name_we_can_give_it\n')
out.write('0\n')

#read in existing taxonomy
children = {}
names = {}
parent = {}
while True :
  f = a.readline()
  if len(f) == 0 : break
  name = a.readline()
  assert(len(name) > 0)
  t = f.split()
  tid = t[0]
  parent[tid] = t[-1]
  names[tid] = name[:-1]
  children[tid] = {}
  for x in t[2:-1] :
    children[tid][x] = 0

#insert leaves for plasmids
for tid in tid_old_to_new.keys() :
  ntid_lst = tid_old_to_new[tid]
  assert ntid_lst != []
  for ntid in ntid_lst :
     if not ntid in ntid_lst:
       print('failed to find name entry for tid=', tid)
       continue
     if ntid != tid :
        parent[ntid] = tid
        names[ntid] = names[tid]
        children[tid][ntid] = 0
        children[ntid] = {}

for tid in parent.keys() :
  tmp = str(len(children[tid].keys()))
  ntids[tid] = 2
  out.write(tid+' ' + tmp+' ')
  for child in children[tid].keys() :
    out.write(child + ' ')
  out.write(parent[tid]+'\n')  
  out.write(names[tid]+'\n')
out.close()

for chk in ntids.keys() :
   if ntids[chk] != 2 :
      print("weird: ",ntids[chk], chk)


node2depth={}
dfile = open('taxonomy/depth_for_ncbi_taxonomy.dat')
for line in dfile :
   line=line.rstrip()
   vals=line.split()
   node2depth[vals[0]]=vals[1]

#remap: ncbi_taxonomy_rank.txt
a = open('taxonomy/ncbi_taxonomy_rank.txt').readlines()
print('>>>>>>>>>>>>>>>>>', len(a))
out = open('ncbi_taxonomy_rank.segment.txt','w')
out3 = open('depth_for_ncbi_taxonomy.segment.dat','w')
found = {}
xx = 0
for line in a :
  #out2.write(line)
  line=line.rstrip()
  xx += 1
  t = line.split(',')
  tt = line.split('\t')
  t2 = t[2].split('=')
  old_tid = t2[1]
  depth = node2depth[old_tid]
  ndepth = int(depth)+1
  out3.write(old_tid+" "+str(depth) + "\n")
  assert old_tid != '-1'
  if not old_tid in tid_old_to_new:
    out.write(line + "\n")
    continue
  ntid_lst = tid_old_to_new[old_tid]
  out.write(line + "\n")
  assert ntid_lst[0]
  for ntid in ntid_lst : 
     if ntid != old_tid :
       assert ntid != '-1'
       out3.write(str(ntid)+" "+str(ndepth) + "\n")
       #print "for read? ", len(t), "update with",ntid
       for it in range(4) :
         if( it == 1 ) :
            out.write("taxid="+str(ntid)+",")
         elif( it == 2 ) :
            out.write("ktaxid="+str(ntid)+",")
         elif (it < 3 ) : 
            out.write(t[it]+",")
         elif (it == 3 ) : 
            out.write(t[it])
       pidx1 = save_hdr[ntid].find('Plasmid')
       pidx2 = save_hdr[ntid].find('plasmid')
       pstr = ''
       if pidx1 != -1 :
          lidx=save_hdr[ntid].find(',',pidx1+1)
          pstr = save_hdr[ntid][pidx1:lidx]
       elif pidx2 != -1 :
          lidx=save_hdr[ntid].find(',',pidx2+1)
          pstr = save_hdr[ntid][pidx2:lidx]
      
       for it1 in range(1,len(tt),1) :
          out.write("\t"+tt[it1]) 
          if it1 == (len(tt)-1) and  pstr != '' :
            out.write(" "+pstr) 
       out.write("\n")
#print 'num new plasmid ID nodes:', num    
print('num lines read:', xx)
out.close()
#out2.close()


print('len(tid_old_to_new):', len(tid_old_to_new))


#for tid in tid_old_to_new.keys() :
  #if not found.has_key(tid) :
    #print 'no found, tid:', tid, ' gids:',
    #for t in tid_to_gid[tid].keys() : print t,
    #print
