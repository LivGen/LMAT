#!/usr/bin/env python

import sys,os
import pprint

fsfile=sys.argv[1] # fastsummary file
rank_lst_str=sys.argv[2] #ranks of interest
odir=sys.argv[3] #output directory
rdir=os.getenv("LMAT_DIR")
taxtree=rdir + "/" + "ncbi_taxonomy.segment.pruned.dat.nohl"
rankfile=rdir + "/" + "ncbi_taxid_to_rank.pruned.txt"
plasfile=rdir + "/" + "low_numid_plasmids.txt"
plasnf=rdir + "/" + "plasmid.names.txt"
plasmids={}

if os.path.isfile(plasfile) :
   a = open(plasfile)
   for lin in a :
      lin = lin.rstrip()
      plasmids[lin]=1

   a = open(plasnf)
   plasname={}
   for lin in a :
      lin = lin.rstrip()
      vals=lin.split('\t')
      plasname[vals[0]]=vals[len(vals)-1]
   
   
def isPlasmid(plasmids,tid) :
   res=False
   if plasmids.has_key(tid) or (int(tid) >= 10000000 and int(tid) < 20000000) :
      res=True
   return res


gsfile=""
min_gene_cnt=2
if len(sys.argv) >= 5 :
   gsfile=sys.argv[4]
   min_gene_cnt=int(sys.argv[5])

ranktable={}

a = open(rankfile)
for lin in a :
   vals=lin.split()
   ranktable[vals[0]] = vals[1]

a = open(taxtree)
a.readline()
a.readline()
a.readline()

#read in existing taxonomy
children = {}
names = {}
parent = {}
parent.setdefault(1,1)
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

rank_lst=rank_lst_str.split(',')

def getRankTid(rank,tid,ranks,parents,plasmids) :
   stid=tid
   rank_tid=-1

   if (ranks.has_key(stid) and rank == ranks[stid]) or (rank == "plasmid" and isPlasmid(plasmids,tid))  :
      rank_tid=tid
   else :
      while parent[stid] != stid :
         if  ranks.has_key(stid) and ranks[stid] == rank :
            rank_tid=stid
            break
         stid = parent[stid]
   return rank_tid
         
orig={}
store={}
a = open(fsfile)
for line in a :
   line = line.rstrip()
   t = line.split('\t')
   wrc = t[0]
   count = t[1]
   taxid = t[2]
   descrip=''
   orig[taxid]=t[3]
   if not parent.has_key(taxid) :
      if taxid != 1 :
         #print 'warning: did not find parent id for node (ignore)', taxid,'for entry:'
         #taxid=-1
         #print line
         parent.setdefault(taxid,1)
   for rank in rank_lst :
      tid=getRankTid(rank,taxid,ranktable,parent,plasmids)
      if tid == -1 :
         continue
      if not store.has_key(rank) :
         store[rank]={}
      if not store[rank].has_key(tid) :
         store[rank][tid]=[]
      store[rank][tid].append( (taxid,wrc,count) )
gene_store={}
gene_cnt={}
if gsfile != "" : 
   a = open(gsfile)
   for line in a :
      line = line.rstrip()
      t = line.split('\t')
      rc = t[1]
      taxid = t[2]
      if taxid == '0' :
         ## means this read was not assigned to a taxid
         continue
      geneid = t[4]
      type=t[7]

      if not parent.has_key(taxid) :
         #print 'warning: no parent node for', taxid,' entry: (rRNA gene not counted)'
         #print line
         #taxid=-1
         #continue
         parent.setdefault(taxid,1)

      for rank in rank_lst :
         tid=getRankTid(rank,taxid,ranktable,parent,plasmids)
         ### track how many reads were mapped to rRNA
         if tid == -1 :
            continue
         if type == "rRNA" :
            if not gene_store.has_key(rank) :
               gene_store[rank]={}
            if not gene_store[rank].has_key(tid) :
               gene_store[rank][tid]=[]
            gene_store[rank][tid].append( (taxid,rc) )

         if int(rc) > min_gene_cnt :
            if not gene_cnt.has_key(rank) :
               gene_cnt[rank]={}
            if not gene_cnt[rank].has_key(tid) :
               gene_cnt[rank][tid]={}
            if not gene_cnt[rank][tid].has_key(geneid) :
               gene_cnt[rank][tid][geneid] = 0
            gene_cnt[rank][tid][geneid] += int(rc)


for rank in store.keys() :
   fsname=os.path.basename(fsfile)
   fsfileout=odir+"/"+fsname+"."+rank
   print "create fastsummary file ",fsfileout, "for rank=",rank 
   fh=open(fsfileout,"w")
   save=[]
   for tid in store[rank].keys() :
      ## these plasmids do not have descriptive names so use the original headers
      if plasmids.has_key(tid) and plasname.has_key(tid) and rank == "plasmid" :
         name_str=plasname[tid]
      elif orig.has_key(tid) :
         name_str=orig[tid]
      else :
         name_str=names[tid]
      idx =  name_str.find(',') 
      if idx != -1 :
         name_str=name_str[idx+1:]
      ### taxid sum
      lst=store[rank][tid]
      ## where plasmids are located at the species node, only report them
      ## in species output when, there is a descendant
      if len(lst) == 1 and isPlasmid(plasmids,tid) and rank != "plasmid" :
         continue
      best_wrc,best_count=-1,-1
      top_strain=-1
      wrc_sum,count_sum=0,0
      for taxid,wrc,count in lst :
         if isPlasmid(plasmids,taxid) :
            ranktable[taxid] = "plasmid"
         wrc_sum += float(wrc)
         count_sum += int(count)
         if rank == "species" and ranktable[taxid] == "strain" :
            if best_wrc < float(wrc) :
               top_strain=taxid
               best_wrc = float(wrc)
               best_count = count
      strain_info=""
      if top_strain != -1 :
         strain_info = "\t"+str(best_wrc) + "\t"+str(best_count)+"\t"+top_strain+"\t"+orig[top_strain]
      
      gene_lst = []
      rrna_csum=0
      if gene_store.has_key(rank) and gene_store[rank].has_key(tid) :
         gene_lst=gene_store[rank][tid]
         for taxid,count in gene_lst :
            rrna_csum += int(count)
      gene_id_lst=[]
      gene_read_cnt=0
      if gene_cnt.has_key(rank) and gene_cnt[rank].has_key(tid) :
         gene_id_lst=gene_cnt[rank][tid].keys()
         for gid in gene_id_lst :
            gene_read_cnt += gene_cnt[rank][tid][gid]

      tup=(wrc_sum,count_sum,tid,name_str,rrna_csum,len(gene_id_lst),gene_read_cnt,strain_info)
      save.append(tup)
   sval=sorted(save, key=lambda val : val[0],reverse=True)
   if gsfile != "" :
      key_str="Average Read Score\tTotal Read Score\tRead Count\tPcnt. rRNA\tNo. Genes\tNo. Gene Reads\tTaxID\tName\tStrain Info"
   else :
      key_str="Average Read Score\tTotal Read Score\tRead Count\tTaxID\tName\tStrain Info"
   
   fh.write(key_str +"\n")
   for val in sval :
      avg = float(val[0])/float(val[1])
      astr="%.4f" % avg
      if gsfile != "" :
         pcnt = float(val[4])/float(val[1])
         fstr ="%.4f" % pcnt
         out_str=astr+"\t"+str(val[0])+"\t"+str(val[1]) + "\t" + fstr + "\t" +str(val[5])+"\t"+str(val[6])+ "\t" + str(val[2])+"\t"+val[3] + val[7]
      else :
         out_str=astr+"\t"+str(val[0])+"\t"+str(val[1]) + "\t" + str(val[2])+"\t"+val[3] + val[7]
      fh.write(out_str +"\n")
