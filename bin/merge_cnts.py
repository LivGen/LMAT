#!/usr/bin/env python

# Python2 compatibility
from __future__ import print_function, unicode_literals

# ensure we have Python 3 semantics from input, even in Python 2
try:
    input = raw_input
except NameError:
    pass

from sys import *
import pdb

def isHuman(id) :
   human_ids={
   9606 : 1 ,
   63221 : 1,
   741158 : 1
   }
   rval=0
   if id in human_ids : 
      rval=1
   return rval

usage = '''
usage: %s threshold lmat_lineage_file taxtree output_strain_select_version output_near_neighbor_map remapped ids ktax2ncbitax (optional)
where: input_taxonomy is in nodes.dmp format
''' % argv[0]

if len(argv) != 8 :
  print("args: ",len(argv))
  print(usage)
  exit(1)

mag_diff=100
tax_obs=argv[1]
ncbi_taxonomy=argv[2]
ncbi_rank =argv[3]
min_obs=int(argv[4])
print("min_obs:",min_obs)
thc_file=argv[5]
output_file=argv[6]
num_bins=int(argv[7])
ranks = {}
tax_hist_cnt={}
ignore_thc=False
try:
   thc_strm = open(thc_file)
   for line in thc_strm :
      line=line.rstrip()
      vals=line.split()
      tid=int(vals[0])
      cnt=int(vals[1])
      tax_hist_cnt.setdefault(tid,cnt)
except:
   ignore_thc=True
    
a = open(ncbi_rank)
for line in a :
  line=line.rstrip()
  tid,rank=line.split()
  ranks[int(tid)] = rank
ranks.setdefault(1,"life")

#get map: tax_id -> parent
a = open(ncbi_taxonomy)
names={}
parents = {}
read_cnt={}
children={}
while True :
  line1=a.readline()
  if not line1 :
     break
  if line1[0] == '#' :
     a.readline()
     a.readline()
     continue
  vals=line1.split()
  tid=int(vals[0])
  num_child=int(vals[1])
  name=a.readline() # ignore
  #names[tid] = name.rstrip()
  parent = int(vals[len(vals)-1].rstrip())
  parents[tid] = parent
  children[tid]={}
  for it in range(0,num_child,1) :
      assert it+2 < len(vals)
      child=vals[it+2] 
      children[tid][child] = 0

tax_obs_strm = open(tax_obs)
recall={}
store_rank_val={}
isEuk={}
cnt1=0
for line in tax_obs_strm :
   line=line.rstrip()
   t = line.split()
   tid = int(t[0])
   if ignore_thc :
      tax_hist_cnt.setdefault(tid,1)
   if tid not in tax_hist_cnt :
      print("really? ",tid)
      print(line)
      continue
   recall.setdefault(tid)
   curr_tid = parents[tid]
   kmer_cnt = tax_hist_cnt[tid]
   
   tmp_i1 = tid
   while True :
      if tmp_i1 == 2759 :
         isEuk.setdefault(tid,1)
         break
      if tmp_i1 == parents[tmp_i1] :
         break
      tmp_i1 = parents[tmp_i1]

   ## don't use plasmids to build null-models
   ## there seem to be a number of bacterial plasmids that don't have a matching chromosome
   ## so they are not specifically identified as plasmids (will need to fix this in next version)
   isIgnore = False
   tmp_i = tid
   while True :
      if tmp_i == 2 or tmp_i == 2157 or tmp_i == 28384 :
         isIgnore=True
         break
      if tmp_i == parents[tmp_i] :
         break
      tmp_i = parents[tmp_i]
   human_yes=False
   if isHuman(tid) :
      human_yes=True
   if (not human_yes and tid >= 10000000) or ( isIgnore and kmer_cnt < 100000) :
      continue

   verbose = False

   t.pop(0)
   while True :
      if (ranks[curr_tid] == "species" and human_yes) or \
         ranks[curr_tid] == "genus" or \
         ranks[curr_tid] == "family" or \
         ranks[curr_tid] == "order" or \
         ranks[curr_tid] == "class" or \
         ranks[curr_tid] == "phylum" or \
         ranks[curr_tid] == "kingdom" or \
         ranks[curr_tid] == "domain" or \
         ranks[curr_tid] == "life" :
   
         if curr_tid in store_rank_val   :
            lst=store_rank_val[curr_tid]
            for obi in range(0,num_bins,2) :
               pcnt = float(t[obi])
               num_obs = int(t[obi+1]) 
               fnd=False
               if verbose :
                  print("check",tid,curr_tid,kmer_cnt,pcnt,num_obs,lst)
               for it in range(len(lst)) :
                  (obs_lst,store_kmer_cnt) = lst[it]
                  for it1 in range(0,len(obs_lst),2):
                     store_pcnt = obs_lst[it1]
                     store_obs = obs_lst[it1+1]
                     chk_diff =kmer_cnt/store_kmer_cnt
                     if verbose :
                        print("qchkeck",chk_diff,kmer_cnt,store_kmer_cnt)
                     if chk_diff < mag_diff and pcnt >= store_pcnt :
                        if verbose :
                           print("did I get a replacement? ",chk_diff,kmer_cnt,store_kmer_cnt,pcnt,store_pcnt)
                        store_rank_val[curr_tid][it] = (pcnt,num_obs+store_obs,kmer_cnt)
                     if chk_diff < mag_diff : 
                        fnd=True
                        break
               if not fnd :
                  if verbose :
                     print("append",tid,curr_tid,kmer_cnt,pcnt,num_obs)
                  store_rank_val[curr_tid].append((t,kmer_cnt))
         else :      
            if verbose :
               print("initiate?",tid,curr_tid,t)
         
            store_rank_val.setdefault(curr_tid, [(t,kmer_cnt)]) 
         break

      if parents[curr_tid] == curr_tid : break
      curr_tid = parents[curr_tid]
   cnt1+=1

print("how much",cnt1)
if 561 in store_rank_val :
   merge_hack=store_rank_val[561]
if 620 in store_rank_val :
   merge_hack.extend(store_rank_val[620])
cnt=0
out_fh=open(output_file,"w")
out_fh.write(str(num_bins)+"\n")
## use 561 as a default value for the Euks, which are haywire now with high random values
## they likely need to be broken down by GC content type
qlst =  list(tax_hist_cnt.keys())
qlst.insert(0,562)
print(qlst)
once={}
save={}
saveid=[]
warn_cnt=0
for tid in qlst :
   if tid in once :
      continue
   once.setdefault(tid,1)
   if tid not in parents :
      continue
   curr_tid = parents[tid]
   use_val = []
   tid_kcnt=0
   if tid in tax_hist_cnt :
      tid_kcnt=tax_hist_cnt[tid]
   human_yes=False
   if isHuman(tid) :
      human_yes=True
   if tid >= 10000000 and not human_yes:  ## use k-mer count for parent
      tid_kcnt=tax_hist_cnt[curr_tid]

   isOther=False
   verbose=False
   while True :
      ## for now set null model of human to e coli equivalent
      ## it seems the random values for human are too high, possibly 
      ## due to the presence of small extreme gc content sections of the genome
      ## which will require a multi-factor model
      if curr_tid == 28384 :
         isOther=True
         break
      if (ranks[curr_tid] == "species" and human_yes) or \
         ranks[curr_tid] == "genus" or \
         ranks[curr_tid] == "family" or \
         ranks[curr_tid] == "order" or \
         ranks[curr_tid] == "class" or \
         ranks[curr_tid] == "phylum" or \
         ranks[curr_tid] == "kingdom" or \
         ranks[curr_tid] == "domain" or \
         ranks[curr_tid] == "life" :

         if curr_tid in store_rank_val :
            use_val = store_rank_val[curr_tid]
            if curr_tid == 561 or curr_tid == 620 :
               use_val = merge_hack 
               if verbose : 
                  print("debug check",tid_kcnt,use_val)

      if use_val != [] :
         if verbose :
            print("set here",tid,curr_tid,tid_kcnt,ranks[curr_tid],use_val)
         break
      if parents[curr_tid] == curr_tid : break
      curr_tid = parents[curr_tid]
   if isOther :
      use_val = merge_hack
   if tid == 9606 :
      use_val=store_rank_val[9606]    
   (rval_pcnt,rval_kcnt,rval_obs)=([0]*num_bins,[0]*num_bins,[0]*num_bins)
   (rval_pcnt1,rval_kcnt1,rval_obs1)=([1.0]*num_bins,[0]*num_bins,[0]*num_bins)
   close_match=[-1]*num_bins
   fndMatch=False
   for (oblst,kcnt) in use_val :
      diff_pcnt = tid_kcnt/kcnt
      if verbose :
         print("examine",diff_pcnt,pcnt,rval_pcnt,kcnt,tid_kcnt,tid)
      for it2 in range(0,len(oblst),2) :
         pcnt = oblst[it2]
         obs = oblst[it2+1]
         it = it2/2
         if diff_pcnt < mag_diff and pcnt > rval_pcnt[it] :
            rval_pcnt[it] = pcnt
            rval_obs[it] = obs
            rval_kcnt[it]= kcnt
            fndMatch=True
         if diff_pcnt < close_match[it] or close_match[it] == -1 :
            rval_pcnt1[it] = pcnt
            rval_obs1[it] = obs
            rval_kcnt1[it]= kcnt
            close_match[it] = diff_pcnt

   ## if there are no candidates within mag_diff of the representatives pick the closest
   if not fndMatch :
      (rval_pcnt,rval_kcnt,rval_obs)=(rval_pcnt1,rval_kcnt1,rval_obs1)


   if isHuman(tid) :
      use_rank = "genus"
   else :
      use_rank = ranks[curr_tid]
   if tid == 562 :
      def_euk_rval_pcnt = rval_pcnt
      def_euk_rval_obs = rval_obs
      def_euk_rval_kcnt = rval_kcnt
   if tid in isEuk and use_rank == "genus" :
      rval_pcnt = def_euk_rval_pcnt
      rval_obs = def_euk_rval_obs
      rval_kcnt = def_euk_rval_kcnt
   if tid == 1 :
      rval_pcnt = [1.0]*num_bins
   str_out=str(tid)+" "+str(use_rank)+"-"+str(curr_tid)
   save_rit,save_fit=-1,-1
   not_min=False
   for it in range(len(rval_pcnt)) :
      if int(rval_obs[it]) < min_obs : 
         not_min=True
         for rit in range(it-1,-1,-1):
            if int(rval_obs[rit]) >= min_obs :
               save_rit=rit
               break
         for fit in range(it+1,len(rval_pcnt),1):
            if int(rval_obs[fit]) >= min_obs :
               save_fit=fit
               break
         if save_rit >= 0 :
            d1=abs(it-save_rit)
         else :
            d1=num_bins+1
         if save_fit >= 0 :
            d2=abs(it-save_fit)
         else :
            d2=num_bins+1
         if d1 <= d2 and save_rit != -1 :
            rval_pcnt[it] = rval_pcnt[save_rit]
         elif save_fit != -1 :
            rval_pcnt[it] = rval_pcnt[save_fit]

   for it in range(len(rval_pcnt)) :
      str_out += " " + str(rval_obs[it])+" "+str(rval_pcnt[it])+" " +str(rval_kcnt[it])
   if save_rit == -1 and save_fit == -1 and not_min :
      warn_cnt+=1
   out_fh.write(str_out+"\n")
   cnt+=1
