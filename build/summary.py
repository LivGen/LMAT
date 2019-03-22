#!/usr/bin/env python

import sys,os
import copy
from operator import itemgetter

summfile=sys.argv[1]
rankfile=sys.argv[2]
fsummfile=sys.argv[3]
plasmidFile=sys.argv[4]
out_base=sys.argv[5]
rank_calls=sys.argv[6]

globalPlasmids={}
def loadPlasmid(file) :
   fh1=open(file) 
   for ln in fh1 :
      ln=ln.rstrip()
      globalPlasmids.setdefault(int(ln),1)

loadPlasmid(plasmidFile)

def isPlasmid(id) :
   if id >= 10000000 or globalPlasmids.has_key(id) :
      return True
   return False

def loadRank(rfile) :
   save={}
   fh1=open(rfile) 
   for ln in fh1 :
      ln=ln.rstrip()
      vl=ln.split()
      save.setdefault(int(vl[0]),vl[1])
   return save

def loadFastSumm(rfile) :
   save={}
   fh1=open(rfile)
   for ln in fh1 :
      ln=ln.rstrip()
      vl=ln.split()
      save.setdefault(int(vl[2]),ln)
   return save

def cntTabs(arr) :
   cnt=0
   for iter in range(len(arr)) :
      if arr[iter] != '' :
         break
      cnt+=1 
   return cnt

def getNode(arr) :
   cnt,iter=0,0
   while iter < len(arr) :
      if arr[iter] != '' :
         break
      iter += 1
   return arr[iter],int(arr[iter+1]),int(arr[iter+2]),float(arr[iter+3])

rankMap=loadRank(rankfile)
fsumMap=loadFastSumm(fsummfile)

wrdcnt,kmercnt,rdcnt,call_kmercnt={},{},{},{}
def loadTree(fh) :
   names={}
   parent={}
   child={}
   parent.setdefault(1,1)
   lines=[(1,0)]
   last_tab_cnt=-1
   for line in fh :
      line=line.rstrip()
      vals=line.split('\t')
      if vals[0]=='Name' :
         continue
      num_tabs=cntTabs(vals)
      pn,cnode,val2,val3= getNode(vals)
      names[cnode] = pn
      rdcnt[cnode]=val2
      wrdcnt[cnode]=val3
      ## find parent for current node
      while lines != [] : 
         pnode,last_tab_cnt=lines[0]
         if num_tabs > last_tab_cnt :
            if not child.has_key(pnode) :
               child.setdefault(pnode,[])
            child[ pnode ].append( cnode )
            parent[ cnode ] = pnode
            break
         else :
            lines.pop(0)
      lines.insert(0,(cnode,num_tabs))
   return child,names

def summNode(cnode,call_type,child) :
   tot_wrdcnt,tot_rdcnt,tot_kmercnt=0,0,0
   if (rankMap.has_key(cnode) and rankMap[cnode] == call_type and not isPlasmid(cnode)) or (isPlasmid(cnode) and call_type == "plasmid") :
      tot_wrdcnt = wrdcnt[cnode]
      tot_rdcnt = rdcnt[cnode]
      the_call=cnode

      save_strains = []
      lopen=[]
      if child.has_key(cnode) :
         lopen=copy.deepcopy(child[ cnode ])
      while lopen != [] :
         alt = lopen.pop(0)
         ## for species calls don't use plasmid reads - report these separately 
         if (call_type == "species" and not isPlasmid(alt)) or (call_type != "species" ) and rdcnt[alt] > 0 :
            tot_wrdcnt += wrdcnt[alt]
            tot_rdcnt += rdcnt[alt]
         if call_type == "species" and rankMap.has_key(alt) and rankMap[alt] == "strain" and not isPlasmid(alt) and rdcnt[alt] > 0:
            save_strains.append( alt )

         if child.has_key(alt) :
            tlst=child[alt]
            for it in range(len(tlst)) :
               nd = tlst[it]
               lopen.append(nd)
         
      if save_strains != [] :
         slst = sorted( save_strains, key=lambda val : wrdcnt[val], reverse=True )
         ## for now this returns just the top strain call
         the_call=slst[0]
   return (cnode,the_call,tot_wrdcnt,tot_rdcnt)
      
def bread_first_traverse(call_type,child) :
   save_calls=[]
   lopen=[1]
   while lopen != []:
      cnode=lopen.pop(0)
      ## some plasmids are ranked as "species", still skip these during species call
      if (call_type == "plasmid" and isPlasmid(cnode)) or (rankMap.has_key(cnode) and call_type == rankMap[cnode] and not isPlasmid(cnode)) :
         result=summNode(cnode,call_type,child)
         if result[3] > 0 :
            save_calls.append(result)
      else :
         ## should be higher order
         lst=[]
         if child.has_key(cnode) :
            lst=child[ cnode ]
         for it in range(len(lst)) :
            nd = lst[it]
            lopen.insert(0,nd)
   return save_calls

def findPeak(arr) :
   fndFirstLocalMin = False
   copyVal=-1
   for it in range(1,len(arr)-1,1) :
      if not fndFirstLocalMin and arr[it-1][1] >= arr[it][1] and arr[it][1] < arr[it+1][1] :
         fndFirstLocalMin = True
      if fndFirstLocalMin and arr[it-1][1] <= arr[it][1] and arr[it][1] > arr[it+1][1] :
         copyVal=arr[it][0] 
         break
   return copyVal

      

def loadKmerStats(ifile,rank) :
   fh=open(ifile)
   saveDistr=False
   distr=[] 
   hold={}
   tot_kcnt,tid,kval,kcnt=-1,-1,-1,-1
   while True :
      ln = fh.readline()
      ln=ln.rstrip()
      if ln == "" or (ln.find("taxid=") != -1 and ln.find("distinct_kmer_cnt=") != -1) :
         if distr != [] :
            val=findPeak(distr)
            hold.setdefault(tid,{})
            hold[tid].setdefault(kval,(val,kcnt,tot_kcnt))
         if ln == "" :
            break
         saveDistr=False
         distr=[]
         vals=ln.split('=')
         tid=int(vals[1].split(' ')[0])
         if rankMap.has_key(tid) and rank == rankMap[tid] :
            kcnt=int(vals[2].split(' ')[0])
            kval=int(vals[3].split(' ')[0])
            tot_kcnt=int(vals[4].split(' ')[0])
            saveDistr=True
      elif saveDistr :
         vals=ln.split(' ')
         distr.append( (int(vals[2]),int(vals[3])) )
 
   return hold 

def doPrn(save_calls,outh,names,kcov) :
   rep=sorted ( save_calls, key = itemgetter(2), reverse=True )
   pstr="% of Reads, Avg Read Score, Weighted Read Count (WRC), Read Count (RC), Original WRC, Original RC, Name, Taxid"
   outh.write(pstr+"\n")

   rc_sum=0
   for val in rep :
      rc_sum+=val[3]
   for val in rep :
      rep_id=val[0]
      call_id=val[1]
      owrc,orc=-1,-1
      if fsumMap.has_key(call_id) : 
         v1=fsumMap[call_id].split('\t')
         prnName=v1[3]
         owrc,orc=v1[0],v1[1]
      else :
         prnName=names[call_id]
      wrc,rc=val[2],val[3]
      if rc == 0 :
         print "which one?",call_id,val
      avg=float(wrc)/float(rc)
      tot_pcnt=float(rc)/float(rc_sum)
      pstr=str(tot_pcnt)+"\t"+str(avg)+"\t"+str(wrc)+"\t"+str(rc)+"\t"+str(owrc)+"\t"+str(orc)+"\t"+prnName+"\t"+str(call_id)
      pstr += "\t"+str(rep_id)
      if kcov.has_key(rep_id) :
         for kv in kcov[rep_id].keys() :
            pstr += "\t"+str(kv)+","+str(kcov[rep_id][kv][0])+","+str(kcov[rep_id][kv][1])+","+str(kcov[rep_id][kv][2])
      outh.write(pstr+"\n")

fh=open(summfile)
lchild,lname=loadTree(fh)
for ranktype in rank_calls.split(',') :
   outfile=open(out_base+"."+ranktype,"w")
   redun_file=summfile+"."+ranktype+"_kmer_cov"
   k_mer_cov={}
   if os.path.isfile(redun_file) :
      k_mer_cov=loadKmerStats(redun_file,ranktype)
   save_calls=bread_first_traverse(ranktype,lchild)
   doPrn(save_calls,outfile,lname,k_mer_cov)
