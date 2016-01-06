#ifndef _RKMER_HPP
#define _RKMER_HPP

#include <tr1/unordered_set>


typedef map<TID_T,unsigned> cmap_t;
typedef map<TID_T,TID_T> hmap_t;
typedef pair<TID_T,uint16_t> tax_elem_t;
typedef set<tax_elem_t> tax_data_t;
typedef pair<int16_t,tax_data_t> label_info_t;

extern my_map tid_rank_map;
extern id_convback_map_t conv_map;
extern bool tid_map_is_strain_species;
extern bool verbose;
extern bool gPERMISSIVE_MATCH;
extern map<TID_T,string> gRank_table;

struct CmpDepth1 {
   CmpDepth1(const hmap_t& imap) : _imap(imap) {}
   bool operator()(TID_T a, TID_T b) const {
      const int adepth = (*_imap.find(a)).second;
      const int bdepth = (*_imap.find(b)).second;
      return adepth > bdepth;
   }
   const hmap_t& _imap;
};

bool badGenomes(TID_T tid) {
   bool isBad=false;
   switch(tid) {
      // comment in NCBI is that these genomes are likely HIV-1 but labeled only HIV
      // Thus their distinct lineage on another branch away from HIV-1 confounds HIV-1 labeling!
      // for now ignore sequences with this taxid.
      case 12721:
      case 693660:
         isBad=true;
         break;
   }
   return isBad;
}


#define ENCODE(t, c, k) \
switch (c) { \
case 'a': case 'A': t = 0; break; \
case 'c': case 'C': t = 1; break; \
case 'g': case 'G': t =2; break; \
case 't': case 'T': t = 3; break; \
default: k = 0; continue; \
}

/*

   Given a sequence of nucleotide characters,

   break it into canonical k-mers in one pass.

   Nucleotides are encoded with two bits in

   the k-mer. Any k-mers with ambiguous characters

   are skipped.

   str:  DNA sequence (read)

   slen: DNA sequence length in nucleotides

   klen: k-mer length in nucleotides

*/

pair<int,int> retrieve_kmer_labels(INDEXDB<DBTID_T>* table, const char* str, const int slen, const kmer_t klen, 
                          vector<label_info_t>& label_vec, list<TID_T>& taxid_lst, hmap_t& tax2idx, hmap_t& idx2tax, 
                          const hmap_t& dmap, const TaxTree<TID_T>& tax_tree, uint16_t max_count) 
   {
    int j; /* position of last nucleotide in sequence */
    int k = 0; /* count of contiguous valid characters */
    int highbits = (klen-1)*2; /* shift needed to reach highest k-mer bits */
    kmer_t mask = ((kmer_t)1 << klen*2)-1; /* bits covering encoded k-mer */
    kmer_t forward = 0; /* forward k-mer */
    kmer_t reverse = 0; /* reverse k-mer */
    kmer_t kmer_id; /* canonical k-mer */
    set<kmer_t> no_dups;
    cmap_t leaf_track;
    int valid_kmers=0, gc_cnt=0;
    for (j = 0; j < slen; j++) {
        register int t;
        const char base=str[j];
        ENCODE(t, base, k);
        gc_cnt = (base == 'g' || base=='G'||base=='c'||base=='C') ? gc_cnt+1 : gc_cnt;
        forward = ((forward << 2) | t) & mask;
        reverse = ((kmer_t)(t^3) << highbits) | (reverse >> 2);
        if (++k >= (signed)klen) {
           valid_kmers++;
           kmer_id = (forward < reverse) ? forward : reverse;
            /* zero based position of forward k-mer is (j-klen+1) */
            /* zero based position of reverse k-mer is (slen-j-1) */
           const int pos = j-klen+1;
           /* kmer_lookup(kmer); do k-mer lookup here... */
           label_vec[pos].first = 0; // marks the position as having a valid k-mer (for case where n-masked reads are used)
           if(verbose) cout<<"debug valid: "<<j<<" "<<pos<<endl;
           if( no_dups.find(kmer_id) != no_dups.end() ) continue;
           no_dups.insert(kmer_id);

           TaxNodeStat<DBTID_T> *h = new TaxNodeStat<DBTID_T>(*table);
           if(verbose) cout<<"lookup kmer at posit:"<<j<<endl;

#if (TID_SIZE == 16)
           h->begin(kmer_id, tid_rank_map,  max_count, tid_map_is_strain_species, &conv_map);
#else
           h->begin(kmer_id, tid_rank_map,  max_count, tid_map_is_strain_species);
#endif

           unsigned dcnt = 0, mtch = 0;
           list<TID_T> obs_tids;
           while( h->next() ) {
              TID_T tid = h->taxid();
              if( tid == 20999999 || badGenomes(tid)) continue;
              uint16_t ng = h->taxidCount();
               if(dcnt==0) {
                  if( ng <= 0 ) {
                     cout<<"Warning unexpected value for ng: "<<ng<<endl;
                     ng=1;
                  }
                  label_vec[pos].first = ng;
               }
               //collect the unique set of tax ids
               //const uint16_t pr_cnt = h->present();
               const uint16_t pr_cnt = 1;
               obs_tids.push_back(tid);
               if(gPERMISSIVE_MATCH) {
                  label_vec[pos].second.insert( make_pair(tid,pr_cnt) );
                  if( tax2idx.find(tid) == tax2idx.end() ) {
                     const unsigned idx = taxid_lst.size();
                     tax2idx[tid] = idx;
                     idx2tax[idx] = tid;
                     taxid_lst.push_back(tid);
                  }
               }
                  // Note: ng should not change with multiple calls to next here
               if(verbose) {
                  if( dcnt == 0 ) cout<<"debug kmer "<<kmer_id<<" gc="<<ng;
                  cout<<" ["<<pos<<" "<<tid<<" "<<pr_cnt<<"]" ;
               }
               dcnt++;
               ++mtch;
           }
           vector<TID_T> obs_tids_vec(obs_tids.size());
           list<TID_T>::const_iterator it = obs_tids.begin();
           const list<TID_T>::const_iterator is = obs_tids.end();
           for(unsigned i =0; it != is; ++it, ++i) {
              obs_tids_vec[i] = *it;
           }
           CmpDepth1 cd(dmap);
           sort(obs_tids_vec.begin(),obs_tids_vec.end(),cd);
           if( gPERMISSIVE_MATCH ) {
              int last_depth = -1;
              for(unsigned i = 0; i < obs_tids_vec.size(); ++i) {
                const TID_T tid=obs_tids_vec[i];
                const int depth = (*dmap.find(tid)).second;
                if( depth == 0 ) {
                  break;
                }
                if( last_depth == depth || last_depth == -1 ) {
                     TaxTree<TID_T>& tax_tree_tmp = const_cast<TaxTree<TID_T>&>(tax_tree);
                     vector<TID_T> path;
                     tax_tree_tmp.getPathToRoot(tid,path);
                     for(unsigned p = 0; p < path.size(); ++p) {
                        const TID_T ptid = path[p];
                        if( verbose ) cout<<"lineage tids added: "<<ptid<<" for "<<tid<<" pos="<<pos<<" "<<p<<endl;
                        label_vec[pos].second.insert( make_pair(ptid,1) );
                        if( tax2idx.find(ptid) == tax2idx.end() ) {
                           const unsigned idx = taxid_lst.size();
                           tax2idx[ptid] = idx;
                           idx2tax[idx] = ptid;
                           if( verbose ) cout<<"lineage tids added: "<<ptid<<" for "<<tid<<" pos="<<pos<<endl;
                           taxid_lst.push_back(ptid);
                        }
                     }
                } else {
                  break;
                }
             }
          } else {
             std::tr1::unordered_set<TID_T> non_leaf; 
             for(unsigned i = 0; i < obs_tids_vec.size(); ++i) {
                const TID_T tid=obs_tids_vec[i];
                if( non_leaf.find(tid) == non_leaf.end() ) {
                   const uint16_t pr_cnt = 1;
                   const int depth = (*dmap.find(tid)).second;
                   if(verbose) cout<<"this is a non leaf tid:"<<tid<<" "<<depth<<endl;
                   label_vec[pos].second.insert( make_pair(tid,pr_cnt) );
                   if( leaf_track.find(tid) != leaf_track.end()) {
                      leaf_track[tid] += 1;
                   } else {
                      leaf_track.insert(make_pair(tid,1));
                   }
                   if( tax2idx.find(tid) == tax2idx.end() ) {
                      const unsigned idx = taxid_lst.size();
                      tax2idx[tid] = idx;
                      idx2tax[idx] = tid;
                      taxid_lst.push_back(tid);
                   }
                   TaxTree<TID_T>& tax_tree_tmp = const_cast<TaxTree<TID_T>&>(tax_tree);
                   vector<TID_T> path;
                   tax_tree_tmp.getPathToRoot(tid,path);
                   //const int depth = (*dmap.find(tid)).second;
                   for(unsigned p = 0; p < path.size(); ++p) {
                     const TID_T ptid = path[p];
                     non_leaf.insert(ptid);
                   }
               } else {
                  if(verbose) cout<<"skip this non leaf tid:"<<tid<<endl;
               }
             }
          }


          if( verbose ) cout<<"num taxids: "<<taxid_lst.size();
          if(dcnt > 0 && verbose ) cout<<" end k-mer lookup"<<endl;
          if(mtch == 0 && verbose ) cout<<" no k-mer matches "<<endl;
          delete h;
       }
   }
   if( ! gPERMISSIVE_MATCH ) {
      map<TID_T,pair<TID_T,unsigned> > save_spec_rep;
      cmap_t::const_iterator cb = leaf_track.begin();
      cmap_t::const_iterator cs = leaf_track.end();
      for(; cb != cs; ++cb) {
          const TID_T stid = (*cb).first;
          const unsigned stid_cnt = (*cb).second;
          if( gRank_table.find(stid) != gRank_table.end() && gRank_table[stid] == "strain" ) {
             vector<TID_T> path;
             TaxTree<TID_T>& tax_tree_tmp = const_cast<TaxTree<TID_T>&>(tax_tree);
             tax_tree_tmp.getPathToRoot(stid,path);
             for(unsigned p = 0; p < path.size(); ++p) {
                const TID_T ptid = path[p];
                if( gRank_table.find(ptid) != gRank_table.end() && gRank_table[ptid] == "species" ) {
                  if( save_spec_rep.find(ptid) == save_spec_rep.end() ) {
                     save_spec_rep.insert(make_pair(ptid,make_pair(stid,stid_cnt)));
                  } else if( stid_cnt > save_spec_rep[ptid].second ) {
                     save_spec_rep[ptid] = make_pair(stid,stid_cnt);
                  }
                  break;
               }
             } 
         }
      }
      std::tr1::unordered_set<TID_T> rep_strain;
      map<TID_T,pair<TID_T,unsigned> >::const_iterator sb1=save_spec_rep.begin(); 
      map<TID_T,pair<TID_T,unsigned> >::const_iterator se1=save_spec_rep.end(); 
      for( ; sb1 != se1; ++sb1) {
         TID_T species_id= (*sb1).first;
         TID_T strain_id= (*sb1).second.first;
         rep_strain.insert( strain_id );
         if(verbose) cout<<"Representative strain: "<<strain_id<<" "<<species_id<<" "<<(*sb1).second.second<<endl;
      }
      for(unsigned pos = 0; pos < label_vec.size(); ++pos) {
         if( label_vec[pos].first >= 0 ) {
               set<tax_elem_t>::const_iterator sb=label_vec[pos].second.begin(); 
               const set<tax_elem_t>::const_iterator se=label_vec[pos].second.end(); 
               for(; sb != se; ++sb) {
                  TID_T tid = (*sb).first;
                  if( rep_strain.find(tid) != rep_strain.end() || gRank_table[tid] != "strain"  ) {
                     TaxTree<TID_T>& tax_tree_tmp = const_cast<TaxTree<TID_T>&>(tax_tree);
                     vector<TID_T> path;
                     tax_tree_tmp.getPathToRoot(tid,path);
                     for(unsigned p = 0; p < path.size(); ++p) {
                        const TID_T ptid = path[p];
                        if( verbose ) cout<<"lineage tids added: "<<ptid<<" for "<<tid<<" pos="<<pos<<" "<<p<<endl;
                        label_vec[pos].second.insert( make_pair(ptid,1) );
                        if( tax2idx.find(ptid) == tax2idx.end() ) {
                           const unsigned idx = taxid_lst.size();
                           tax2idx[ptid] = idx;
                           idx2tax[idx] = ptid;
                           if( verbose ) cout<<"lineage tids added: "<<ptid<<" for "<<tid<<" pos="<<pos<<endl;
                           taxid_lst.push_back(ptid);
                        }
                     }
                  }
               }
         }
      }
   }

   float gc_pcnt = ((float)gc_cnt / (float)j) * 100.0;
   int bin_sel = gc_pcnt / 10; // better to not hard code this 
   if(verbose) cout<<"GC info: "<<gc_pcnt<<" "<<bin_sel<<endl;
   return make_pair(valid_kmers,bin_sel);
}

#endif
