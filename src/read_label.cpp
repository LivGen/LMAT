#include <unistd.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <list>
#include <stdlib.h>
#include <omp.h>
#include <queue>
#include "all_headers.hpp"
#include "TaxNodeStat.hpp"
#include "tid_checks.hpp"
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include <gzstream.h>

#include <version.h>

#define MMAP_SIZE 0
#define TID_T uint32_t

#define HUMAN_TAXID 9606


#define QUEUE_SIZE_MAX 2000000000


using std::fstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::istringstream;
using namespace metag;

static bool verbose=false;
static bool screenPhiXGlobal = true;
size_t perm_bytes_allocd;
static map<TID_T,string> gRank_table;

typedef map<TID_T,unsigned> cmap_t;
typedef pair<TID_T,float> ufpair_t;
typedef pair<TID_T,uint16_t> tax_elem_t;
typedef set<tax_elem_t> tax_data_t;
typedef pair<int16_t,tax_data_t> label_info_t;
typedef map<TID_T,TID_T> hmap_t;
typedef std::tr1::unordered_map<TID_T,float> ufmap_t;
typedef std::tr1::unordered_map<TID_T,vector<float> > uvfmap_t;
typedef std::tr1::unordered_map<TID_T,string> usmap_t;
typedef std::tr1::unordered_map<uint16_t, uvfmap_t> u_ufmap_t;
typedef std::tr1::unordered_map<uint16_t, usmap_t> u_usmap_t;
typedef pair<string,string> read_pair;


static bool gPERMISSIVE_MATCH = false;
vector <int> read_len_vec(1,0);
vector <int> read_len_avgs(1,0);

static std::tr1::unordered_map<int,string> gNum2rank;
static std::tr1::unordered_map<string,int> gRank2num;
static std::tr1::unordered_set<int> gLowNumPlasmid;

#define _USE_KPATH_IDS 0

#define isPlasmid(tid) (((tid >=10000000 && tid < 11000000) || (gLowNumPlasmid.find(tid) != gLowNumPlasmid.end())) ? true : false)



my_map tid_rank_map;
id_convback_map_t conv_map;


// this determines which style of tid reduction we want, default is
// the "comphrehensive reduction; requires the tid to rank mapping
// table.  -w option enables strain to species only reduction.
bool tid_map_is_strain_species = false;

bool badGenomes(TID_T tid) {
   bool isBad=false;
#if _USE_KPATH_IDS == 1
   switch(tid) {
      case 1154764:
      case 1218173:
         isBad=true;
         break;
      }
#else 
   switch(tid) {
      // comment in NCBI is that these genomes are likely HIV-1 but labeled only HIV
      // Thus their distinct lineage on another branch away from HIV-1 confounds HIV-1 labeling!
      // for now ignore sequences with this taxid.
      case 12721:
      case 693660:
         isBad=true;
         break;
   }
   
#endif
      return isBad;
}

// vec is assumed to be sorted
int closest(int value)
{
  unsigned i;

  for (i =0; i<read_len_avgs.size(); i++) {

    if (value <= read_len_avgs[i]) 
      return read_len_vec[i];

  }

  return (read_len_vec[i]);
  
}



static int getReadLen(int rl)
{
  int len = closest(rl);

  if (len > 0)
    return len;

  return 80;

}


typedef std::pair<TaxNode<TID_T>*,TID_T> tncpair_t;

static bool isAncestor(const TaxTree<TID_T>& tax_tree, TID_T prev_taxid /*ancestor */, TID_T curr_taxid /* descendant */ ) {
   bool isOk=false;
   vector<TID_T> path;
   TaxTree<TID_T>& tax_tree_tmp = const_cast<TaxTree<TID_T>&>(tax_tree);
   tax_tree_tmp.getPathToRoot(curr_taxid,path);
   for(unsigned p = 0; p < path.size(); ++p) {
      if(path[p] == prev_taxid) {
         isOk=true;
         break;
      }
   }
   return isOk;
}


struct SimpleCmp {
   bool operator()(const pair<TID_T,float>& a, const pair<TID_T,float>& b) const {
      return a.second > b.second;
   }
};

struct CmpDepth {
   CmpDepth(const hmap_t& imap) : _imap(imap) {}
   bool operator()(const pair<TID_T,float>& a, const pair<TID_T,float>& b) const {
      const int adepth = (*_imap.find(a.first)).second;
      const int bdepth = (*_imap.find(b.first)).second;
      return adepth > bdepth;
   }
   const hmap_t& _imap;
};

struct CmpDepth1 {
   CmpDepth1(const hmap_t& imap) : _imap(imap) {}
   bool operator()(TID_T a, TID_T b) const {
      const int adepth = (*_imap.find(a)).second;
      const int bdepth = (*_imap.find(b)).second;
      return adepth > bdepth;
   }
   const hmap_t& _imap;
};



enum nomatch_t { eReadTooShort=0, eNoDbHits, eLowScore };

static string nomatch2str(nomatch_t match) {
   string str = "error";
   switch(match) {
   case eReadTooShort :
      str="ReadTooShort";
      break;
   case eNoDbHits :
      str="NoDbHits";
      break;
   case eLowScore :
      str="LowScore";
      break;
   default:
      str="Error";
      break;
   }
   return str;
}

enum match_t { eDirectMatch, eMultiMatch, ePartialMultiMatch, eNoMatch, eNoLCAError };
static string match2str(match_t match) {
   string str = "error";
   switch(match) {
   case eDirectMatch : 
      str="DirectMatch";
      break;
   case eMultiMatch : 
      str="MultiMatch";
      break;
   case ePartialMultiMatch : 
      str="PartialMultiMatch";
      break;
   case eNoMatch : 
      str="NoMatch";
      break;
   case eNoLCAError : 
      str="LCA_ERROR";
      break;
   }
   return str;
}

static bool addToCandLineage(const ufpair_t cand, list<ufpair_t>& lineage, const hmap_t& dmap, const TaxTree<TID_T>& tax_tree) {
   bool addNode = false;
   if( lineage.empty() ) {
      addNode=true;
   } else {
      const hmap_t::const_iterator mtch0 = dmap.find(cand.first);
      unsigned cand_depth = 0;
      if( mtch0 != dmap.end()) {
         cand_depth = (*mtch0).second;
      }
      addNode = true;
      list<ufpair_t>::iterator it = lineage.begin();
      const list<ufpair_t>::iterator is = lineage.end();
      for(; it != is; ++it) {
         const TID_T taxid = (*it).first;
         const hmap_t::const_iterator mtch = dmap.find(taxid);
         unsigned chk_depth = 0;
         if( mtch != dmap.end() ) {
            chk_depth = (*mtch).second;
         }
         if( chk_depth > cand_depth && !isAncestor(tax_tree,cand.first,taxid) ) {
            addNode=false;
            break;
         } else if( chk_depth < cand_depth && !isAncestor(tax_tree,taxid,cand.first)) {
            addNode=false;
            break;
         } else if( chk_depth == cand_depth ) {
            addNode=false;
            break;
         }
      }
   }
   if( addNode ) {
      if(verbose) cout<<"Add to Lineage: "<<cand.first<<" "<<cand.second<<endl;
      lineage.push_back(cand);
   }
   return addNode; 
}

static bool cmpCompLineage(ufpair_t cand, const vector<ufpair_t>& lineage, set<TID_T>& no_good, float diff_thresh, const TaxTree<TID_T>& tax_tree) {
   const float undef = -10000;
   bool keep_going=true;
   for(unsigned i = 0; i < lineage.size(); ++i) {
      if( isAncestor(tax_tree, lineage[i].first, cand.first )) {
         break;
      } 
      if( lineage[i].second != undef &&  (lineage[i].second - cand.second) > diff_thresh  ) {
         if(verbose) cout<<"competing lineage is too far to care: "<<cand.first<<" "<<cand.second<<" "<<lineage[i].first<<" "<<lineage[i].second<<" "<<diff_thresh<<endl;
         keep_going=false;
         break;
      }
      if( (lineage[i].second - cand.second) <= diff_thresh ) {
         if(verbose) cout<<"competing lineage is too close: "<<cand.first<<" "<<cand.second<<" "<<lineage[i].first<<" "<<lineage[i].second<<" "<<diff_thresh<<endl;
         no_good.insert(lineage[i].first);
      }
   }
   return keep_going;
}

static pair<ufpair_t,match_t> 
findReadLabelVer2(const vector<ufpair_t>& rank_label, float diff_thresh, const TaxTree<TID_T>& tax_tree, const hmap_t& taxid2idx, 
                  list<ufpair_t>& cand_lin, const hmap_t& dmap, const ufmap_t& all_cand_set, const float topScore) {
   match_t match = eNoMatch;
   TID_T savePlasmidId;
   bool plasmidTopHit =false;
   unsigned lowest_depth = 0, highest_depth = 0;
   ufpair_t lowest=make_pair(0,0), highest = make_pair(0,0);
   signed lidx = -1;
   assert( cand_lin.empty() );
   bool linDone=false;
   for(signed i = rank_label.size()-1; i >= 0; --i) {
      if(verbose) {
         cout<<"huh? plasmid: "<<rank_label[i].first<<" "<<rank_label[i].second<<" "<< topScore<<endl;
         if( isPlasmid(rank_label[i].first) ) {
            cout<<"found plasmid: "<<rank_label[i].first<<" "<<rank_label[i].second<<" "<< topScore<<endl;
         }
      }
      if( rank_label[i].second >= topScore && isPlasmid(rank_label[i].first) ) {
         plasmidTopHit=true;
         savePlasmidId=rank_label[i].first;
         if(verbose) cout<<"Found top hit plasmid: "<<savePlasmidId<<endl;
      } 
      if( !linDone && !addToCandLineage(rank_label[i],cand_lin,dmap, tax_tree) )  {
         lidx = i;
         linDone=true;
      } else if( !linDone ) {
         const hmap_t::const_iterator mtch = dmap.find(rank_label[i].first);
         if( (*mtch).second> lowest_depth || (i == (signed)rank_label.size()-1) ) {
            lowest = rank_label[i];
            lowest_depth = (*mtch).second;
         }
         if( (*mtch).second < highest_depth || i == (signed)rank_label.size()-1 ) {
            highest = rank_label[i];
            highest_depth = (*mtch).second;
         }
      }
      //stick around to make sure we've checked for plasmids
      if( linDone && rank_label[i].second < topScore ) {
         break;
      }
   }
   set<TID_T> add_set;
   if( highest_depth != 0 ) {
      vector<TID_T> path;
      TaxTree<TID_T>& tax_tree_tmp = const_cast<TaxTree<TID_T>&>(tax_tree);
      tax_tree_tmp.getPathToRoot(highest.first,path);
      for(unsigned i = 0; i < path.size(); ++i) {
         add_set.insert( path[i] );
         ufmap_t::const_iterator mtch = all_cand_set.find( path[i] );
         if( mtch != all_cand_set.end() ) {
            ufpair_t val = make_pair( path[i], (*mtch).second);
            cand_lin.push_back(val);
         } else {
            //compute log_odds
            const float undef = -10000;
            cand_lin.push_back( make_pair(path[i], undef) );
         }
      }
   } 
   vector<ufpair_t> cand_lin_vec(cand_lin.size());
   list<ufpair_t>::const_iterator it = cand_lin.begin();
   const list<ufpair_t>::const_iterator is = cand_lin.end();
   for(unsigned i= 0; i < cand_lin_vec.size(); ++it, ++i) {
      cand_lin_vec[i] = *it;
   }
   CmpDepth cd(dmap);
   sort(cand_lin_vec.begin(),cand_lin_vec.end(),cd);

   set<TID_T> no_good;

   for(signed i = lidx; i >= 0; --i) {
      if( add_set.find( rank_label[i].first ) == add_set.end() ) {
         if( !cmpCompLineage(rank_label[i],cand_lin_vec,no_good,diff_thresh,tax_tree) )  {
            if(verbose) cout<<" competing taxid too low scoring, quit search "<<rank_label[i].first<<" "<<rank_label[i].second<<endl;
            break;
         }
      }
   }
   ufpair_t taxid_call;
   if( cand_lin.empty() && no_good.empty() ) {
      match = eNoMatch;
   } else if( !cand_lin.empty() && no_good.empty() ) {
      taxid_call = lowest; 
      match = eDirectMatch;
   } else {
      vector<ufpair_t> cand_vec(cand_lin.size());
      list<ufpair_t>::const_iterator it = cand_lin.begin(); 
      const list<ufpair_t>::const_iterator is = cand_lin.end(); 
      unsigned cnt= 0;
      for(cnt= 0; it != is; ++it, ++cnt) {
         if(verbose) cout<<"merging cand_lst "<<cnt<<" "<<(*it).first<<" "<<(*it).second<<endl;
         cand_vec[cnt] = *it;
      }  
      CmpDepth cd(dmap);
      sort(cand_vec.begin(),cand_vec.end(),cd);
      float sval = 0, max_val = -10000;
      pair<TID_T,bool> res = make_pair(0,false);
      int root_idx = -1;
      for(unsigned i = 0; i < cand_vec.size(); ++i) {
         const TID_T tax_i = cand_vec[i].first;
         max_val = std::max(cand_vec[i].second,max_val);
         sval += cand_vec[i].second;
         ++cnt; 
         if( no_good.find(tax_i) == no_good.end() ) {
            res = make_pair( cand_vec[i].first, true );
            root_idx = i;
            break;
         }
      }
      if(!res.second) {
         taxid_call = make_pair(0,-1);
         match = eNoLCAError; 
      }  else {
         const TID_T lca_tid = res.first;
         match = eMultiMatch; 
         if( all_cand_set.find(lca_tid) != all_cand_set.end() ) {
            assert(root_idx != -1);
            if( max_val < cand_vec[root_idx].second ) {
               match = ePartialMultiMatch;
               max_val = cand_vec[root_idx].second;
            }
         }
         taxid_call = make_pair(lca_tid, max_val);
      }
   }
   if( plasmidTopHit ) {
      if(verbose) cout<<"Check top hit plasmid to see if consistent with call: "<<endl;
      if( isAncestor(tax_tree,taxid_call.first,savePlasmidId) ) {
         taxid_call.first=savePlasmidId; 
         if(verbose) cout<<"YES: siwthc to plasmid: "<<taxid_call.first<<endl;
      }
   }
   if(verbose) cout<<"I'm confused: "<<taxid_call.first<<" "<<taxid_call.second<<endl;
   return make_pair(taxid_call,match);
}

void
fill_in_labels(const TaxTree<TID_T>& taxtree, vector<TID_T>& row, const tax_data_t& taxids, const hmap_t& idx2taxid, const hmap_t& taxid2idx) {
   for(unsigned tax_idx = 0; tax_idx < row.size(); ++tax_idx) {
      const hmap_t::const_iterator mtch = idx2taxid.find(tax_idx);
      assert(mtch != idx2taxid.end());
      const TID_T col_tax_val = (*mtch).second;
      const TaxTree<TID_T>::const_iterator col_tax_val_it = taxtree.find(col_tax_val); 
      if(col_tax_val_it == taxtree.end()) {
         cout<<"we have a taxtree mapping problem: "<<col_tax_val<<endl;
         continue;
      }
      TaxNode<TID_T>* col_tax_val_node = (*col_tax_val_it).second;
      // see if this id's genome count needs to be update
      //  based on its children   - can skip if its a leaf
      if( !col_tax_val_node->isLeaf() && row[tax_idx] == 0 ) {
         tax_data_t::const_iterator it = taxids.begin();
         const tax_data_t::const_iterator is = taxids.end();
         unsigned has_cnt_sum = 0;
         for(; it != is; ++it) {   
            const TID_T tax_id = (*it).first;
            const TaxTree<TID_T>::const_iterator tax_val_it = taxtree.find(tax_id); 
            if(tax_val_it == taxtree.end()) {
               cout<<"we have a taxtree mapping problem: "<<col_tax_val<<endl;
               continue;
            }
            TaxNode<TID_T>* tax_val_node = (*tax_val_it).second;
            // test if tax_id is a descendant of tax_val
            // save max value
            const hmap_t::const_iterator idx_mtch = taxid2idx.find(tax_id);
            assert(idx_mtch != taxid2idx.end());
            const TID_T tax_val_idx = (*idx_mtch).second;
            
            if( col_tax_val != tax_id && row[tax_val_idx] > 0 && tax_val_node->isAncestor(col_tax_val)) {
               // want to capture the descendant tax nodes closest to col_tax_val_node 
               has_cnt_sum = 1;
               break;
            }
         }
         if( has_cnt_sum > 0 ) {
            row[tax_idx] = has_cnt_sum;
            if(verbose) cout<<"Final Tax Node Val "<<col_tax_val<<" final transfer cnt for taxid="<<row[tax_idx]<<endl;
         }
      }
   }
   if(verbose) cout<<"fill_row: "; 
   for(unsigned tax_idx = 0; tax_idx < row.size(); ++tax_idx) {
      const hmap_t::const_iterator mtch = idx2taxid.find(tax_idx);
      assert(mtch != idx2taxid.end());
      const TID_T tax_val = (*mtch).second;
      if(verbose) cout<<"["<<tax_val<<","<<row[tax_idx]<<"]";
   }
   if(verbose) cout<<endl;
}

struct TCmp { 
   TCmp(const hmap_t& imap) : _imap(imap) {}
   bool operator()(const pair<TID_T,float>& a, const pair<TID_T,float>& b) const {
	if( fabs(a.second- b.second) < 0.001 ) {
		const int adepth = (*_imap.find(a.first)).second;
		const int bdepth = (*_imap.find(b.first)).second;
		return adepth < bdepth;
        } 
	return a.second < b.second; }
   const hmap_t& _imap;
};

struct ScoreOptions {
   ScoreOptions(hmap_t& imap) : _equal_kmer_vote(false), _strict_kmer_match(false), _prn_all(false), _imap(imap), _diff_thresh(1.0), _diff_thresh2(3.0) {}
   bool _equal_kmer_vote;
   bool _strict_kmer_match;
   bool _prn_all;
   hmap_t& _imap;
   float _diff_thresh, _diff_thresh2;
   u_ufmap_t _rand_hits; 
   u_usmap_t _rand_class; 
   bool _comp_rand_hits;
};

static void loadLowNumPlasmids(const string& file) {
   ifstream ifs_lst(file.c_str());
   if( !ifs_lst) {
      cerr<<"Unexpected reading error: "<<file<<endl;
      return;
   }
   TID_T pid;
   while(ifs_lst>>pid) {
      cout<<"debug: "<<pid<<endl;
      gLowNumPlasmid.insert(pid);
   }
}

static void loadRandHits(const string& file_lst, u_ufmap_t& rand_hits_all, u_usmap_t& rand_class_all) {
   ifstream ifs_lst(file_lst.c_str());
   if( !ifs_lst) {
      cerr<<"Unexpected reading error: "<<file_lst<<endl;
      return;
   }

   unsigned cnt = 0;
   gRank2num.insert( make_pair("no_rank",cnt) );

   gRank2num.insert( make_pair("ethnic",cnt++) );
   gRank2num.insert( make_pair("region",cnt++) );
   gRank2num.insert( make_pair("species",cnt++) );

   gRank2num.insert( make_pair("genus",cnt++) );
   gRank2num.insert( make_pair("family",cnt++) );
   gRank2num.insert( make_pair("order",cnt++) );
   gRank2num.insert( make_pair("class",cnt++) );
   gRank2num.insert( make_pair("phylum",cnt++) );
   gRank2num.insert( make_pair("kingdom",cnt++) );
   gRank2num.insert( make_pair("depth=0",cnt++) );

   cnt = 0;
   gNum2rank.insert( make_pair(cnt,"no_rank") );

   gNum2rank.insert( make_pair(cnt++,"ethnic") );
   gNum2rank.insert( make_pair(cnt++,"region") );
   gNum2rank.insert( make_pair(cnt++,"species") );

   gNum2rank.insert( make_pair(cnt++,"genus") );
   gNum2rank.insert( make_pair(cnt++,"family") );
   gNum2rank.insert( make_pair(cnt++,"order") );
   gNum2rank.insert( make_pair(cnt++,"class") );
   gNum2rank.insert( make_pair(cnt++,"phylum") );
   gNum2rank.insert( make_pair(cnt++,"kingdom") );
   gNum2rank.insert( make_pair(cnt++,"depth=0") );


   int read_len;
   string file;
   unsigned num_elements=0;
   while(ifs_lst>>read_len>>file) {
      bool first=true;
      const char* path = getenv("LMAT_DIR");
      if(path ) {
         file = string(path) + string("/") + file;
      } 
      cout<<"load: "<<read_len<<" "<<file<<endl;
      read_len_vec.push_back(read_len);

      ifstream ifs_pre(file.c_str());
      if( !ifs_pre) {
         cerr<<"Unexpected reading error: "<<file<<endl;
         continue;
      }
      ifs_pre.close();
      igzstream ifs(file.c_str());
      uvfmap_t& rand_hits = rand_hits_all[read_len];
      usmap_t& rand_class= rand_class_all[read_len];
      if(!first) {
         rand_hits.rehash(num_elements);
         rand_class.rehash(num_elements);
      }
      const int buff_size=20004;
      char buff[buff_size];
      ifs.getline(buff,buff_size);
      istringstream istrm(buff);
      int num_bins=0;
      istrm>>num_bins;
      assert(num_bins > 0);
      vector<float> save_ecoli(num_bins,0.5);
      while(ifs.getline(buff,buff_size)) {
         istringstream istrm(buff);
         TID_T taxid;
         float max_val=0;
         string class_str;
         istrm>>taxid>>class_str;
         size_t pos = class_str.find("-");
         assert(pos != string::npos);
         string val = class_str.substr(0,pos);
         if( val.size() >= 3 ) {
            // try to save cost of a string search
            // vectors got labeled this way (no_rank)
            // not really sure how best to handle these 
            if( val[0] == 'n' && val[1] == 'o' && val[2] == '_') {
               val="genus";
            }
         }
         list<unsigned> revisit;
         vector<float> cutoff(num_bins,0);
         for(unsigned bin=0; (signed)bin < num_bins; ++bin) {
            int num_obs,kmer_cnt;// not used at the moment
            istrm>>num_obs>>max_val>>kmer_cnt;
            // when random model failed to generate a single random match
            // it implies that GC content of the reference genome is likely too 
            // different, thus even "true" matches should not match this GC profile
            // for now raise the threshold before considering this match
            if( num_obs == 0 && kmer_cnt >= 100000 ) {
               max_val = 0.5;  
               cutoff[bin] = max_val;
            // for small genomes 0 values may be due to to lack of observations
            // so fill with closest observed bin
            }  else if( num_obs == 0 && kmer_cnt < 100000) {
               revisit.push_back(bin);
            }
            if( num_obs > 0 ) {
               cutoff[bin] = max_val;
               if( taxid == 562 ) {
                  save_ecoli[bin]=cutoff[bin];
               }
            }
            if( taxid == 28384 ) {
               val="genus";
               cutoff=save_ecoli;
            }
         }
         if( !revisit.empty() ) {
            list<unsigned>::const_iterator it = revisit.begin();
            const list<unsigned>::const_iterator is = revisit.end();
            for(; it != is; ++it) {
               signed j = *it - 1;
               unsigned i = *it + 1;
               while(j >= (signed)0 || i < cutoff.size() ) {
                  float a_val=0.0,b_val=0.0;
                  if( j >= 0 ) {
                     a_val=cutoff[j];
                  }
                  if( i < cutoff.size() ) {
                     b_val=cutoff[i];
                  }
                  if( a_val > 0 && b_val > 0 ) {
                     cutoff[*it] = max(a_val,b_val);
                  } else if( a_val > 0 ) {
                     cutoff[*it] = a_val;
                  } else if( b_val > 0 ) { 
                     cutoff[*it] = b_val;
                  }
                  if( cutoff[*it] > 0 ) {
                     break;
                  } 
                  --j;
                  ++i;
               }
               // wow no observations
               // should really catch this earlier.
               if( cutoff[*it] <= 0) {
                  cutoff[*it] = 0.5;
                  //if(verbose) cout<<"WarningDebug: "<<*it<<" "<<cutoff[*it]<<endl;
               }
            }
         }
         rand_hits[taxid]=cutoff;
         rand_class[taxid]=val;
         if(first) ++num_elements;
      }
      first=false;
   } 
   sort(read_len_vec.begin(),read_len_vec.end());
   int i; 
   read_len_avgs.resize(0);
   for (i = 1; i < (signed)read_len_vec.size(); i++) {
     read_len_avgs.push_back((read_len_vec[i-1] + read_len_vec[i]) / 2);
   }
}

static float log_odds_score(float label_prob, float random_prob) {
      //if( label_prob >= 1 ) label_prob = 0.9999;
      //if( random_prob >= 1 ) random_prob = 0.9999;
      //const float denom = (1-label_prob)*random_prob; // must be > zero
      //const float numer =  label_prob*(1-random_prob); // must be > zero
      //if(verbose)  cout<<"log_odds_calc: "<<numer<<" "<<denom<<endl;
      const float numer = label_prob;
      const float denom = random_prob <= 0 ? 0.00001 : random_prob;
      const float log_odds = log( numer / denom );
      return log_odds;
}

pair<ufpair_t,match_t>
construct_labels(const TaxTree<TID_T>& tax_tree, const vector<label_info_t>& label_vec, const list<TID_T>& taxid_lst, const hmap_t& tax2idx, const hmap_t& idx2taxid, 
                 ofstream& ofs, size_t mer_len, const ScoreOptions& sopt, int bin_sel, int min_valid_kmers, int min_fnd_kmers) {
   const unsigned num_tax_ids = taxid_lst.size();
   vector<bool> any_kmer_match(label_vec.size(),false) ;
   unsigned cnt_fnd_kmers=0;
   vector<vector<TID_T> > label_matrix(label_vec.size());
   uint16_t cand_kmer_cnt = 0;
   uint16_t debug_bad_cand_kmer_cnt = 0;
   for(unsigned pos = 0; pos < label_vec.size(); ++pos) {
      if(label_vec[pos].first >= 0 ) ++cand_kmer_cnt;
      else {
         if(verbose) cout<<"Debug bad kmers pos: "<<pos<<" "<<label_vec[pos].first<<endl;
         ++debug_bad_cand_kmer_cnt;
      }
      label_matrix[pos].resize(num_tax_ids,0);
      tax_data_t::const_iterator it = label_vec[pos].second.begin();
      const tax_data_t::const_iterator is = label_vec[pos].second.end();
      unsigned any=0;
      for(; it != is; ++it) {
         const TID_T tax_id = (*it).first;
         const uint16_t has_cnt = (*it).second;
         if( tax2idx.find(tax_id) == tax2idx.end() ) {
            cout<<"HOUSTON WE HAVE A PROBLEM: "<<tax_id<<endl;
         }
         const unsigned idx = (*(tax2idx.find(tax_id))).second;
         label_matrix[pos][idx] = has_cnt;
         ++any;
         if(verbose) cout<<"check: "<<pos<<" "<<idx<<" "<<tax_id<<" "<<has_cnt<<endl;
      }   
      any_kmer_match[pos] = !label_vec[pos].second.empty();
      if( any_kmer_match[pos] ) {
         ++cnt_fnd_kmers;
      }
   }
   if( (int)cnt_fnd_kmers < min_fnd_kmers ) { 
      return make_pair(make_pair(0,-1),eNoMatch);
   }

   //don't let LMAT print out labels for reads that are below the threshold, otherwise
   //counting won't add up when you compare individual read labels and the summary counts
   if( cand_kmer_cnt < min_valid_kmers ) return make_pair(make_pair(0,-1),eNoMatch);

   float defRandMod=0.1;
   const  int cand_kmer_cnt_match = getReadLen(cand_kmer_cnt);
   vector<float> rank_first(num_tax_ids,0); 
   u_ufmap_t::const_iterator mtch = sopt._rand_hits.find(cand_kmer_cnt_match);
   const bool useRandMod = (mtch == sopt._rand_hits.end()) ? false : true;
   const uvfmap_t& rand_hits = (*mtch).second;
   u_usmap_t::const_iterator e_mtch = sopt._rand_class.find(cand_kmer_cnt_match);
   const usmap_t& equiv_class = (*e_mtch).second;

   ufmap_t all_cand_set;
   bool hasHuman=false;
   std::tr1::unordered_map<string,float> track;
   std::tr1::unordered_map<string,float> counter;
   for(unsigned tax_idx = 0; tax_idx < num_tax_ids; ++tax_idx) {
      float found_genome_cnt = 0;
      const TID_T taxid = (*(idx2taxid.find(tax_idx))).second;
      if( isHuman(taxid) ) hasHuman=true;
      if(verbose) cout<<"col: "<<tax_idx<<" "<<taxid;
      bool noMatch=false;
      for(unsigned pos = 0; pos < label_vec.size(); ++pos) {
         if( label_matrix[pos][tax_idx] > 0 ) {
            found_genome_cnt += 1;
            if(verbose) cout<<" pos="<<pos<<" "<<label_vec[pos].first<<" "<<label_matrix[pos][tax_idx]<<" taxid="<<taxid<<" "<<found_genome_cnt<<endl;
         }
      }
      if( !noMatch ) {
         const float label_prob = (float)found_genome_cnt / (float) cand_kmer_cnt;
         rank_first[tax_idx] = label_prob; 
         if(verbose) cout<<" "<<label_prob<<endl;
      }
      float random_prob = 0.5;
      if( !useRandMod ) {
         random_prob = 0.1; // arbitray threshold, lets one ignore need for null model
      } else if( rand_hits.find(taxid) != rand_hits.end()) {
         const vector<float>& val_vec = (*(rand_hits.find(taxid))).second;
         const float val = val_vec[bin_sel];
         random_prob = val+0.0001;
      } else {
         cerr<<"ERROR, ALL TAXIDS MUST HAVE NULL MODELS: "<<taxid<<" "<<cand_kmer_cnt_match<<" "<<cand_kmer_cnt<<endl;
         random_prob = 1.0; // basically try to ignore these tax ids for now.
      } 
      if( useRandMod ) {
         usmap_t::const_iterator chk = equiv_class.find( taxid );
         assert(chk != equiv_class.end());
         const string cval = (*chk).second;
         if( track.find(cval) == track.end() ) {
            if(verbose) cout<<"Confirm first max: "<<taxid<<" "<<cval<<" "<<random_prob<<" "<<track[cval]<<" "<<endl;
            track[ cval ] = random_prob;
            ///make sure that null model for rank is not smaller than the lower rank ones 
            ///in theory this should not happen, but with pruning of taxids in the db, you can
            // get unexpected combinations, where a higher order rank is represented for one brank but not another
            const int cval_rank = gRank2num[cval];
            for(signed ti = cval_rank - 1; ti >= 0; ti--) {
               const string cval_lower = gNum2rank[ti];
               track[ cval] = std::max(track[cval],track[cval_lower]);
            }
         } else {
            if(verbose) cout<<"Confirm max: "<<taxid<<" "<<cval<<" "<<random_prob<<" <? "<<track[cval]<<" "<<endl;
            track[ cval] = std::max(random_prob,track[cval]);
            const int cval_rank = gRank2num[cval];
            for(signed ti = cval_rank - 1; ti >= 0; ti--) {
               const string cval_lower = gNum2rank[ti];
               track[ cval] = std::max(track[cval],track[cval_lower]);
            }
         }
      }
      if(verbose) cout<<" cand kmer_cnt "<<cand_kmer_cnt_match<<" in="<<cand_kmer_cnt<<" "<<debug_bad_cand_kmer_cnt<<" "<<label_matrix.size()<<endl;
   }
   vector<ufpair_t> rank_label(num_tax_ids,make_pair(0,0)); 
   bool fndPhiX=false;
   float log_sum=0.0, pos_log_sum = 0.0, top_score = 0.0, phiXscore=0.0;
   unsigned sig_hits = 0, pos_sig_hits=0;
   for(unsigned tax_idx = 0; tax_idx < num_tax_ids; ++tax_idx) {
         const TID_T taxid = (*(idx2taxid.find(tax_idx))).second;
         const float label_prob = rank_first[tax_idx];
         float random_prob = defRandMod;
         string cval;
         if( useRandMod ) {
            usmap_t::const_iterator chk = equiv_class.find( taxid );
            assert(chk != equiv_class.end());
            cval = (*chk).second;
            random_prob = track[cval] ;
         }
         const float log_odds = useRandMod ? log_odds_score(label_prob,random_prob) : label_prob;
         if(verbose) cout<<"track score: "<<taxid<<" "<<cval<<" "<<label_prob<<" "<<random_prob<<" "<<log_odds<<endl;
         rank_label[tax_idx] = make_pair(taxid,log_odds);
         all_cand_set.insert( rank_label[tax_idx] );
         log_sum += log_odds;
         sig_hits++;
         if( log_odds > 0) {
            pos_sig_hits++;
            pos_log_sum += log_odds;
         }
         if(verbose) cout<<"Sig match chec: "<<taxid<<" "<<rank_label[tax_idx].second<<" "<<label_prob<<" "<<random_prob<<" "<<endl;
         if( screenPhiXGlobal && isPhiX(taxid) ) {
            phiXscore=log_odds;
            fndPhiX=true;
         }

         if( tax_idx == 0 || log_odds > top_score ) {
            top_score = log_odds;
         }
   }
   match_t mtype;
   ufpair_t best_guess = make_pair(0,0);
   // bypass normal analysis and set this read to phiX
   if( screenPhiXGlobal && phiXscore >= top_score && fndPhiX) {
      best_guess = make_pair(ART_SEQ_TID,phiXscore);
      mtype = eDirectMatch;
      ofs<<(-1)<<" "<<(-1)<<" "<<cand_kmer_cnt<<"\t";
      ofs<<best_guess.first<<" "<<best_guess.second;
      ofs<<"\t";
      ofs<<best_guess.first<<" "<<best_guess.second<<" "<<match2str(mtype);
      ofs<<endl;
   } else {
      list<ufpair_t> valid_cand;
      pair<ufpair_t,match_t> res = make_pair(make_pair(0,0),eNoMatch);   
      string matchType=match2str(res.second);
      unsigned use_sig_hits = 0;
      float log_avg;
      const unsigned min_pos_examples=3;
      if( pos_sig_hits > min_pos_examples ) {
         use_sig_hits = pos_sig_hits;
         log_avg = pos_log_sum / (float)pos_sig_hits;
      } else {
         use_sig_hits = sig_hits;
         log_avg = sig_hits > 0 ? log_sum / (float)sig_hits : 0;
      }
      float log_std=0, max_score = 0;
      bool first=true;
      for(unsigned tax_idx = 0; tax_idx < num_tax_ids; ++tax_idx) {
         if( rank_label[tax_idx].second > 0 ) { 
            if( pos_sig_hits > min_pos_examples) {
               const float val = log_avg - rank_label[tax_idx].second;
               log_std += (val*val);
            }
            if( rank_label[tax_idx].second > max_score || first ) {
               max_score = rank_label[tax_idx].second;
               first=false;
            }
         }
         if( pos_sig_hits <= min_pos_examples ) {
            const float val = log_avg - rank_label[tax_idx].second;
            log_std += (val*val);
         }
      }
      float stdev1 = use_sig_hits > 1 ? sqrt(log_std/(use_sig_hits-1)) : 0;
      if( use_sig_hits > 0 ) {
         if( hasHuman ) {
            for(unsigned tax_idx = 0; tax_idx < num_tax_ids; ++tax_idx) {
               //float found_genome_cnt = 0;
               const TID_T taxid = (*(idx2taxid.find(tax_idx))).second;
               if( isHuman(taxid) ) {
                     rank_label[tax_idx].second += (sopt._diff_thresh2*stdev1);
               }
            } 
         }
         TCmp tcmp(sopt._imap);
         sort(rank_label.begin(),rank_label.end(),tcmp);
         ofs<<log_avg<<" "<<stdev1<<" "<<cand_kmer_cnt<<"\t";
         stdev1 *= sopt._diff_thresh;
         
         res = findReadLabelVer2(rank_label,stdev1,tax_tree,tax2idx,valid_cand,sopt._imap,all_cand_set,top_score);    
         if(sopt._prn_all) {
            bool prn=false;
            for(signed i = rank_label.size()-1; i >= 0; --i) {
               if( rank_label[i].second >= 0  || verbose ) {
                  ofs<<" "<<rank_label[i].first<<" "<<rank_label[i].second;
                  prn=true;
               }
            }
            if(!prn) {
               ofs<<"-1 -1";
            }
            ofs<<"\t";
         }
         matchType=match2str(res.second);
      }
      if( res.second == eDirectMatch) {
         best_guess = res.first;
         ofs<<best_guess.first<<" "<<best_guess.second<<" "<<matchType;
      }  else if( res.second == eMultiMatch || res.second == ePartialMultiMatch ) {
         if( !sopt._prn_all) {
            list<ufpair_t>::const_iterator it = valid_cand.begin();
            const list<ufpair_t>::const_iterator is = valid_cand.end();
            for( ; it != is; ++it) {
               ofs<<" "<<(*it).first<<" "<<(*it).second;
            }
            if(valid_cand.empty() ) {
               ofs<<"-1 -1";
            }
            ofs<<"\t";
         }
         best_guess = res.first;
         ofs<<best_guess.first<<" "<<best_guess.second<<" "<<matchType;
      }
      else if( res.second == eNoMatch ) {
         ofs<<-1<<" "<<-1<<" "<<matchType;
      } else {
         cerr<<"Unexpected match type"<<endl;
         ofs<<-1<<" "<<-1<<" "<<"Unmatched";
      }
      ofs<<endl;
      mtype = res.second;
   }
   return make_pair(best_guess,mtype);
}

#define ENCODE(t, c, k,gc_cnt,tot_cnt) \
switch (c) { \
case 'a': case 'A': t = 0; break; \
case 'c': case 'C': t = 1; break; \
case 'g': case 'G': t =2; break; \
case 't': case 'T': t = 3; break; \
default: k = 0; gc_cnt=0; tot_cnt=0; continue; \
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

static

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
    int valid_kmers=0, gc_cnt=0, valid_gc_cnt = 0, valid_tot_cnt = 0, tot_cnt = 0;
    for (j = 0; j < slen; j++) {
        register int t;
        const char base=str[j];
        ENCODE(t, base, k,gc_cnt,tot_cnt);
        forward = ((forward << 2) | t) & mask;
        reverse = ((kmer_t)(t^3) << highbits) | (reverse >> 2);
        if(base == 'g' || base=='G'||base=='c'||base=='C') {
         ++gc_cnt;
         ++tot_cnt;
        } else if(base == 'a' || base=='A'||base=='t'||base=='T') {
         ++tot_cnt;
        } 
      
        //gc_cnt = (base == 'g' || base=='G'||base=='c'||base=='C') ? gc_cnt+1 : gc_cnt;
        if (++k >= (signed)klen) {
           valid_kmers++;
           valid_gc_cnt += gc_cnt; 
           valid_tot_cnt += tot_cnt; 
           if(verbose) cout<<"check j:"<<valid_gc_cnt<<" "<<valid_tot_cnt<<" "<<gc_cnt<<" "<<tot_cnt<<" "<<j<<" "<<base<<endl;
           gc_cnt=0;
           tot_cnt=0;
           kmer_id = (forward < reverse) ? forward : reverse;
           if( no_dups.find(kmer_id) != no_dups.end() ) continue;
            /* zero based position of forward k-mer is (j-klen+1) */
            /* zero based position of reverse k-mer is (slen-j-1) */
           const int pos = j-klen+1;
           /* kmer_lookup(kmer); do k-mer lookup here... */
           label_vec[pos].first = 0; // marks the position as having a valid k-mer (for case where n-masked reads are used)
           if(verbose) cout<<"debug valid: "<<j<<" "<<pos<<endl;
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
           bool seenHuman=false; 
           while( h->next() ) {
              TID_T tid = h->taxid();
              if(isHuman(tid) && seenHuman) continue;
              else if(isHuman(tid) && !seenHuman) {
                tid=HUMAN_TAXID;
                seenHuman=true;
              }
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
         } else {
            //save_spec_rep.insert(make_pair(stid,make_pair(stid,stid_cnt)));
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
   float gc_pcnt = ((float)valid_gc_cnt / (float)valid_tot_cnt) * 100.0;
   int bin_sel = gc_pcnt / 10; // better to not hard code this 
   if(verbose) cout<<"GC info: "<<gc_pcnt<<" "<<bin_sel<<" "<<valid_gc_cnt<<" "<<valid_tot_cnt<<endl;
   return make_pair(valid_kmers,bin_sel);
}

void proc_line(const TaxTree<TID_T>& tax_tree, int ri_len, string &line, int k_size, INDEXDB<DBTID_T> *table, ofstream &ofs, float threshold, const ScoreOptions& sopt, uint16_t max_count
               ,map<TID_T,int>& track_taxids, map<nomatch_t,int>& track_nomatch, map<TID_T,float>& track_tscores, float min_label_score, int min_kmer, int min_fnd_kmer) {

     if(ri_len < 0 || ri_len > (signed)line.length()) {
         cout<<"unexpected ri_len value: "<<ri_len<<endl;
         return;
     } else if( ri_len < k_size ) {
         ofs<<"-1 -1 -1"<<"\t-1 -1\t"<<ri_len<<" "<<k_size<<" ReadTooShort"<<endl;
        if( track_nomatch.find(eReadTooShort) != track_nomatch.end()) {
            track_nomatch[eReadTooShort] += 1;
        } else {
            track_nomatch[eReadTooShort] = 1;
        }
     } else {
        vector<label_info_t> label_vec(ri_len-k_size+1,make_pair(-1,tax_data_t()));
        list<TID_T> taxid_lst; 
        hmap_t tax2idx, idx2tax;
        const pair<int,int> res = retrieve_kmer_labels(table, line.c_str(), ri_len, k_size,label_vec,taxid_lst,tax2idx,idx2tax, sopt._imap, tax_tree, max_count);
        const int valid_kmers = res.first;
        const int bin_sel = res.second;
        // needed to check if this is empty due to short read first!!!
        if( valid_kmers < min_kmer ) {
           ofs<<"-1 -1 -1"<<"\t-1 -1\t"<<valid_kmers<<" "<<min_kmer<<" ReadTooShort"<<endl;
           if( track_nomatch.find(eReadTooShort) != track_nomatch.end()) {
               track_nomatch[eReadTooShort] += 1;
           } else {
               track_nomatch[eReadTooShort] = 1;
           }
        } else if( !taxid_lst.empty() ) {  
           pair< ufpair_t, match_t> mtch = construct_labels(tax_tree,label_vec,taxid_lst,tax2idx,idx2tax,ofs,k_size,sopt,bin_sel,min_kmer,min_fnd_kmer); 
           if(mtch.second == eNoMatch && valid_kmers < min_kmer ) {
              ofs<<"-1 -1 -1"<<"\t-1 -1\t"<<valid_kmers<<" "<<min_kmer<<" ReadTooShort"<<endl;
              if( track_nomatch.find(eReadTooShort) != track_nomatch.end()) {
                  track_nomatch[eReadTooShort] += 1;
              } else {
                  track_nomatch[eReadTooShort] = 1;
              }
           } else if(mtch.second == eNoMatch ) {
              if( track_nomatch.find(eNoDbHits) != track_nomatch.end()) {
                  track_nomatch[eNoDbHits] += 1;
              } else {
                  track_nomatch[eNoDbHits] = 1;
              }
           } else if( mtch.first.second >= min_label_score && valid_kmers >= min_kmer) {
               if( track_taxids.find(mtch.first.first) == track_taxids.end()) {
                  track_taxids[mtch.first.first] = 1;
                  track_tscores[mtch.first.first] = mtch.first.second;
               } else {
                  track_taxids[mtch.first.first] += 1;
                  track_tscores[mtch.first.first] += mtch.first.second;
               }
           } else if( mtch.first.second < min_label_score ) {
              if( track_nomatch.find(eLowScore) != track_nomatch.end()) {
                  track_nomatch[eLowScore] += 1;
              } else {
                  track_nomatch[eLowScore] = 1;
              }
           } 
           
        } else {
           ofs<<"-1 -1 "<<valid_kmers<<"\t-1 -1\t"<<ri_len<<" "<<k_size<<" NoDbHits"<<endl;
           if( track_nomatch.find(eNoDbHits) != track_nomatch.end()) {
               track_nomatch[eNoDbHits] += 1;
           } else {
               track_nomatch[eNoDbHits] = 1;
           } 
        }
     }
  }


   size_t *split_file(int n_threads, ifstream &file)
   {
     size_t * arr = new size_t[n_threads+1];

     arr[0] = 0;

     file.seekg(0, ios::end);
  size_t end = file.tellg();

  for (size_t i = 1; (signed)i<n_threads; i++)  {

    file.seekg( i * (end / n_threads ));
    
    string junk;
    getline(file, junk);
    
    arr[i] = file.tellg();
  }

  arr[n_threads] = end;

  return arr;

}

void usage(char *execname)
{
  cout << "LMAT version " << LMAT_VERSION  << "\n";
  cout << "Usage:\n" ;
  cout << execname << " -d <input db file (list)> -i <query fasta file> -t <number of threads> -o <output path> [-l <human bias>]\n";
  cout << "[-v <pct cutoff>]  -c <tax tree file> -k <kmer size> [-z:verbose] [-r:rank table]\n";
  cout << "[-m <rank/tid-map-file>] [-g <tid-cutoff>] [-w <with-strain-species-map> (affects -m option)]\n";
  cout << "[-h:turn phiX screening off]\n";
 
}


int main(int argc, char* argv[]) 
{
   char c = '\0';
   int k_size=-1;
   int n_threads = 0;

   bool restore = true;

   float threshold = 0.0, min_score = 0.0;
   int min_kmer = 35, min_fnd_kmer = 1;

   string rank_map_file,rank_ids, kmer_db_fn, query_fn, ofname, ofbase, tax_tree_fn, tax_tree_options, depth_file, rand_hits_file, rank_table_file, id_bit_conv_fn;
   string low_num_plasmid_file;
   hmap_t imap;	
   ScoreOptions sopt(imap);

   size_t mmap_size = 0;
   bool fastq=false; 
   uint16_t max_count = ~0;
   bool prn_read = true;



   while ((c = getopt(argc, argv, "u:ahn:j:b:ye:w:pk:c:v:k:i:d:l:t:r:sm:o:x:f:g:z:qV")) != -1) {

      switch(c) {
      case 'h':
        screenPhiXGlobal=false;
         break;
      case 'r':
        low_num_plasmid_file = optarg;
        break;
      case 'f':
        id_bit_conv_fn = optarg;
        break;
      case 'j':
         min_kmer = atoi(optarg);
         break;
      case 'z':
         min_fnd_kmer = atoi(optarg);
         break;
      case 'u':
         rank_ids = optarg;
         break;
      case 'x':
         min_score = atof(optarg);
         break;
      case 'a':
         prn_read=false;
         break;
      case 'w':
	      rank_map_file = optarg;
	      break;
      case 's':
         gPERMISSIVE_MATCH=true;
         break;
      case 'n':
         rand_hits_file=optarg;
         break;
      case 'b':
         sopt._diff_thresh = atof(optarg);
         break;
      case 'l':
         sopt._diff_thresh2 = atof(optarg);
         break;
      case 'y':
         verbose = true;
         break;
      case 'e':
         depth_file = optarg;
         break;
      case 'q':
         fastq=true;
         break;
      case 'p':
         sopt._prn_all = true;
         break;   
      case 'm':
	      rank_table_file = optarg;
        break;
      case 't':
        n_threads = atoi(optarg);
        omp_set_num_threads(n_threads);
        break;
      case 'v':
         threshold = atof(optarg);
         break;
      case 'c':
         tax_tree_fn = optarg;
	      break;
      case 'k':
         k_size = atoi(optarg);
         break;
      case 'g':
	      max_count = atoi(optarg);
         break;
      case 'i':
         query_fn = optarg;
         break;
      case 'd':
         kmer_db_fn = optarg;
         break;
      case 'o':
         ofbase = optarg;
         break;
    case 'V':
      cout << "LMAT version " << LMAT_VERSION  << "\n";
	exit(0);
      default:
         cout << "Unrecognized option: "<<c<<", ignore."<<endl;
      }
   }
   if (depth_file == "") cout << "depth_file\n";
   if (ofbase == "") cout << "ofbase\n";
   if (n_threads == 0) cout << "n_threads\n";
   if (kmer_db_fn == "") cout << "kmer_db_fn\n";
   if (query_fn == "") cout << "query_fn\n";

   if (depth_file == "" ||  ofbase == "" || n_threads == 0 || kmer_db_fn == "" || query_fn == "")  {
     cout<<ofbase<<" "<<n_threads<<" "<<kmer_db_fn<<" "<<query_fn<<" "<<depth_file<<endl; 
     usage(argv[0]);
     return -1;

   }
   if (!restore && k_size == -1)  {
     cout << "missing kmer size!\n";
     usage(argv[0]);
     return -1;

   }

   cout << "Start kmer DB load\n";
   INDEXDB<DBTID_T> *taxtable;


#if USE_BOOST == 1
     bip::managed_mapped_file mfile(bip::open_only, kmer_db_fn.c_str());
     std::pair<INDEXDB<DBTID_T>*, std::size_t> ret = mfile.find<INDEXDB<TID_T>>("KmerDB");

     taxtable = ret.first;
     k_size = taxtable->get_kmer_length();
     cout << "k size:  " << k_size  <<   endl ;
     taxtable->conv_ptrs();
#else


#if WITH_PJMALLOC == 1

   if (restore) {
      
     perm(&taxtable, sizeof(taxtable));
     if( mopen(kmer_db_fn.c_str(), "r", mmap_size) != 0 ) {
         cerr<<"Error: unable to open kmer db ["<<kmer_db_fn<<"]"<<endl;
         return -1;
     }
     if (k_size < 1)
       k_size = taxtable->get_kmer_length();

     cout << "num kmers: " << taxtable->size() << " - " << k_size  <<   endl ;


   } else 
    
#endif

{


#if (USE_SORTED_DB == 0)
     taxtable = new INDEXDB<DBTID_T>;

     ifstream qifs(kmer_db_fn.c_str());
     if( !qifs ) {
       cerr<<"Unable to open: "<<kmer_db_fn<<endl;
       return -1;


     }




     string fname;
     
     while(qifs>>fname) {
       cout<<"register file: "<<fname<<endl;
       taxtable->registerFile(fname.c_str());
     }
     taxtable->ingest();
#endif


   }

#endif
   if( low_num_plasmid_file.length() > 0 ) {
      loadLowNumPlasmids(low_num_plasmid_file);
   }

   sopt._comp_rand_hits = (rand_hits_file.length() == 0);
   if( !sopt._comp_rand_hits ) {
      loadRandHits(rand_hits_file,sopt._rand_hits,sopt._rand_class);
   }
   if( k_size <= 0 ) {
      cerr<<"Unable to read database, k_size="<<k_size<<endl;
      return -1;
   }



   string line;

   omp_lock_t buffer_lock;
   std::queue < read_pair > read_buffer_q;
   omp_init_lock(&buffer_lock);

   bool finished;

   ofstream ofs;

   if (rank_table_file.length() > 0) {
     if(max_count == ~0) {
       cout << "Need to set -h <tid-cutoff> to use rank file map!\n";
     } else {
       FILE * rmfp = fopen(rank_table_file.c_str(), "r");
     
       uint32_t src, dest;
       
       while (fscanf(rmfp,"%d%d", &src, &dest) > 0) {
	      tid_rank_map[src] = dest;
       }
       
       fclose(rmfp);
       
     }
	 
   }
   if( rank_map_file.size() > 0) {
      ifstream ifs1(rank_map_file.c_str());
      TID_T tid;
      string rank;
      while(ifs1>>tid>>rank) {
         gRank_table.insert(make_pair(tid,rank));
      }
   }




   cout<<"Read taxonomy tree: "<<tax_tree_fn<<endl;
   TaxTree<TID_T> tax_tree(tax_tree_fn.c_str());

   cout<<"Read taxonomy depth: "<<depth_file<<endl;
   ifstream ifs1(depth_file.c_str());
   if(!ifs1) {
	cerr<<"unable to open: "<<depth_file<<endl;
	return -1;
   }
   TID_T taxid,depth;
   while(ifs1>>taxid>>depth) {
	sopt._imap[taxid] = depth;
   }
   

   if (id_bit_conv_fn.length() > 0) {
     cout << "Loading map file,\n";
     FILE * tfp = fopen(id_bit_conv_fn.c_str(), "r");
     if( !tfp ) {
         cout<<"Unable to read 16-bit map file:"<<id_bit_conv_fn<<endl;
         return -1;
     }

     uint32_t src;
     uint16_t dest;

     while (fscanf(tfp,"%d%hd", &src, &dest) > 0) {

       conv_map[dest] = src;
     }
     fclose(tfp);
   }



   StopWatch clock;
   clock.start();
   vector<map<TID_T,float> > track_tscoreall(n_threads); 
   vector<map<TID_T,int> > track_matchall(n_threads); 
   vector<map<nomatch_t,int> > track_nomatchall(n_threads); 

   size_t read_count_in =0;
   size_t read_count_out = 0;


   string read_buff, hdr_buff, save_hdr;

   bool in_finished = false;

   ifstream tmpstream;

   istream ifs(cin.rdbuf());
   if (query_fn != "-") {
     tmpstream.open(query_fn.c_str());

     if(!tmpstream) {
	  cerr<<"did not open for reading: "<<query_fn<<endl;
	  
	  exit(-1);
	  
     }  
     ifs.rdbuf(tmpstream.rdbuf());
   }

   




#pragma omp parallel shared(k_size, query_fn, ofbase, taxtable, tax_tree, sopt,  prn_read,track_matchall,track_nomatchall,track_tscoreall,min_score,min_kmer, in_finished, read_count_in, read_count_out, min_fnd_kmer, ifs)  private(finished, ofs, ofname, line, read_buff, hdr_buff, save_hdr)
  {

    finished = false;

    ofname = ofbase;
    std::stringstream outs;
    outs << omp_get_thread_num();
    ofname += outs.str();
    ofname += ".out" ;
    ofs.open(ofname.c_str());
    bool eof = false;

    while (!finished)   {
      if ((in_finished == false) && (omp_get_thread_num() == 0)) {
	      int j = 0 ;
	      int queue_size = 0;
	      omp_set_lock(&buffer_lock);
	      queue_size = read_buffer_q.size();
	      omp_unset_lock(&buffer_lock);
	      string last_hdr_buff;
	      while (queue_size < QUEUE_SIZE_MAX && j< 2* n_threads && (!in_finished)) {
            eof = !getline(ifs, line);
   	      if (eof) {
	            in_finished = true;
	    
		    if(verbose) cout << line.size() << " line length\n";
		    line = "";
	  }


	  if (line[0] == '>' || (fastq && line[0] == '@') ) {

	    last_hdr_buff = hdr_buff;
	    // skip the ">"                                                        
	    hdr_buff=line.substr(1,line.length()-1);
	  }

	  if (line[0] != '>' && line.length() > 1 && !fastq) {
	    read_buff += line;
	    line = "";
	  }

	  if( fastq && line[0] != '@' && line[0] != '+' && line[0] != '-' ) {
	    read_buff += line;
	    line = "";
	  }
	  if( ((line[0] == '>' || in_finished) || (fastq && (line[0] == '+' ||
	     line[0] == '-'))) && read_buff.length() > 0 ) {

	    omp_set_lock(&buffer_lock);
	    
	    if (in_finished)
		read_buffer_q.push(read_pair(read_buff, hdr_buff));
	    else
		read_buffer_q.push(read_pair(read_buff, last_hdr_buff));

	    read_count_in++;
	    omp_unset_lock(&buffer_lock);

	    read_buff="";

	    j ++;
	    
	    if(fastq) eof = !getline(ifs, line); // skip quality values for now       

	  }

	  if (in_finished) {
	    
	    cout << read_count_in << " reads in\n";
	    break;

	  }
	}
	
      }
      read_buff = "";
      save_hdr = "";
      omp_set_lock(&buffer_lock);
      if (!read_buffer_q.empty()) {
         read_pair in_pair = read_buffer_q.front();
         read_buff = in_pair.first;
         save_hdr = in_pair.second;
         read_buffer_q.pop();
         read_count_out++;

      }

      omp_unset_lock(&buffer_lock);
      if (read_buff.length() > 0) {
         if(save_hdr[0] == '\0') {
           ostringstream ostrm;
           ostrm<<"unknown_hdr:"<<read_count_out;
           save_hdr=ostrm.str();
         }
         ofs<<save_hdr<<"\t";
         if( prn_read ) {
           ofs<<read_buff<<"\t";
         } else {
           ofs<<"X"<<"\t";
         }

         int thread = omp_get_thread_num();

         map<TID_T,int>& track_match = track_matchall[thread];
         map<nomatch_t,int>& track_nomatch = track_nomatchall[thread];
         map<TID_T,float>& track_tscore = track_tscoreall[thread];

        proc_line(tax_tree, read_buff.length(), read_buff, k_size, taxtable, 
		  ofs, threshold,sopt, max_count, track_match, track_nomatch, 
		  track_tscore,min_score, min_kmer, min_fnd_kmer);
	     read_buff="";
	
      }
      if ((read_count_in == read_count_out) && in_finished)
	   finished = true;
    }
    ofs.close();
  }



   map<TID_T,int>  merge_count;
   map<TID_T,float>  merge_score;
   map<nomatch_t,int>  nomatch_merge_count;
   for(unsigned thread = 0; thread < track_tscoreall.size();  ++thread) {
      map<TID_T,float>& gt = track_tscoreall[thread];
      map<TID_T,float>::const_iterator it = gt.begin();
      const map<TID_T,float>::const_iterator is = gt.end();
      for(; it != is; ++it) {
         TID_T tid = (*it).first;
         float score = (*it).second;
         if( merge_score.find(tid) == merge_score.end() ) {
            merge_score.insert( make_pair(tid,score));
         } else {
            merge_score[tid] += score;
         }
      }
      map<TID_T,int>& gtcnt = track_matchall[thread];
      map<TID_T,int>::const_iterator it1 = gtcnt.begin();
      const map<TID_T,int>::const_iterator is1 = gtcnt.end();
      for(; it1 != is1; ++it1) {
         TID_T tid = (*it1).first;
         int cnt = (*it1).second;
         if( merge_count.find(tid) == merge_count.end() ) {
            merge_count.insert( make_pair(tid,cnt));
         } else {
            merge_count[tid] += cnt;
         }
      }
      map<nomatch_t,int>& nogtcnt = track_nomatchall[thread];
      map<nomatch_t,int>::const_iterator it2 = nogtcnt.begin();
      const map<nomatch_t,int>::const_iterator is2 = nogtcnt.end();
      for(; it2 != is2; ++it2) {
         nomatch_t mt = (*it2).first;
         int cnt = (*it2).second;
         if( nomatch_merge_count.find(mt) == nomatch_merge_count.end() ) {
            nomatch_merge_count.insert( make_pair(mt,cnt));
         } else {
            nomatch_merge_count[mt] += cnt;
         }
      }
   }
   set<TID_T> cand_tid;
   vector< pair<TID_T,float> > sort_val(merge_score.size());
   map<TID_T,float>::const_iterator it = merge_score.begin();
   map<TID_T,float>::const_iterator is = merge_score.end();
   for(unsigned i = 0; it != is; ++it, ++i) {
      sort_val[i] = (*it);
      const TID_T tid = sort_val[i].first;
      if( cand_tid.find( tid ) ==  cand_tid.end()) {
         cand_tid.insert(tid);
      }
   }
   const unsigned buff_size = 200000;
   char buff[buff_size];
   map<TID_T,string> save_id;
   if(verbose) cout<<"rank read: "<<rank_ids<<endl;
   ifstream tax_strm(rank_ids.c_str());
   while(tax_strm.getline(buff,buff_size)) {
      string proc(buff);
      char* val = strtok(buff,"=,");
      while( val != NULL ) {
         if( strcmp(val,"taxid")==0 ) {
            val = strtok(NULL,"=,");
            istringstream istrm(val);
            TID_T cid;
            istrm>>cid;
            if( cand_tid.find(cid) != cand_tid.end() ) {
               size_t pos = proc.rfind('\t');
               string id=proc.substr(pos+1,proc.length()-pos);
               save_id.insert( make_pair(cid,id) );
            }
            break;
         }
         val = strtok(NULL,"=,");
      }
   }
   ostringstream ostrm;
   ostrm<<ofbase<<"."<<min_score<<"."<<min_kmer<<".fastsummary";
   ofstream sum_ofs(ostrm.str().c_str());
   if( !sum_ofs ) {
      cout<<"Could not open for writing "<<ostrm.str()<<endl;
      return -1;
   }
   sort(sort_val.begin(),sort_val.end(),SimpleCmp());
   for(unsigned i= 0; i < sort_val.size(); ++i) {
      const TID_T tid = sort_val[i].first;
      assert( merge_count.find( tid ) != merge_count.end());
      const int cnt = (*merge_count.find(tid)).second;
      const float wght_cnt = sort_val[i].second;
      string str_id = save_id[tid];
      sum_ofs<<wght_cnt<<"\t"<<cnt<<"\t"<<tid<<"\t"<<str_id<<endl;
   }
   ostringstream nomostrm;
   nomostrm<<ofbase<<"."<<min_score<<"."<<min_kmer<<".nomatchsum";
   ofstream nom_ofs(nomostrm.str().c_str());
   if( !nom_ofs ) {
      cout<<"Could not open for writing "<<nomostrm.str()<<endl;
      return -1;
   }
   map<nomatch_t,int>::const_iterator it3 = nomatch_merge_count.begin();
   map<nomatch_t,int>::const_iterator is3 = nomatch_merge_count.end();
   for(; it3 != is3; ++it3) {
      const string id = nomatch2str( (*it3).first );
      const int count = (*it3).second;
      nom_ofs<<id<<"\t"<<count<<endl;
   }
   cout << "query time: " << clock.stop() << endl;

   return 0; 
}
