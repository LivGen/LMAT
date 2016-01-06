#ifndef __TAXNODESTAT__
#define __TAXNODESTAT__

#include "metag_typedefs.hpp"

#if (USE_SORTED_DB == 1)
#define INDEXDB SortedDb
#include "SortedDb.hpp"
#else
#define INDEXDB TaxTable
#include "TaxTable.hpp"
#endif

#include <ext/hash_map>
#include <vector>
#include <queue>

typedef  __gnu_cxx::hash_map<uint16_t,uint32_t> id_convback_map_t;



#define BAD_TAXID 3
			     
namespace metag {

using namespace std;

template<class tid_T>
class TaxNodeStat {
public:




  TaxNodeStat(INDEXDB<tid_T> &table)
    : m_table(&table), m_offset(0), m_taxid_count(0), m_calls_to_next(1)  {} 



  void begin(kmer_t kmer,  id_convback_map_t  *p_map = NULL) {

    m_calls_to_next = 0; 


    mp_idmap=p_map;
    
    // load map file

    m_kmer = kmer;

    bool good = m_table->begin_(kmer, m_taxid_count, m_offset, m_page);
    if (!good) {
      m_calls_to_next = 1;
      m_taxid_count = 0;
    }

  }

  void begin(kmer_t kmer,  my_map &p_map, int tid_cut = 0, bool strainspecies = 0,  id_convback_map_t  *pp_map = NULL) { 

    m_calls_to_next = 0;
    
    // load map file
    mp_idmap = pp_map;

    m_kmer = kmer;
    
    bool good = m_table->begin_(kmer, m_taxid_count, m_offset, m_page);
    if (!good) {
      m_calls_to_next = 1;
      m_taxid_count = 0;
    }
    else { 

      if (tid_cut > 0 && m_taxid_count > tid_cut)  {

	if (p_map.size() == 0) {
	  m_filtered_list.push_back(1);
	  m_taxid_count = 1;
	}
	else  /* p_map.size() nonzero */ {

	  if (strainspecies) {
	    
	    set<uint32_t> write_set;
	    
	    for  (int count = 0; count < m_taxid_count; count ++) {

	      tid_T in_taxid;
	      m_table->next(m_offset, m_page, in_taxid);
	      m_taxid = in_taxid;
	    
	      if (p_map.find(m_taxid) == p_map.end() ) {
		write_set.insert(m_taxid);
	      } else {
		uint32_t mapped_tid = p_map[m_taxid];
		write_set.insert(mapped_tid);
	      
		
	      }
	    }
	    
	    // if not enough reduction , then we revert to the old method   
	    if (write_set.size() > tid_cut ) {
	      m_filtered_list.push_back(1);
	      m_taxid_count = 1;
	    } else /* write_set.size() <= tid_cut */ {
	      for (set<uint32_t>::const_iterator it = write_set.begin(); it != write_set.end(); it++)	      {
		m_filtered_list.push_back(*it);
	      }
	      
	      
	    }
	    m_taxid_count = m_filtered_list.size();
	    
	    
	  }  else /* not strain speces */ {
	    
      	    // lookup the rank value and put all taxids in priority queue
	  
	    for  (int count = 0; count < m_taxid_count; count ++) {


	      tid_T in_taxid;



	      m_table->next(m_offset, m_page, in_taxid);

	      /*	      if (in_taxid == BAD_TAXID)
		{
		  cout << "in_taxid (prefilter) " << BAD_TAXID <<  " , kmer = " << m_kmer << "\n"; 

		  } */


	      if (pp_map) {
		m_taxid = (*pp_map)[in_taxid];
		if (m_taxid == 0)
		  {
		    cout << "bad taxid: " <<  in_taxid << " kmer: " << m_kmer << " count: " << m_taxid_count << "\n";
		    assert(0);		    
		  }
	      }
	      else
		m_taxid = in_taxid;
	      
	      

	      const MyPair  pp(p_map[m_taxid], m_taxid);
	      taxid_q.push(pp);
	    }
	    
	    //	    m_filtered_list.clear();	  

	  // iterate on taxid priority queue in batches of taxons with
	  // the same rank

	    //	    cout << "kmer pruning ";
	    while(!taxid_q.empty()) {

	    // get current priorty
	      int cur_priority = taxid_q.top().first;
	      
	    // pull out elements that match top priority

	      //	      cout << " priority: " << cur_priority << " count " << taxid_q.size();
	      while(taxid_q.top().first == cur_priority) {
	    
		const MyPair res = taxid_q.top(); 
		taxid_q.pop();

		// FIX - don't add to qu
		//		m_filtered_list.push_back(res.second);
		if (taxid_q.empty())
		  break;
	      }


	    // if we are under the cut, then we can stop
	      if (taxid_q.size() <= tid_cut) {
		//	      cout << "Cut to rank: " << cur_priority << " org count  " << m_taxid_count << " new count" <<  m_filtered_list.size()  << "\n";
		m_taxid_count = taxid_q.size();
		//		cout << " priority: " << taxid_q.top().first << " count " << taxid_q.size();
		break;
		
	      } 
	      
	    }

	    cout <<"\n";
	    if (taxid_q.size() == 0) {
	      m_taxid_count = 1;
	      const MyPair pp(1, 1);
	      taxid_q.push(pp);
	      
	      //	    cout << "Cut to rank 1\n";
	    }
	  
	  }
	}
      }
    }

  }

  bool next() {
    
    if (m_calls_to_next >= m_taxid_count) {
      return false;
    }

    if (taxid_q.size() > 0) {

      const MyPair res = taxid_q.top(); 
      taxid_q.pop();

      m_taxid = res.second;

      //      m_taxid = m_filtered_list[m_calls_to_next];
      /*      if (m_taxid == BAD_TAXID)
	{
	  cout << "pruned taxid " << BAD_TAXID  <<  " , kmer = " << m_kmer << "\n"; 

	  } */
    }
    else if (mp_idmap) {
      // read from db then convert bit
      tid_T tid_16bit;
      m_table->next(m_offset, m_page, tid_16bit);

      m_taxid = (*mp_idmap)[tid_16bit];

      if (m_taxid == 0) {
	cout << "bad taxid: " << tid_16bit << " kmer: " << m_kmer << " count: " << m_taxid_count << "\n";
	assert(0);
      }
      
    } else {
      // read straight from db
      tid_T tid_in;
      m_table->next(m_offset, m_page, tid_in);

      /*      if (tid_in == BAD_TAXID)
	{
	  cout << "taxid "<< BAD_TAXID <<", kmer = " << m_kmer << "\n"; 

	  } */
      m_taxid = tid_in;
    }

    ++m_calls_to_next;
	
    return true;
  }

  uint32_t taxid() const {
    return m_taxid;
  }

  uint16_t taxidCount() const {
    return m_taxid_count;
  }

private:
  INDEXDB<tid_T> *m_table;
  page_t m_page;
  offset_t m_offset;
  uint16_t m_taxid_count;
  uint32_t m_taxid;
  uint64_t m_kmer;
  uint16_t m_calls_to_next;

  vector <tid_T> m_filtered_list;
  priority_queue<MyPair> taxid_q;

  id_convback_map_t *mp_idmap;

};



}
#endif
