#ifndef __KMER_ITERATOR__
#define __KMER_ITERATOR__

#include "KmerFileMetaData.hpp"
#include "metag_typedefs.hpp"
#include "KmerNode.hpp"
#include <set>
#include <cassert>
#include <cstdio>
#include <stdint.h>

namespace metag {

/* I always seem to be writing the same code to iterate over the
 * kmerDB files (the output from jellylist/kmerPrefixCounter),
 * so it makes sense to but this in a class
 */

template<class tid_T>
class KmerIterator {
public :
  KmerIterator(const char *fn) : m_num_read(0), m_sanity(~0) {
    m_fp = fopen(fn, "r");
    assert(m_fp);
    m_metadata.read(m_fp); 
    m_kmer_count = m_metadata.size();
  }

  uint64_t size() {
    return m_kmer_count;
  }

  //! returns the kmer 
  bool next(kmer_t &kmer) {
    if (!getNext()) {
      return false;
    }
    kmer = m_kmer;
    return true;
  }  

  //! returns the kmer and the set of tax IDs
  bool next(kmer_t &kmer, const std::set<tid_T> *tax_ids) {
    if (!getNext()) {
      return false;
    }
    kmer = m_kmer;
    tax_ids = m_taxids;
    return true;
  }  

  //! returns the tax ID count
  bool next(uint16_t &tid_count) {
    if (!getNext()) {
      return false;
    }
    tid_count = m_taxids->size();
    return true;
  }  
   
private :
  KmerFileMetaData m_metadata;
  FILE *m_fp;
  uint64_t m_num_read;
  uint64_t m_kmer_count;
  uint64_t m_sanity, m_test;
  KmerNode<tid_T> m_node;
  kmer_t m_kmer;
  const std::set<tid_T> *m_taxids;

  bool getNext() {
    if (! m_fp) {
      return false;
    }
    m_node.read(m_fp);
    ++m_num_read;
    if (m_num_read == m_kmer_count) {
      fclose(m_fp);
      m_fp = 0;
    } else if (m_num_read % 1000 == 0) {
      assert(fread(&m_test, sizeof(uint64_t), 1, m_fp) == 1);
      assert(m_test == m_sanity);
    }
    m_kmer = m_node.getKmer();
    m_taxids = &m_node.getTaxIDs();
    return true;
  }
};


}

#endif
