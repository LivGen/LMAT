#ifndef __KMER_NODE_WITH_TAX_ID_HPP__
#define __KMER_NODE_WITH_TAX_ID_HPP__

# define UINT16_MAX   (65535)

#include <iostream>
#include <ostream>
#include <fstream>
#include <set>
#include <exception>
#include <cassert>
#include "perm.h"
#include "TaxTree.hpp"
#include "metag_typedefs.hpp"

namespace metag {

/**
 class for encapsulating information about a kmer, including
 taxonomy nodes associated with the kmer; this version does not
 implement genome positional information; this class is not thread
 safe, due to the static member s_work.
 */

template<class tid_T>
class KmerNode {
public:

  //! ctor
  KmerNode() : m_taxid_offset(0) {}

  //! dtor; needs to be redone to free memory from TaxidMetaData*
  ~KmerNode() { }

  //! returns the kmer encoded as a 64 bit integer
  kmer_t getKmer() const {
    return m_kmer;
  }

  //! insert a 64bit encoded kmer
  void setKmer(kmer_t kmer) {
    m_kmer = kmer;
  }

  char & getPage() {
    return m_taxid_page;
  }

  uint32_t & getOffset() {
    return m_taxid_offset;
  }

  const set<tid_T> & getTaxIDs() {
    return m_tids;
  }


  void read(FILE *fp) {
    uint32_t tid_ct, tid;
    m_tids.clear();
    assert(fread(&m_kmer, sizeof(kmer_t), 1, fp) == 1);
    assert(fread(&tid_ct, sizeof(uint32_t), 1, fp) == 1);

    for (uint32_t h=0; h<tid_ct; h++) {
      assert(fread(&tid, sizeof(uint32_t), 1, fp) == 1);
      m_tids.insert(tid);
    }
  }


  //! struct used for sorting a vector of KmerNodes
  struct NodeCmp {
    bool operator() (KmerNode* a, KmerNode* b) {
      return a->getKmer() < b->getKmer();
    }
  };


private:

  //! binary encoded kmer
  kmer_t m_kmer;

  uint32_t m_taxid_offset;
  char     m_taxid_page;

  std::set<tid_T> m_tids;
};

}

#endif
