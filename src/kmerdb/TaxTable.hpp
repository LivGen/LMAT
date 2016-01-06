#ifndef __TAX_TABLE__
#define __TAX_TABLE__

#include <ext/hash_map>
#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <stdint.h>
#include <cstdlib>
#include <string.h>
#include "TaxTree.hpp"
#include "StopWatch.hpp"
#include <metag_const.h>
#include "KmerFileMetaData.hpp"


#undef PAGE_SIZE
#define PAGE_SIZE 3221225472

#undef MAX_PAGE
#define MAX_PAGE 255  // would be 35535 for 16 bit field

#ifdef USE_BOOST

#include <boost/interprocess/offset_ptr.hpp>
#include <boost/interprocess/allocators/allocator.hpp>

#include <boost/unordered_map.hpp>



#define mcpyout(dest, page, off, len) memcpy(&dest, m_data[page].get()+off, len)
#define mcpyin(src, len) memcpy(reinterpret_cast<void*>(m_data[m_cur_page].get()+m_cur_offset), &src, len)


#else


#define mcpyout(dest, page, off, len) memcpy(&dest, m_data[page]+off, len)
#define mcpyin(src, len) memcpy(m_data[m_cur_page]+m_cur_offset, &src, len)
#endif

#include <perm.h>
namespace metag {


/**
A somewhat clumsy class for storing and retrieving
tax stat data
*/
//typedef uint64_t kmer_t;
typedef std::pair<uint32_t, uint8_t> table_val_t;


#if USE_BOOST == 1

namespace bip = boost::interprocess;

//typedef uint64_t kmer_t;
typedef std::pair<uint32_t, uint8_t> table_val_t;
typedef bip::allocator<std::pair<uint64_t, table_val_t>,  bip::managed_mapped_file::segment_manager>  kmer_pair_allocator_t;

typedef boost::unordered_map<kmer_t, table_val_t, boost::hash<kmer_t>, std::equal_to<kmer_t>,   kmer_pair_allocator_t>
TaxTableHash;

#elif defined _MMAP_MALLOC

typedef __gnu_cxx::hash_map<uint64_t, std::pair<uint32_t, uint8_t>, __gnu_cxx::hash<uint64_t>, __gnu_cxx::equal_to<uint64_t>,  PERM_NS::allocator<std::pair<uint32_t, uint8_t> > > TaxTableHash;

#elif WITH_PJMALLOC == 1

typedef __gnu_cxx::hash_map<uint64_t, std::pair<uint32_t, uint8_t>, __gnu_cxx::hash<uint64_t>, __gnu_cxx::equal_to<uint64_t>,  PERM_NS::allocator<std::pair<uint32_t, uint8_t> > > TaxTableHash;

#else
typedef __gnu_cxx::hash_map<uint64_t, std::pair<uint32_t, uint8_t> > TaxTableHash;
#endif

template<class tid_T>
class TaxTable : public TaxTableHash {
public :

  //! ctor

#if USE_BOOST == 1
  TaxTable(bip::managed_mapped_file const&in_heap) : TaxTableHash(kmer_pair_allocator_t(in_heap.get_segment_manager())),  m_cur_offset(1), m_cur_page(0), m_kmer_length(0) {
    m_heap = const_cast<bip::managed_mapped_file*>(&in_heap);
#else
  TaxTable() : m_cur_offset(0), m_cur_page(-1), m_kmer_length(0) {
#endif
    addStorage();
  }

#if USE_BOOST == 1
  TaxTable(bip::managed_mapped_file const&in_heap, size_t sz) : TaxTableHash(sz, boost::hash<kmer_t>(), std::equal_to<kmer_t>() , kmer_pair_allocator_t(in_heap.get_segment_manager())),  m_cur_offset(0), m_cur_page(-1), m_kmer_length(0) {
    m_heap = const_cast<bip::managed_mapped_file*>(&in_heap);
#else
  TaxTable(size_t sz) : TaxTableHash(sz), m_cur_offset(0), m_cur_page(-1), m_kmer_length(0)  {
#endif
    addStorage();
  }


#if USE_BOOST ==1
  TaxTable(kmer_pair_allocator_t const&alloc, char kl) : TaxTableHash(alloc),  m_cur_offset(0), m_cur_page(0), m_kmer_length(kl)
#else
  TaxTable(char kl) :  m_cur_offset(0), m_cur_page(0), m_kmer_length(kl)
#endif
  {
    addStorage();
  }

  //! dtor
  ~TaxTable() {}

  //! call once for each file being loaded; after registering
  //! all files call ingest();
  void registerFile(const char *fn) {
    std::cout << "registering file: " << fn << endl;
    m_filenames.push_back(fn);
  }



  void set_kmer_length(char c) {
    m_kmer_length = c;
  }

  char get_kmer_length() {
    return m_kmer_length;

  }


  //! Users should not call this method; instead, see TaxNodeStat
  //!
  //! on return count_out will
  //! contain the number of triples <tax_id, present>
  //! and genome_count_out will contain the total number of genomes
  //! in which the kmer appears (ignoring multiplicity).  If the kmer
  //! appears in more than 65536 genomes, returns false, in which case
  //! count_out and genome_count_out are meaningless, and you should
  //! not iterate over next(); if the kmer is not in the database, also
  //! return false.
  //!
  //! input: kmer; all other args are outputs
  bool begin_(kmer_t kmer_in, uint16_t &taxid_count_out,  offset_t &offset_out, page_t &page_out) {
    if ((*this).find(kmer_in) == (*this).end()) {
      return false;
    }

    offset_out = (*this)[kmer_in].first;
    page_out = (*this)[kmer_in].second;

    if (page_out == MAX_PAGE) {
      taxid_count_out = 1;
    } else {

      //sanity check
      uint64_t km;
      if (kmer_in % 4096 == 0) {
        mcpyout(km,page_out,offset_out, 8);
        assert(km == kmer_in);
        offset_out += 8;
      }

      mcpyout(taxid_count_out,page_out,offset_out,2);
      offset_out += 2;
    }
    return true;
  }

  //! get data for the next tuple <taxid, present, distance_to_root>;
  //! it is the user's responsibility to call next() not more than the value
  //! returned in count_out from the call to begin()
  void next(uint32_t &offset_out_in, uint8_t &page_in, tid_T &taxid_out) {
    if (page_in == MAX_PAGE) {
      taxid_out = offset_out_in;
    } else {
      mcpyout(taxid_out,page_in,offset_out_in,sizeof(tid_T));
      offset_out_in += sizeof(tid_T);
    }
  }


  //! load all data files
  //! tid_T should be uint16_t or uint32_t; if using uint16_t, then
  //! fn is a mapping file: uint32_t -> uint16_t
  void  ingest(bool use_tax_histo_files = true, bool use_16 = false) {
    cout << "starting TaxTable<tid_T>::ingest\n";
    cout << "use_tax_histo_files: " << use_tax_histo_files << endl;
    StopWatch clck;
    clck.start();
    size_t sanity_count = 0;
    uint64_t num_tid = 0;

    //loop over the files to be ingested
    for (size_t j=0; j<m_filenames.size(); j++) {
      StopWatch clck2;
      clck2.start();
      cout << "ingesting: " << m_filenames[j].c_str() << endl;

      FILE *in = fopen(m_filenames[j].c_str(), "r");
      assert(in != NULL);
      fseek(in, 0, SEEK_END);
      long f = ftell(in);
      fseek(in, 0, SEEK_SET);

      //sanity check
      KmerFileMetaData metadata;
      metadata.read(in);
      if (use_tax_histo_files) {
        assert(metadata.version() == TAX_HISTO_VERSION);
      }
      uint64_t kmer_ct = metadata.size();

      uint64_t sanity = ~0, test;

      kmer_t kmer;
      tid_T    tid;
      uint16_t tid_count;
      uint32_t tid_count_32;
      uint32_t tid_32;
      uint16_t tid_16;

      for (uint64_t i=0; i<kmer_ct; i++) {

        //loop exit condition
        if (ftell(in) == f) break;

        //read kmer and tid count
        assert(fread(&kmer, 8, 1, in) == 1);
        if (use_tax_histo_files) {
          assert(fread(&tid_count, 2, 1, in) == 1);
        } else {
          assert(fread(&tid_count_32, 4, 1, in) == 1);
          tid_count = (uint16_t)tid_count_32;
        }

        if (tid_count == 1) {
          ++num_tid;
          assert(fread(&tid_32, 4, 1, in) == 1);
          if (use_16) {
            tid_16 = (uint16_t)tid_32;
            (*this)[kmer] = pair<tid_T, uint8_t>(tid_16 , MAX_PAGE);
          } else {
            (*this)[kmer] = pair<tid_T, uint8_t>(tid_32 , MAX_PAGE);
          }
          (*this)[kmer] = pair<tid_T, uint8_t>(tid , MAX_PAGE);

        }


        else {

          //add another page to storage?
          if (16+m_cur_offset+tid_count*(2+sizeof(tid_T)) > PAGE_SIZE) {
            addStorage();
          }

          //set entry in hash: kmer -> <offset, page>
          //if (tid_count < discard_if_more_than) {
            (*this)[kmer] = pair<uint32_t, uint8_t>(m_cur_offset, m_cur_page);
          //}

          //sanity check
          if (kmer % 4096 == 0) {
            sanity_count ++;
            mcpyin(kmer, 8);
            m_cur_offset += 8;
          }

          //write taxid count
          mcpyin(tid_count, 2);
          m_cur_offset += 2;
          num_tid += tid_count;

          //write the tuples
          for (uint16_t k=0; k<tid_count; k++) {
            assert(fread(&tid_32, 4, 1, in) == 1);
            if (use_16) {
              tid_16 = (uint16_t)tid_32;
              mcpyin(tid, 2);
              m_cur_offset += 2;
            } else {  
              mcpyin(tid_32, 4);
              m_cur_offset += 4;
            }
          }

          //sanity check
          if (use_tax_histo_files) {
            if ((i+1) % TAX_HISTO_SANITY_COUNT == 0) {
              assert(fread(&test, sizeof(uint64_t), 1, in) == 1);
              assert(test == sanity);
            }
          } else {
            if ((i+1) % KMER_SANITY_COUNT == 0) {
              assert(fread(&test, sizeof(uint64_t), 1, in) == 1);
              assert(test == sanity);
            }
          }
        }
        }
        if ((j+1)%10 == 0) {
          cout << j << " files read!\n";
        }

        fclose(in);
        cout << "ingest time: " << clck2.stop() << endl;
      }

      cout << "time to load TaxTable: " << clck.stop() << endl;
      cout << "tid count: " << num_tid << endl;
    }

private:

#if USE_BOOST == 1
    bip::managed_mapped_file *m_heap;
    std::vector<std::string> m_filenames;
    bip::offset_ptr<char> m_data[256];

#elif defined _MMAP_MALLOC
    char *m_data[256];
    int m_heap_handle;

#elif WITH_PJMALLOC == 1
    std::vector<char*, PERM_NS::allocator<char*> > m_data;
    std::vector<std::string> m_filenames;
#else
    std::vector<char*> m_data;
    std::vector<std::string> m_filenames;
#endif

    //std::vector<std::vector<char> > m_data;
    offset_t m_cur_offset;
    page_t m_cur_page;

    char m_kmer_length;

    void addStorage() {

      ++m_cur_page;

#if USE_BOOST == 1

      m_data[m_cur_page] = new(m_heap->allocate(PAGE_SIZE*sizeof(char))) char[PAGE_SIZE];

#elif WITH_PJMALLOC == 1
      m_data.resize(m_data.size()+1);
      m_data.back() = new(JEMALLOC_P(malloc)(PAGE_SIZE*sizeof(char))) char[PAGE_SIZE];
#else
      m_data.resize(m_data.size()+1);
      m_data.back() = new char[PAGE_SIZE];
#endif

      std::cout << "TaxTable::addStorage; current number of pages: " << m_cur_page+1 << std::endl;

      m_cur_offset = 0;
    }


  };

}

#endif
