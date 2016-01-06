#ifndef __SORTED_DB_HPP__
#define __SORTED_DB_HPP__



#include <iostream>
#include <set>
#include <ext/hash_map>
#include <ext/hash_set>

#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <stdint.h>
#include <cstdlib>
#include <string.h>


#include <metag_const.h>
#include "metag_typedefs.hpp"


#if (WITH_PJMALLOC == 0)
#define JEMALLOC_P(x) x
#endif


#undef PAGE_SIZE 
#undef MAX_PAGE

#define PAGE_SIZE 4294701056
#define MAX_PAGE 255  // would be 35535 for 16 bit field



// 18-mer



// 20-mer

#define TT_BLOCK_COUNT_18 134217728
#define BITS_PER_2ND_18 9
#define MASK_2ND_18 0x00000000000001ff
#define LENGTH_MAX_2ND_18 513

#define TT_BLOCK_COUNT_20  134217728
#define BITS_PER_2ND_20 13
#define MASK_2ND_20 0x0000000000001fff
#define LENGTH_MAX_2ND_20 8193


#if (IDX_CONFIG == 2031)


#define TT_BLOCK_COUNT 2147483648 
#define BITS_PER_2ND 9
#define MASK_2ND 0x00000000000001ff
#define LENGTH_MAX_2ND 513

#elif (IDX_CONFIG == 2030)

#define TT_BLOCK_COUNT 1073741824 
#define BITS_PER_2ND 10
#define MASK_2ND 0x00000000000003ff
#define LENGTH_MAX_2ND 1025


#elif (IDX_CONFIG == 2029)

#define TT_BLOCK_COUNT 536870912 
#define BITS_PER_2ND 11
#define MASK_2ND 0x00000000000007ff
#define LENGTH_MAX_2ND 2049

#elif (IDX_CONFIG == 2028)

#define TT_BLOCK_COUNT 268435456
#define BITS_PER_2ND 12
#define MASK_2ND 0x0000000000000fff
#define LENGTH_MAX_2ND 4097

#elif (IDX_CONFIG == 2027)

#define TT_BLOCK_COUNT  134217728
#define BITS_PER_2ND 13
#define MASK_2ND 0x0000000000001fff
#define LENGTH_MAX_2ND 8193

#elif (IDX_CONFIG == 2026)

#define TT_BLOCK_COUNT  67108864
#define BITS_PER_2ND 14
#define MASK_2ND 0x0000000000003fff
#define LENGTH_MAX_2ND 16385

#elif (IDX_CONFIG == 2025)

#define TT_BLOCK_COUNT  33554432
#define BITS_PER_2ND 15
#define MASK_2ND 0x0000000000007fff
#define LENGTH_MAX_2ND 32769


#elif (IDX_CONFIG == 2024)

#define TT_BLOCK_COUNT  16777216
#define BITS_PER_2ND 16
#define MASK_2ND 0x000000000000ffff
#define LENGTH_MAX_2ND 65537


#elif (IDX_CONFIG == 1827)

#define TT_BLOCK_COUNT 134217728
#define BITS_PER_2ND 9
#define MASK_2ND 0x00000000000001ff
#define LENGTH_MAX_2ND 513

#endif

#include <perm.h>


typedef  __gnu_cxx::hash_map<uint32_t,uint32_t> my_map;
typedef  __gnu_cxx::hash_map<uint32_t,uint16_t> bitreduce_map_t;
typedef   __gnu_cxx::hash_set<uint64_t> kmer_set_t;

#define mcpyoutsdb(dest, page, off, len) memcpy(&dest, m_storage_space+(PAGE_SIZE*page)+off, len) 
#define mcpyinsdb(src, len) memcpy(m_storage_space+(m_cur_page*PAGE_SIZE)+m_cur_offset, &src, len) 

#define mcpyinsdbt(src, len, i) memcpy(m_storage_space+(m_cur_page_arr[i]*PAGE_SIZE)+m_cur_offset_arr[i], &src, len) 



namespace metag {

class MyPair
{
public:
  MyPair(unsigned int f, uint32_t s) : first(f), second(s) {}

  const bool operator < (const MyPair &mp ) const {
    return (first < mp.first);
  }

  unsigned int first;
  uint32_t second;
};



typedef struct {
  uint16_t kmer_lsb;
  uint16_t page_id;
  uint32_t page_offset;
	  
} kmer_record;

  int kmer_rec_comp(const void *a, const void *b);


  template<class tid_T>
  class SortedDb {

public:



    SortedDb(size_t n_kmers, size_t space_size) 

  {

    top_tier_block = new (JEMALLOC_P(malloc)(sizeof(uint64_t)*TT_BLOCK_COUNT)) uint64_t[TT_BLOCK_COUNT];
    kmer_table = new (JEMALLOC_P(malloc)(sizeof(kmer_record)*n_kmers)) kmer_record[n_kmers];
    m_storage_space = new (JEMALLOC_P(malloc)(sizeof(char)*space_size)) char[space_size];
    assert (m_storage_space);

    m_n_kmers = 0;

    m_cur_offset = 0;
    m_cur_page = 0;
    m_list_offset = 0;

    bzero((void*)top_tier_block, sizeof(uint64_t) *TT_BLOCK_COUNT);
    
    std::cout << " init db with " << n_kmers << " kmers and " << space_size << " for storage\n" ;

    
  }

    


    void add_data(const char *, size_t, bool, bitreduce_map_t *, my_map &, int, bool, FILE *, FILE *, uint32_t);


    bool begin_(uint64_t kmer_in, uint16_t &taxid_count_out,  uint32_t &offset_out, uint8_t &page_out) {
      
      switch ( m_kmer_length) {
      case 20:
	return    begin_20(kmer_in, taxid_count_out,  offset_out, page_out);
      case 18:	
	return    begin_18(kmer_in, taxid_count_out,  offset_out, page_out);
      default:
	std::cerr << "K size " << m_kmer_length << " not supported by this application version!\n";
	assert( 0);
      }
      
    }
    
    bool begin_18(uint64_t kmer_in, uint16_t &taxid_count_out,  uint32_t &offset_out, uint8_t &page_out) {

      // perform first level lookup
      size_t top_index =  (kmer_in >> BITS_PER_2ND_18) ; //  & 0x0000000007ffffff;


      uint64_t tt_val = top_tier_block[top_index];
      // no kmers present for that prefix
      if (tt_val == 0)
	return false;
	
      // there are kmers present so we now search in the second level
      uint16_t k_count =   (tt_val >> 48) ;
      long long int kmer_offset = tt_val & 0x0000ffffffffffff;

      // if this fails there probably is db corruption

      assert (k_count < LENGTH_MAX_2ND_18);

      if (k_count == 1) {
	// simple case of one kmer - no search needed to but need to
	// check the suffix bits, bail if not matching
	if (kmer_table[kmer_offset].kmer_lsb != (kmer_in & MASK_2ND_18))
	  return false;

	offset_out = kmer_table[kmer_offset].page_offset;
	page_out = kmer_table[kmer_offset].page_id;
      }
      // this case we search the array and bail if it is not found
      else {

	kmer_record rec_in;
	rec_in.kmer_lsb=  kmer_in & MASK_2ND_18; 

	kmer_record *res = reinterpret_cast<kmer_record*>(bsearch(&rec_in, kmer_table + kmer_offset ,  k_count, sizeof(kmer_record), kmer_rec_comp));
	if (!res)
	  return false;

	offset_out = res->page_offset;
	page_out = res->page_id;
	
      }
      
      // we found the kmer so now we can move onto the taxid list info
      // if we got this far we return true at the end
      
      // special case for one taxid
      if (page_out == MAX_PAGE) {
	taxid_count_out =1;
	
      } else {
	
	if (kmer_in % 4096 == 0) {
	  uint64_t km;
	  mcpyoutsdb(km, page_out,offset_out, 8);
        //    memcpy(&km, m_data[page_out]+offset_out, 8);
	  assert(km == kmer_in);
	  offset_out += 8;
	}
      
	//    memcpy(&taxid_count_out, m_data[page_out]+offset_out, 2);
	//	mcpyoutsdb(taxid_count_out,page_out,offset_out,2);
	//	offset_out += 2;
	//    memcpy(&genome_count_out, m_data[page_out]+offset_out, 2);
	//mcpyoutsdb(genome_count_out,page_out,offset_out,2);
	//offset_out += 2;
	//    memcpy(&tuple_count_out, m_data[page_out]+offset_out, 2);
	mcpyoutsdb(taxid_count_out,page_out,offset_out,2);
	offset_out += 2;



      }
      return true;

    }

bool begin_20(uint64_t kmer_in, uint16_t &taxid_count_out,  uint32_t &offset_out, uint8_t &page_out) {

      // perform first level lookup
  size_t top_index =  (kmer_in >> BITS_PER_2ND_20) ; //  & 0x0000000007ffffff;


      uint64_t tt_val = top_tier_block[top_index];
      // no kmers present for that prefix
      if (tt_val == 0)
	return false;
	
      // there are kmers present so we now search in the second level
      uint16_t k_count =   (tt_val >> 48) ;
      long long int kmer_offset = tt_val & 0x0000ffffffffffff;

      // if this fails there probably is db corruption

      assert (k_count < LENGTH_MAX_2ND_20);

      if (k_count == 1) {
	// simple case of one kmer - no search needed to but need to
	// check the suffix bits, bail if not matching
	if (kmer_table[kmer_offset].kmer_lsb != (kmer_in & MASK_2ND_20))
	  return false;

	offset_out = kmer_table[kmer_offset].page_offset;
	page_out = kmer_table[kmer_offset].page_id;
      }
      // this case we search the array and bail if it is not found
      else {

	kmer_record rec_in;
	rec_in.kmer_lsb=  kmer_in & MASK_2ND_20; 

	kmer_record *res = reinterpret_cast<kmer_record*>(bsearch(&rec_in, kmer_table + kmer_offset ,  k_count, sizeof(kmer_record), kmer_rec_comp));
	if (!res)
	  return false;

	offset_out = res->page_offset;
	page_out = res->page_id;
	
      }
      
      // we found the kmer so now we can move onto the taxid list info
      // if we got this far we return true at the end
      
      // special case for one taxid
      if (page_out == MAX_PAGE) {
	taxid_count_out =1;
	
      } else {
	
	if (kmer_in % 4096 == 0) {
	  uint64_t km;
	  mcpyoutsdb(km, page_out,offset_out, 8);
        //    memcpy(&km, m_data[page_out]+offset_out, 8);
	  assert(km == kmer_in);
	  offset_out += 8;
	}
      
	//    memcpy(&taxid_count_out, m_data[page_out]+offset_out, 2);
	//	mcpyoutsdb(taxid_count_out,page_out,offset_out,2);
	//	offset_out += 2;
	//    memcpy(&genome_count_out, m_data[page_out]+offset_out, 2);
	//mcpyoutsdb(genome_count_out,page_out,offset_out,2);
	//offset_out += 2;
	//    memcpy(&tuple_count_out, m_data[page_out]+offset_out, 2);
	mcpyoutsdb(taxid_count_out,page_out,offset_out,2);
	offset_out += 2;



      }
      return true;

    }
    /*
    bool exists(uint64_t kmer_in) {
      size_t top_index =  (kmer_in >> BITS_PER_2ND) ; //  & 0x0000000007ffffff;
      uint64_t tt_val = top_tier_block[top_index];
      if (tt_val == 0) {
        return false;
      }
      return true;
    }
    */	
	
   void next(uint32_t &offset_out_in, uint8_t &page_in, tid_T &taxid_out) {

     

    if (page_in == MAX_PAGE) {
      taxid_out = offset_out_in;

    } else {

    //    memcpy(&taxid_out, m_data[page_in]+offset_out_in, 4);
      mcpyoutsdb(taxid_out,page_in,offset_out_in,sizeof(tid_T));
      //      assert (taxid_out <= MAX_TID && taxid_out != INVALID_TID_2);
      
      offset_out_in += sizeof(tid_T);
    //    memcpy(&present_out, m_data[page_in]+offset_out_in, 2);
    //  mcpyoutsdb(present_out,page_in,offset_out_in,2);
      // offset_out_in += 2;
    }

   }


   void get_values(uint32_t *in_arr, uint32_t size)

   {

     int factor = TT_BLOCK_COUNT / size;

     size_t i;

     for (i = 0 ; i < TT_BLOCK_COUNT; i++) {
       int idx = i / factor;

       if (i % size == 0) 
	 in_arr[idx] = 0;
       
       in_arr[idx] += (top_tier_block[i] >> 48);


     }
     

   }

#ifdef USE_BOOST
  void conv_ptrs() 
  {
    //    if (strcmp(boost_version, BOOST_LIB_VERSION) != 0)
    //  std::cout << "Warning! This application built with boos libraries v. " << BOOST_LIB_VERSION << ", database built with " << boost_version << ". No guarantees of compatibility." ;


    top_tier_block = m_top_tier_block_op.get();
    kmer_table = m_kmer_table_op.get();
    m_storage_space = m_storage_space_op.get();

  }
#endif

  bool check_config() {
    return  (idx_config == IDX_CONFIG);
  }
  void set_kmer_length(char c)
  {
    m_kmer_length = c;
    idx_config = IDX_CONFIG;
  }

  char get_kmer_length() {
    return m_kmer_length;
    
  }
         
  size_t size() {
    return m_list_offset;
  }

private:

	   //  char *allocator_str;
	   //  char *lmat_version;

#ifdef USE_BOOST
  //  char *boost_version;
  //#elif (WITH_PJMALLOC == 1) 
  //  char *perm_je_version;
#endif

  int idx_config;

  size_t m_n_kmers;

  uint8_t m_kmer_length;

#ifdef USE_BOOST  
  bip::offset_ptr<char> m_storage_space_op;
  bip::offset_ptr<kmer_record> m_kmer_table_op;
  bip::offset_ptr<uint64_t> m_top_tier_block_op;
#endif

  char *m_storage_space;
  kmer_record *kmer_table; 
  uint64_t *top_tier_block;


  // used for taxid lists



  size_t  m_list_offset ;
  uint16_t m_cur_page;
  uint32_t m_cur_offset;


  size_t * list_offset_arr ;
  uint16_t * m_cur_page_arr;
  uint32_t * m_cur_offset_arr;


  };

  }
#endif
