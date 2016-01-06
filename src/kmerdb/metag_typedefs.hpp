#ifndef __METAG_TYPEDEF__
#define __METAG_TYPEDEF__

#define TAX_HISTO_VERSION 999

#define KMER_SANITY_COUNT 1000

//write a sanity count after every TAX_HISTO_SANITY_COUNT entries in a tax_histo file
#define TAX_HISTO_SANITY_COUNT 1500

#define TAX_TABLE_SANITY_COUNT 4096

typedef uint64_t kmer_t;
typedef uint32_t offset_t;
typedef uint8_t page_t;

#endif
