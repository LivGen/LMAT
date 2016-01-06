#include <stdint.h>
#include <string>
#include <assert.h>

#ifndef KENCODE_H
#define KENCODE_H

using namespace std;


namespace kencode_ns {
  /**
Class to encode and decode 64-bit representation of a kmer string. 

For sliding window kmers of a genome or a read, it saves the 64-bit 
representation of the previous call,
and just shifts in the bit representation of the next character.
Therefore, the sliding window use is not thread-safe. Each thread
should have its own instance of the class.
The class is parameterized by the kmer length, the bit size of the
alphabet, and the encode/decode functions.
The kmer size * bit size of the alphabet  must be <= 64 bits. 
  */


  int tokbits_default(char c) {
    c = toupper(c);
    switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default: 
      cout << "character: '" << (int)c << "'\n" ; 
      //      assert(0);
      return 0;
    }
    
  }

  char tokchar_default(int k) {
    static char tab[4]= {'A', 'C', 'G', 'T'};
    assert (k>=0 || k<4);
    return tab[k];
  }

  typedef uint64_t  u64int;
  typedef int (*tokbits_t)(char c);
  typedef char (*tokchar_t)(int k);

class kencode_c {
private:
  u64int curr;
  int ksize;  // number of units in a kmer
  int asize;  // size of each alphabet unit in bits, eg. 2 for dna with 4 bases
  u64int kmask;  // mask of valid kmer bits - zero out unused high order bits
  int amask;     // mask of a unit (2 bits for the 4-base DNA)
  tokbits_t tokbits;
  tokchar_t tokchar;


public:
  inline kencode_c(int ks=6, int as=2, 
		   tokbits_t bp = &tokbits_default,
		   tokchar_t cp = &tokchar_default): 
    curr(0),ksize(ks), asize(as), tokbits(bp), tokchar(cp)
  {
    kmask = 0; amask = 0;
    for (int i=0;i<ksize*asize;i++)
      kmask = (kmask << 1) |1;
    for (int i=0;i<asize; i++)
      amask = (amask << 1) | 1;
  }


  uint64_t kencode(const char *str) {

    curr = 0;
    for (int i=0; i< ksize; i++)
      curr = (curr << asize) | tokbits(str[i]);
    return curr;
  }




  uint64_t kencode(const string &str) {
    assert (str.size()>=(unsigned) ksize);
    curr = 0;
    for (int i=0; i< ksize; i++)
      curr = (curr << asize) | tokbits(str[i]);
    return curr;
  }

  uint64_t kencode(char c) {
    curr = ((curr << asize) | tokbits(c)) & kmask;
    return curr;
  }
  
  void kdecode(uint64_t kid, string& s) {
    for (int i = ksize-1; i>=0; i--) {
      s[i]=tokchar(kid & amask);
      kid = kid >> asize;
    }
  };
};

}

#endif
