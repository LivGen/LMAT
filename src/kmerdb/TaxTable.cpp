#include "TaxTable.hpp"
#include "StopWatch.hpp"
#include "KmerFileMetaData.hpp"
#include <cassert>
#include <fstream>
#include <sstream>

#include <metag_const.h>

using namespace std;
using namespace metag;

//for debugging, so we can stop after reading in a portion of a file
int stop_early = 0;
int discard_if_more_than = 65535;


//version for 16-bit tax IDs, with optional remapping
#if 0
template<class tid_T>
void TaxTable<tid_T>::ingest(const char *tid_remap_fn) {
  cout << "starting TaxTable<tid_T>::ingest\n";
  if (tid_remap_fn) {
    if (!sizeof(tid_T) == 2) {
      cerr << "you must use uint16_t for tid_T, if you pass in a remapping file\n";
      exit(1);
    }
  }
  StopWatch clck;
  clck.start();
  size_t sanity_count = 0;
  uint64_t num_tid = 0;

  //read the tid mapping file
  map<uint32_t, tid_T> mp;
  if (tid_remap_fn) {
    ifstream in(tid_remap_fn);
    assert(in);
    uint32_t t1;
    tid_T t2;
    while (in >> t1 >> t2) {
      mp[t1 = t2];
    }
    in.close();
  }

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
    assert(metadata.version() == TAX_HISTO_VERSION);
    uint64_t kmer_ct = metadata.size();

    uint64_t sanity = ~0, test;

    kmer_t kmer;
    tid_T    tid;
    uint16_t tid_count;

    for (uint64_t i=0; i<kmer_ct; i++) {

      //loop exit condition
      if (ftell(in) == f) break;

      //sanity check
      if (i+1 % TAX_HISTO_SANITY_COUNT == 0) {
        assert(fread(&test, sizeof(uint64_t), 1, in) == 1);
        assert(test == sanity);
      }

      //read kmer and tid count
      assert(fread(&kmer, 8, 1, in) == 1);
      assert(fread(&tid_count, 2, 1, in) == 1);

      if (tid_count == 1) {
        ++num_tid;
        assert(fread(&tid, sizeof(tid_T), 1, in) == 1);

        if (!tid_remap_fn) {
          (*this)[kmer] = pair<tid_T, uint8_t>(tid , MAX_PAGE);
        } else {
          if (mp.find(tid) == mp.end()) {
            cout << "failed to find 32->16bit mapping for tax ID: " << tid << endl;
          } else {
            tid = (uint32_t)mp[tid];
            (*this)[kmer] = pair<tid_T, uint8_t>(tid , MAX_PAGE);
          }
        }
      }

      else {

        //add another page to storage?
        if (16+m_cur_offset+tid_count*(2+sizeof(tid_T)) > PAGE_SIZE) {
          addStorage();
        }

        //set entry in hash: kmer -> <offset, page>
        if (tid_count < discard_if_more_than) {
          (*this)[kmer] = pair<offset_t, page_t>(m_cur_offset, m_cur_page);
        }

        //sanity check
        if (kmer % TAX_TABLE_SANITY_COUNT == 0) {
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
          assert(fread(&tid, 4, 1, in) == 1);
          if (!tid_remap_fn) {
            (*this)[kmer] = pair<tid_T, uint8_t>(tid , MAX_PAGE);
          } else {
            if (mp.find(tid) == mp.end()) {
              cout << "failed to find 32->16bit mapping for tax ID: " << tid<< endl;
            } else {
              tid = mp[tid];
              (*this)[kmer] = pair<tid_T, uint8_t>(tid, MAX_PAGE);
            }
          }
        }

        mcpyin(tid, sizeof(tid_T));
        m_cur_offset += 4;
      }
    }
    if (j%10 == 9) {
      cout << j << " files read!\n";
    }

    fclose(in);
    cout << "ingest time: " << clck2.stop() << endl;
  }

  cout << "time to load TaxTable: " << clck.stop() << endl;
  cout << "tid count: " << num_tid << endl;
}
#endif

