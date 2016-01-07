#include "SortedDb.hpp"

#include <cassert>
#include <fstream>
#include <sstream>
#include <set>
#include <queue>

#include <metag_const.h>
#include "KmerFileMetaData.hpp"

#include <kencode.hpp>

using namespace std;
using namespace metag;
using namespace kencode_ns;

int metag::kmer_rec_comp(const void *a, const void *b)
{
    const kmer_record *kr_a = reinterpret_cast<const kmer_record*>(a);
    const kmer_record *kr_b = reinterpret_cast<const kmer_record*>(b);
    
    return kr_a->kmer_lsb - kr_b->kmer_lsb;
    
}

size_t ext_taxids = 0 ;
size_t singletons = 0;
size_t doubles = 0;

size_t reduced_kmers = 0;
size_t cut_kmers = 0;

size_t new_human = 0;
size_t matched_in = 0;
size_t new_isect = 0;


uint64_t read_encode(FILE *f, kencode_c &ken) {
    
    char buf[33];
    
    
    int rc = fscanf(f,"%s", buf);
    
    if (rc == EOF)
        return ~0;
    
    if (strlen(buf) > 0) {
        
        uint64_t kmer = ken.kencode(buf);
        //    cout << kmer << "\n";
        
        return kmer;
    }
    else
        return ~0;
    
}

kmer_set_t *get_kmer_set(FILE *in_kmers_fp, kencode_c &ken )
{
    
    uint64_t next_kmer = read_encode(in_kmers_fp, ken);
    
    kmer_set_t *kmer_set_tmp = new kmer_set_t;
    
    int i=0;
    
    while (next_kmer != (~0)) {
        kmer_set_tmp->insert(next_kmer);
        next_kmer = read_encode(in_kmers_fp, ken);
        i++;
    }
    
    cout << i << " adaptor k-mers \n";
    return kmer_set_tmp;
    
}


/* This should really be refactored!! */
template <class tid_T>
void SortedDb<tid_T>::add_data(const char *filename, size_t stopper, bool use_tax_histo_format, bitreduce_map_t *p_br_map, my_map &species_map, int tid_cutoff, bool strainspecies, FILE *human_kmers_fp, FILE * illum_kmers_fp, uint32_t adaptor_tid)
{
    
    
    static uint64_t last_kmer = 0;
    
    FILE *in = fopen(filename, "r");
    assert(in != NULL);
    //assert(fread(&kmer_count, 8, 1, fp) == 1);
    //for (uint64_t j=0; j<kmer_count; j++)
    fseek(in, 0, SEEK_END);
    long f = ftell(in);
    fseek(in, 0, SEEK_SET);
    
    KmerFileMetaData metadata;
    metadata.read(in);
    if (use_tax_histo_format)
        assert(metadata.version() == TAX_HISTO_VERSION);
    
    kencode_c ken(m_kmer_length);
    
    static uint64_t last_human = 0;
    
    if (human_kmers_fp)
        
        last_human=read_encode(human_kmers_fp, ken);
    else {
        last_human=~0;
    }
    
    static kmer_set_t *p_adaptor_set = NULL;
    
    if (!p_adaptor_set && illum_kmers_fp)
        
        p_adaptor_set = get_kmer_set(illum_kmers_fp, ken);
    
    
    uint64_t kmer_ct = metadata.size();
    
    uint64_t sanity = ~0, test;
    
    if (stopper == 0)
        stopper = ~0;
    
    uint64_t kmer;
    uint32_t tid;
    uint16_t tid_count; //, genome_count, p_count, tuple_count;
    
    uint32_t tid_count_32;
    
    ////uint16_t genome_count;
    
    static long long int start_count;
    static long long int start_offset;
    
    static uint16_t count_marker = 0;
    
    uint16_t two_count =  2;
    
    //  cout << "stopper set to: " << stopper << "\n";
    
    uint16_t HUMAN_16 = 0;
    uint16_t ADAPTOR_16 = 0;
    
    if (p_br_map) {
        HUMAN_16 = (*p_br_map)[9606];
        ADAPTOR_16 = (*p_br_map)[adaptor_tid];
    }
    
    for  ( uint64_t i=0; i<kmer_ct; i++)  {
        
        if (i > stopper)
            break;
        
        //loop exit condition
        if (ftell(in) == f) break;
        
        //read kmer, taxid count, and tuple count
        assert(fread(&kmer, 8, 1, in) == 1);
        
        if ((last_kmer > 0) &&  (kmer <= last_kmer)) {
            cout << "Kmers arriving out of order.  New: " << kmer << " last: " << last_kmer << "\n";
            exit(1);
        }
        
        
        while (last_human < kmer) {
            
            // This is a new human k_mer
            
            size_t top_index =  (last_human >> BITS_PER_2ND);  // & 0x0000000007ffffff;
            
            // check the slot if it has been written to yet
            
            if (top_tier_block[top_index] == 0) {
                
                
                
                // if empty write the current offset - this will be the start of the short list within the second tier
                start_offset = m_list_offset;
                start_count = 1;
            }
            
            uint16_t kmer_lsb_in = MASK_2ND & last_human;
            
            top_tier_block[top_index] = ((uint64_t) start_count << 48) | start_offset;
            
            kmer_table[m_list_offset].kmer_lsb = kmer_lsb_in;
            
            assert (top_tier_block[top_index] >> 48 == start_count);
            kmer_table[m_list_offset].page_id = MAX_PAGE;
            
            
            if (p_adaptor_set && p_adaptor_set->find(last_human) != p_adaptor_set->end()) {
                cout << "adaptor human k-mer " << last_human << "\n";
                
                if (ADAPTOR_16)
                    kmer_table[m_list_offset].page_offset = ADAPTOR_16;
                else
                    kmer_table[m_list_offset].page_offset = adaptor_tid;
                
                
            }
            
            
            else {
                if (HUMAN_16)
                    kmer_table[m_list_offset].page_offset = HUMAN_16;
                else
                    kmer_table[m_list_offset].page_offset = 9606;
            }
            new_human++;
            m_list_offset++;
            start_count++;
            
            
            last_human=read_encode(human_kmers_fp, ken);
            
        }
        
        bool add_human = false;
        
        if (last_human == kmer) {
            
            matched_in ++;
            
            add_human = true;
            
            last_human=read_encode(human_kmers_fp, ken);
        }
        
        
        if (use_tax_histo_format) {
            assert(fread(&tid_count, 2, 1, in) == 1);
        } else {
            assert(fread(&tid_count_32, 4, 1, in) == 1);
            
            tid_count = (uint16_t)tid_count_32;
            
        }
        
        // First get the mapping - the msb for the kmer - assuming 27 for now
        // TODO: make this configurable at runtime
        
        size_t top_index =  (kmer >> BITS_PER_2ND);  // & 0x0000000007ffffff;
        
        
        // check the slot if it has been written to yet
        
        if (top_tier_block[top_index] == 0) {
            
            
            
            // if empty write the current offset - this will be the start of the short list within the second tier
            start_offset = m_list_offset;
            start_count = 1;
        }
        
        uint16_t kmer_lsb_in = MASK_2ND & kmer;
        
        top_tier_block[top_index] = ((uint64_t) start_count << 48) | start_offset;
        
        kmer_table[m_list_offset].kmer_lsb = kmer_lsb_in;
        
        assert (top_tier_block[top_index] >> 48 == start_count);
        
        uint16_t tmp_tid_count = tid_count;
        
        //    set<uint32_t> write_set;
        priority_queue<MyPair> taxid_q;
        
        if (p_adaptor_set && p_adaptor_set->find(kmer) != p_adaptor_set->end()) {
            
            cout << "adaptor tax-histo k-mer " << kmer << "\n" ;
            for (uint16_t k=0; k<tid_count; k++)
                
                assert(fread(&tid, 4, 1, in) == 1);
            
            kmer_table[m_list_offset].page_id = MAX_PAGE;
            if (ADAPTOR_16)
                kmer_table[m_list_offset].page_offset = ADAPTOR_16;
            else
                kmer_table[m_list_offset].page_offset = adaptor_tid;
            
            m_list_offset++;
            start_count++;
            
            assert(start_count <= LENGTH_MAX_2ND );
            
        } else {
            
            
            if (tid_cutoff > 0 && tid_count > tid_cutoff) {
                
                if (species_map.size() == 0) {
                    tmp_tid_count = 0;
                    
                    for (uint16_t k=0; k<tid_count; k++) {
                        // scan through and dump
                        assert(fread(&tid, 4, 1, in) == 1);
                    }
                } else {
                    
                    if (strainspecies) {
                        
                        cout << "functionality disabled!\n";
                        exit(1);
                        
                        /*
                         for (uint16_t k=0; k<tid_count; k++) {
                         
                         
                         assert(fread(&tid, 4, 1, in) == 1);
                         
                         if (species_map.find(tid) == species_map.end()) {
                         // alternate version - try to not write the tid because
                         // doing so might be keeping the accuracy down
                         write_set.insert(tid);
                         
                         } else {
                         
                         uint32_t mapped_tid = species_map[tid];
                         write_set.insert(mapped_tid);
                         
                         }
                         }
                         
                         tmp_tid_count = write_set.size();
                         
                         // if not enough reduction , then we revert to the old method
                         if (tmp_tid_count > tid_cutoff)
                         tmp_tid_count = 0;
                         */
                    }
                    else {  // NOT STRAIN SPECIES
                        
                        
                        for (uint16_t k=0; k<tid_count; k++) {
                            
                            assert(fread(&tid, 4, 1, in) == 1);
                            
                            if (add_human && tid == 9606) {
                                add_human = false;
                            }
                            
                            const MyPair  pp(species_map[tid], tid);
                            taxid_q.push(pp);
                            
                            
                        }
                        
                        if (add_human) {
                            
                            
                            new_isect++;
                            const MyPair  pp(species_map[9606], 9606);
                            taxid_q.push(pp);
                            matched_in--;
                        }
                        
                        
                        // iterate on taxid priority queue in batches of taxons with
                        // the same rank
                        
                        while(!taxid_q.empty()) {
                            
                            // get current priorty
                            int cur_priority = taxid_q.top().first;
                            
                            // pull out elements that match top priority
                            while(taxid_q.top().first == cur_priority) {
                                
                                const MyPair res = taxid_q.top();
                                taxid_q.pop();
                                
                                //	      write_set.insert(res.second);
                                
                                if (taxid_q.empty())
                                    break;
                            }
                            // if we are under the cut, then we can stop
                            if (taxid_q.size() <= tid_cutoff) {
                                //            cout << "Cut to rank: " << cur_priority << " org
                                // cout  " << m_taxid_count << " new count" <<  m_filtered_list.size()  << "\n";
                                tmp_tid_count = taxid_q.size();
                                
                                break;
                                
                            }
                            //else {
                            // prepare for next batch of element poping
                            //write_set.clear();
                            //}
                            
                        }
                        if (taxid_q.size() == 0) {
                            tmp_tid_count = 1;
                            const MyPair pp(1, 1);
                            taxid_q.push(pp);
                            cut_kmers++;
                        }
                        
                    }
                }
                
            }
            
            
            if (tmp_tid_count > 1) {
                
                if (16+m_cur_offset+tmp_tid_count*4 > PAGE_SIZE) {
                    
                    m_cur_page ++;
                    m_cur_offset = 0;
                    
                }
                
                kmer_table[m_list_offset].page_id = m_cur_page;
                kmer_table[m_list_offset].page_offset = m_cur_offset;
            }
            else if (tid_count == 1) {
                
                assert(fread(&tid, 4, 1, in) == 1);
                
                
                // we need to handle the case where a single tid matches a human k-mer
                
                if (add_human && tid != 9606) {
                    
                    doubles++;
                    matched_in--;  // modify the counters
                    
                    if (m_cur_offset+32 > PAGE_SIZE) {
                        m_cur_page++;
                        m_cur_offset = 0;
                    }
                    
                    
                    
                    kmer_table[m_list_offset].page_id = m_cur_page;
                    kmer_table[m_list_offset].page_offset = m_cur_offset;
                    
                    if (kmer % 4096 == 0) {
                        mcpyinsdb(kmer, 8);
                        m_cur_offset += 8;
                        
                    }
                    
                    
                    
                    mcpyinsdb(two_count, 2);
                    m_cur_offset += 2;
                    
                    
                    // add the first tid from tax_histo input
                    
                    if (p_br_map) {
                        
                        uint16_t tid_16 = (*p_br_map)[tid];
                        if (tid_16 > p_br_map->size()+1 || tid_16 == 0) {
                            cout << "bad read: " << tid << " " << tid_16 << "\n";
                            assert(0);
                        }
                        
                        
                        //	  cout << "adding-tid " ;
                        
                        mcpyinsdb(tid_16, 2);
                        m_cur_offset += 2;
                        
                        
                        mcpyinsdb(HUMAN_16, 2);
                        m_cur_offset += 2;
                        
                        
                    } else {
                        
                        mcpyinsdb(tid, 4);
                        m_cur_offset += 4;
                        tid = 9606;
                        mcpyinsdb(tid, 4);
                        m_cur_offset += 4;
                        
                        
                    }
                    
                    
                    
                } else {
                    
                    
                    singletons++;
                    
                    
                    //	assert (tid <= MAX_TID && tid != INVALID_TID_2 );
                    
                    
                    kmer_table[m_list_offset].page_id = MAX_PAGE;
                    
                    if (p_br_map) {
                        uint16_t tid_16 = (*p_br_map)[tid];
                        if ((tid_16 > p_br_map->size()+1) || tid_16 ==0)  // pad for starting the count at 2
                        {
                            cout << "bad single: kmer " << kmer << " - "  << tid << " " << tid_16 << "\n";
                            assert(0);
                        }
                        kmer_table[m_list_offset].page_offset = tid_16;
                    }
                    else {
                        
                        kmer_table[m_list_offset].page_offset = tid;
                    }
                    
                }
            } 
            else if (tmp_tid_count == 1) {
                
                
                tid = taxid_q.top().second;
                kmer_table[m_list_offset].page_id = MAX_PAGE;
                
                if (p_br_map) {
                    uint16_t tid_16 = (*p_br_map)[tid];
                    if (tid_16 > p_br_map->size()+1 || tid_16 == 0)
                    {
                        cout << "bad single: " << tid << " " << tid_16 << "\n";
                        assert(0);
                    }
                    kmer_table[m_list_offset].page_offset = tid_16;
                } else {
                    kmer_table[m_list_offset].page_offset = tid;
                    
                }
                reduced_kmers++;
            } else {
                assert(tmp_tid_count == 0);      
                kmer_table[m_list_offset].page_id = MAX_PAGE;
                kmer_table[m_list_offset].page_offset = 1;
                cut_kmers++;
            }
            
            
            m_list_offset++;
            start_count++;
            
            assert(start_count <= LENGTH_MAX_2ND );
            
            //write the kmer
            //      memcpy(m_data[m_cur_page]+m_cur_offset, &kmer, 8);
            if (taxid_q.size() > 1 && tmp_tid_count > 1)	{
                
                // check for no reduction
                
                reduced_kmers++;
                if (kmer % 4096 == 0) {
                    mcpyinsdb(kmer, 8);
                    m_cur_offset += 8;
                    
                }	
                
                mcpyinsdb(tmp_tid_count, 2);
                m_cur_offset += 2;
                
                //      if (add_human) 
                //cout << "kmer-match: " << kmer;
                
                for (int i=0; i < tmp_tid_count; i++) {
                    
                    tid = taxid_q.top().second;
                    
                    //	if (add_human) 
                    //cout << " " << taxid_q.top().first << " " << tid;
                    
                    taxid_q.pop();
                    
                    if (p_br_map) {
                        
                        uint16_t tid_16 = (*p_br_map)[tid];
                        if (tid_16 > p_br_map->size()+1 || tid_16 == 0) {
                            cout << "bad set: " << tid << " " << tid_16 << "\n";
                            assert(0);
                        }
                        //	  assert (tid_16 > 0);
                        mcpyinsdb(tid_16,2);
                        m_cur_offset += 2;
                    }  else {
                        
                        mcpyinsdb(tid,4);
                        m_cur_offset += 4;
                    }
                    
                    ext_taxids++;
                } // end for 
                //      if (add_human) 
                //cout << "\n";
                
            }
            // no attempt to reduce list; copy in the taxid list
            else if (taxid_q.size() == 0 && tmp_tid_count > 1) {
                
                if (kmer % 4096 == 0) {
                    mcpyinsdb(kmer, 8);
                    m_cur_offset += 8;
                }
                
                
                uint32_t save_offset = m_cur_offset;
                
                mcpyinsdb(tid_count, 2);
                
                
                
                // save the offset in case we have a human match to something else      
                
                
                m_cur_offset += 2;
                
                char convbuf[12];
                
                string outbuf;
                
                //write the tuples
                for (uint16_t k=0; k<tid_count; k++) {
                    
                    ext_taxids++;
                    assert(fread(&tid, 4, 1, in) == 1);        
                    
                    if (tid == 9606)
                        add_human = false;
                    
                    if (add_human) {
                        
                        
                        sprintf(convbuf, " %d", tid);
                        outbuf += convbuf;
                    }
                    
                    if (p_br_map) {
                        
                        
                        uint16_t tid_16 = (*p_br_map)[tid];
                        if ( tid_16 > p_br_map->size()+1 || tid_16 == 0) {
                            cout << "bad read: " << tid << " " << tid_16 << "\n";
                            assert(0);
                        }
                        
                        
                        mcpyinsdb(tid_16, 2);
                        m_cur_offset += 2;
                    } else {
                        
                        mcpyinsdb(tid, 4);
                        m_cur_offset += 4;
                    }
                    
                }
                
                
                
                if (add_human)  {
                    
                    tid = 9606;
                    
                    if (p_br_map) {
                        
                        
                        uint16_t tid_16 = (*p_br_map)[tid];
                        if (tid_16 > p_br_map->size()+1 || tid_16 == 0) {
                            cout << "bad read: " << tid << " " << tid_16 << "\n";
                            assert(0);
                        }
                        
                        
                        
                        //	  cout << "adding-tid " ; 
                        
                        mcpyinsdb(tid_16, 2);
                        m_cur_offset += 2;
                    } else {
                        
                        mcpyinsdb(tid, 4);
                        m_cur_offset += 4;
                    }
                    
                    
                    // here we have to go back and increment the count
                    
                    uint32_t tmp_offset = m_cur_offset;
                    
                    m_cur_offset = save_offset;
                    
                    uint16_t inc_tid_count = tid_count + 1;
                    
                    
                    mcpyinsdb(inc_tid_count, 2);
                    
                    m_cur_offset = tmp_offset;
                    /*
                     cout << "kmer-match no pruning: " << kmer;
                     cout << " inc " << inc_tid_count;
                     cout << outbuf.c_str();
                     cout << "\n"; */
                    
                    new_isect++;
                }
                
                
            } else {
                assert (tid_count == 1 || tmp_tid_count < 2);
                
            }
        }
        
        if (use_tax_histo_format) {
            
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
        last_kmer = kmer;
    }
    
    m_n_kmers += kmer_ct;
    
    
    fclose(in);
    
    cout << "storage used for tax ids, counts, etc.."  << m_cur_page << " - "  << m_cur_offset << "\n";
    cout << "kmer count: " << m_n_kmers << "\n";
    cout << "offset at: " << m_list_offset << "\n";
    cout << "singletons: " << singletons << "\n";
    cout << "doubles: " << doubles << "\n";
    cout << "taxids in storage: " << ext_taxids << "\n";
    cout << "kmers reduced: " << reduced_kmers << "\n";
    cout << "kmers cut to 1: " << cut_kmers << "\n";
    cout << "new human k-mers: " << new_human << "\n";
    cout << "matched human k-mers: " << matched_in << "\n";
    cout << "new human + other k-mers: " << new_isect << "\n";
    
    
    
}



template class metag::SortedDb<DBTID_T>;

