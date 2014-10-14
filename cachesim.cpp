#include "cachesim.hpp"
#include <iostream>
using namespace std;
//Global Variable Declaration
uint64_t **tag_store;
bool **dirty;
uint64_t **age;      //lru tracker
bool **p_bit;
bool **valid;

int offset_size;
int index_size;
int tag_size;
int asso;       //2^s

uint64_t *vc_tag;
uint64_t *vc_index;
bool *vc_dirty;
int *vc_age;
int vc_size;    //v
bool *vc_p_bit;
bool *vc_valid;

int p_degree;   //k
uint64_t LastMissBlockAddress = 0;
uint64_t PendingStride = 0;

//Debug
uint64_t temp1 = 0;


void init_stats(cache_stats_t *stats, uint64_t s) {
	stats->accesses = 0;
    stats->reads = 0;
    stats->read_misses = 0;
    stats->read_misses_combined = 0;
    stats->writes = 0;
    stats->write_misses = 0;
    stats->write_misses_combined = 0;
    stats->misses = 0;
	stats->write_backs = 0;
	stats->vc_misses = 0;
	stats->prefetched_blocks = 0;
	stats->useful_prefetches = 0;
	stats->bytes_transferred = 0;

	stats->hit_time = 2 + (float)s/5;
    stats->miss_penalty = 200;
    stats->miss_rate = 0;
    stats->avg_access_time = 0;
}
/**
 * Subroutine for initializing the cache. You many add and initialize any global or heap
 * variables as needed.
 * XXX: You're responsible for completing this routine
 *
 * @c The total number of bytes for data storage is 2^C
 * @b The size of a single cache line in bytes is 2^B
 * @s The number of blocks in each set is 2^S
 * @v The number of blocks in the victim cache is V
 * @k The prefetch distance is K
 */
void setup_cache(uint64_t c, uint64_t b, uint64_t s, uint64_t v, uint64_t k) {

	tag_size = 64-c+s;
	index_size = c-b-s;
	offset_size = b;
	asso = (1 << s);
	vc_size = v;
	p_degree = k;
	const int sets = asso;
	const int total_lines = 1 << (index_size);			// number of rows
    const int vc_lines = vc_size;
	// Set Cache Rows
	tag_store = new uint64_t*[total_lines];
	dirty = new bool*[total_lines];
	age = new uint64_t*[total_lines];
	p_bit = new bool*[total_lines];
	valid = new bool*[total_lines];
	// Set Cache Columns
	for (int i=0; i < total_lines; i++) {
		tag_store[i] = new uint64_t[sets];
		dirty[i] = new bool[sets];
		age[i] = new uint64_t[sets];
		p_bit[i] = new bool[sets];
		valid[i] = new bool[sets];

		for(int j=0; j<sets; j++) {
            dirty[i][j] = false;
            age[i][j] = 0;
            p_bit[i][j] = false;
            valid[i][j] = false;
		}
	}

	//set victom cache blocks and inititalize dirty bits to 0
	vc_tag = new uint64_t[vc_lines];
	vc_index = new uint64_t[vc_lines];
	vc_dirty = new bool[vc_lines];
	vc_age = new int[vc_lines];
	vc_p_bit = new bool[vc_lines];
	vc_valid = new bool[vc_lines];
	for(int i=0; i<vc_lines; i++) {
        vc_age[i] = 0;
        vc_tag[i] = 0;
        vc_index[i] = 0;
        vc_dirty[i] = false;
        vc_p_bit[i] = false;
        vc_valid[i] = false;
    }
}
/**
 * Subroutine that simulates the cache one trace event at a time.
 * XXX: You're responsible for completing this routine
 *
 * @rw The type of event. Either READ or WRITE
 * @address  The target memory address
 * @p_stats Pointer to the statistics structure
 */
void cache_access(char rw, uint64_t address, cache_stats_t* p_stats) {

	p_stats->accesses++;
	uint64_t tag, index, offset, d;
	bool hit, vc_hit;
	int hit_set;
	hit = false;
	vc_hit = false;
	//separate offset, index and tag
	offset = (address & ((1 << offset_size) - 1));
	index = (address & (((1 << index_size) - 1) << offset_size)) >> offset_size;
	tag = (address & (((1 << tag_size) - 1)<< (offset_size + index_size))) >> (offset_size+index_size);
	//cout << hex << address << "\t" << hex << tag << "\t" << hex << index << "\t" << hex << offset << endl;

    //Check tag store for tag
    for (int i = 0; i< asso; i++) {
        if (tag == tag_store[index][i] && valid[index][i]) {
            hit = true;
            hit_set = i;
            break;
        }
    }

    if(hit) {
        // IT'S A HIT
        // update block to MRU
        age[index][hit_set] = p_stats->accesses;
        // Check if block was prefetched
        if(p_bit[index][hit_set]) {
			p_stats->useful_prefetches++;
			p_bit[index][hit_set] = false;
		}
		// Update Stats
        if(rw == WRITE) {
            dirty[index][hit_set] = true;
            p_stats->writes++;
        }
        else { p_stats->reads++; }
    }
    else {
        //IT'S A MISS!
        p_stats->misses++;
        //Find lru block to be replaced
        uint64_t lru_age = p_stats->accesses;
        int lru_block;
        for(int i=0; i<asso; i++) {
            //Check for invalidity
            if(!valid[index][i]) {
                lru_block = i;
                break;
            }
            if(age[index][i] <= lru_age) {
                lru_age = age[index][i];
                lru_block = i;
            }
        }
        if(vc_size > 0) {    //if vc exists
        //Victim cache Exists
            //Check VC for hit
            int vc_block;
            for(int i=0; i<vc_size; i++) {
                if(vc_valid[i] && index == vc_index[i] && tag == vc_tag[i]) {
                    vc_hit = true;
                    vc_block = i;
                    break;
                }
            }
            if(vc_hit) {
            //If VC hits
                // Check is block was prefetched
                if( vc_p_bit[vc_block]) {
                    p_stats->useful_prefetches++;
                    vc_p_bit[vc_block] = false;
                }
                uint64_t temp_tag = tag_store[index][lru_block];
                bool temp_dirty = dirty[index][lru_block];
                bool temp_p_bit = p_bit[index][lru_block];
                bool temp_valid = valid[index][lru_block];

                tag_store[index][lru_block] = vc_tag[vc_block];
                dirty[index][lru_block] = vc_dirty[vc_block];
                p_bit[index][lru_block] = vc_p_bit[vc_block];
                valid[index][lru_block] = vc_valid[vc_block];

                vc_tag[vc_block] = temp_tag;
                vc_dirty[vc_block] = temp_dirty;
                vc_index[vc_block] = index;
                vc_p_bit[vc_block] = temp_p_bit;
                vc_valid[vc_block] = temp_valid;

                // update block to MRU
                age[index][lru_block] = p_stats->accesses;

                // update dirty bit if write access
                // Update Stats
                if(rw == WRITE){
                    dirty[index][lru_block] = true;
                    p_stats->writes++;
                    p_stats->write_misses++;
                }
                else {
                    p_stats->reads++;
                    p_stats->read_misses++;
                }
            }
            else {
            //VC misses
                p_stats->vc_misses++;
                uint64_t min_acc = 0;
                //Find block to be replaced in vc (first_in block will be evicted)
                int first_in, max_age = 0;
                for(int i = 0; i<vc_size;i++) {
                    if(!vc_valid[i]) {
                        first_in = i;
                        break;
                    }
                    if(vc_age[i] >= max_age) {
                        max_age = vc_age[i];    //max_age => first in
                        first_in = i;
                    }
                }
                // Evict block in VC and writeback data if dirty
                if(vc_valid[first_in] && vc_dirty[first_in]) {
                    p_stats->write_backs++;
                }
                // Replace block in victim cache with cache lru block
                vc_tag[first_in] = tag_store[index][lru_block];
                vc_index[first_in] = index;
                vc_dirty[first_in] = dirty[index][lru_block];
                vc_p_bit[first_in] = p_bit[index][lru_block];
                vc_valid[first_in] = valid[index][lru_block];

                // Update vc FIFO
                for(int i = 0; i<vc_size;i++){
                    vc_age[i]++;
                }
                vc_age[first_in] = 0;

                // Put new block in cache
                tag_store[index][lru_block] = tag;
                dirty[index][lru_block] = false;
                p_bit[index][lru_block] = false;
                valid[index][lru_block] = true;

                // Update cache block to MRU
                age[index][lru_block] = p_stats->accesses;

                // Update Stats
                if(rw == WRITE) {
                    dirty[index][lru_block] = true;    //Due to write allocate
                    p_stats->writes++;
                    p_stats->write_misses++;
                    p_stats->write_misses_combined++;
                }
                else {
                    p_stats->reads++;
                    p_stats->read_misses++;
                    p_stats->read_misses_combined++;
                }
            } //endelse vc_misses
            // Check prefetcher
            if(p_degree!=0) {
            //Is prefetching defined?
                uint64_t now = (address & ~((1<<offset_size) -1));
                uint64_t p_offset,p_index,p_tag;
                int block;
                d = now - LastMissBlockAddress;
                LastMissBlockAddress = now;

                if(d == PendingStride) {
                    for(int i=1; i<=p_degree; i++) {
                        p_stats->prefetched_blocks++;
                        uint64_t p_addr;
                        bool p_vc_hit = false;
                        bool p_hit = false;
                        p_addr = now + (i*PendingStride);
                        p_offset = (p_addr & ((1 << offset_size) - 1));
                        p_index = (p_addr & (((1 << index_size) - 1) << offset_size)) >> offset_size;
                        p_tag = (p_addr & (((1 << tag_size) - 1)<< (offset_size + index_size))) >> (offset_size+index_size);
                        uint64_t min_acc = p_stats->accesses;
                        // Check for prefetched block in cache
                        for( int j=0; j<asso; j++) {
                            if( valid[p_index][j] && p_tag == tag_store[p_index][j] ) {
                                p_hit = true;
                                break;
                            }
                        }
                        // If block present, then do nothing
                        if(!p_hit) {
                        // if prefetched block is not in cache
                            // find lru block to be replaced
                            min_acc = p_stats->accesses;
                            for(int j=0; j<asso; j++) {
                                if(!valid[p_index][j]) {
                                    block = j;
                                    break;
                                }
                                if(age[p_index][j]<=min_acc) {
                                    min_acc = age[p_index][j];
                                    block = j;
                                }
                            }
                            // Check if prefetched block exists in VC
                            for(int j=0; j<vc_size; j++) {
                                if(vc_valid[j] && p_index == vc_index[j] && p_tag == vc_tag[j]) {
                                    p_vc_hit = true;
                                    vc_block = j;
                                    break;
                                }
                            }
                            if(p_vc_hit) {
                                //swap blocks
                                uint64_t temp_tag = tag_store[p_index][block];
                                bool temp_dirty = dirty[p_index][block];
                                bool temp_p_bit = p_bit[p_index][block];
                                bool temp_valid = valid[p_index][block];

                                tag_store[p_index][block] = vc_tag[vc_block];
                                dirty[p_index][block] = vc_dirty[vc_block];
                                p_bit[p_index][block] = vc_p_bit[vc_block];
                                valid[p_index][block] = vc_valid[vc_block];

                                vc_tag[vc_block] = temp_tag;
                                vc_index[vc_block] = p_index;
                                vc_dirty[vc_block] = temp_dirty;
                                vc_p_bit[vc_block] = temp_p_bit;
                                vc_valid[vc_block] = temp_valid;

                                //No lru update due to prefetch
                                //No rw due to prefetch
                            }
                            else {
                                uint64_t max_age = 0;
                                int first_in;
                                for(int j = 0; j < vc_size; j++) {
                                    if(!vc_valid[j]) {
                                        first_in = j;
                                        break;
                                    }
                                    if(vc_age[j] >= max_age) {
                                        max_age = vc_age[j];    //max_age => first in
                                        first_in = j;
                                    }
                                }
                                // Evict Block from Cache and write back to memory if dirty
                                if(vc_dirty[first_in]) p_stats->write_backs++;

                                // Replace block in victim cache
                                vc_tag[first_in] = tag_store[p_index][block];
                                vc_index[first_in] = p_index;
                                vc_dirty[first_in] = dirty[p_index][block];
                                vc_valid[first_in] = valid[p_index][block];

                                // Update vc FIFO
                                for(int j = 0; j < vc_size; j++){
                                    vc_age[j]++;
                                }
                                vc_age[first_in] = 0;

                                // Put new block in cache
                                tag_store[p_index][block] = p_tag;
                                dirty[p_index][block] = false;
                                p_bit[p_index][block] = true;
                                valid[p_index][block] = true;
                                //No rw due to prefetch
                            }
                        }
                    }
                }
                else {
                    PendingStride = d;
                }
            } //endif (prefetch!=0)
        } //end if vc defined
        else {
        //vc undefined
            // Evict block from cache and write back if dirty
            if(dirty[index][lru_block]) p_stats->write_backs++;

            // Put new block in cache
            tag_store[index][lru_block] = tag;
            dirty[index][lru_block] = false;
            p_bit[index][lru_block] = false;
            valid[index][lru_block] = true;

            // Update block to MRU
            age[index][lru_block] = p_stats->accesses;

            // Update stats
            if(rw == WRITE) {
                dirty[index][lru_block] = true; //Write Allocate
                p_stats->writes++;
                p_stats->write_misses++;
                p_stats->write_misses_combined++;
            }
            else {
                p_stats->reads++;
                p_stats->read_misses++;
                p_stats->read_misses_combined++;
            }
            // Check Prefetcher
            if(p_degree != 0) {
                uint64_t now = (address & ~((1<<offset_size) -1));
                uint64_t p_offset,p_index,p_tag;
                int block;
                d = now - LastMissBlockAddress;
                LastMissBlockAddress = now;
                //Now prefetch
                if(d == PendingStride) {
                    for(int i=1; i<=p_degree; i++) {
                        p_stats->prefetched_blocks++;
                        uint64_t p_addr;
                        bool p_hit = false;
                        p_addr = now + (i*PendingStride);
                        p_offset = (p_addr & ((1 << offset_size) - 1));
                        p_index = (p_addr & (((1 << index_size) - 1) << offset_size)) >> offset_size;
                        p_tag = (p_addr & (((1 << tag_size) - 1)<< (offset_size + index_size))) >> (offset_size+index_size);
                        //Does prefetched block exist in cache?
                        for( int j=0; j<asso; j++) {
                            if( valid[p_index][j] && p_tag == tag_store[p_index][j] ) {
                                p_hit = true;
                                break;
                            }
                        }
                        if(!p_hit) {
                        //if prefetched block is not in cache
                            //find lru block to be replaced
                            uint64_t min_acc = p_stats->accesses;
                            for(int j=0; j<asso; j++) {
                                if(!valid[j]) {
                                    block = j;
                                    break;
                                }
                                if(age[p_index][j] <= min_acc) {
                                    min_acc = age[p_index][j];
                                    block = j;
                                }
                            }
                            if(dirty[p_index][block] == 1) {
                                p_stats->write_backs++;
                                //Write Back
                            }
                            // Put prefetched block in cache
                            tag_store[p_index][block] = p_tag;
                            dirty[p_index][block] = false;
                            p_bit[p_index][block] = true;
                            valid[p_index][block] = true;
                            //Update Age to lru-1
                            if(age[p_index][block] > 0) {
                                age[p_index][block]--;
                            }
                        }
                    }
                }
                else {
                    PendingStride = d;
                }
            } //endif (prefetch!=0)
        } //endelse vc undefined
    }//endelse Miss
}//end function cache_access


/*l*
 * Subroutine for cleaning up any outstanding memory operations and calculating overall statistics
 * such as miss rate or average access time.
 * XXX: You're responsible for completing this routine
 *
 * @p_stats Pointer to the statistics structure
 */
void complete_cache(cache_stats_t *p_stats) {
    if(vc_size == 0) { p_stats->vc_misses = p_stats->misses;  }
    p_stats->bytes_transferred = (p_stats->vc_misses + p_stats->write_backs + p_stats->prefetched_blocks)*(1<<offset_size);
    p_stats->miss_rate = (p_stats->read_misses + p_stats->write_misses)/(double)p_stats->accesses;
    p_stats->avg_access_time = p_stats->hit_time + (((p_stats->read_misses_combined + p_stats->write_misses_combined)/(double)p_stats->accesses) * p_stats->miss_penalty);

    //cout << temp1 <<endl;
}
