#include "cachesim.hpp"
#include <iostream>
using namespace std;
//Global Variable Declaration
uint64_t **tag_store;   // Cache Tag Store
int **dirty;            // Block Dirty Bit
uint64_t **age;         // Block Age Tracker (minimum is LRU)
int **p_bit;            // Useful Prefetch Chcek Bit
int **valid;            // Block Valid Bit

int offset_size;        // b
int index_size;         // c-s-b
int tag_size;           // 64-c-s
int asso;               // 2^s

uint64_t *vc_tag;       // Victim Cache Tag Store
uint64_t *vc_index;     // Victim Cache Index Store
int *vc_dirty;          // Victim Cache Block Dirty Bit
int *vc_age;            // Victim Cache Block FIFO Tracker (maximum is FI)
int vc_size;            // Victim Cache Size (v)
int *vc_p_bit;          // Victim Cache Block Useful Prefetch Bit
int *vc_valid;          // Victim Cache Block Valid Bit

int prefetch;           // Prefetch Degree k
uint64_t LastMissBlockAddress = 0;
uint64_t PendingStride = 0;

/** Subroutine to initialize the stats **/
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
 * Subroutine for initializing the cache.
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
	prefetch = k;
	const int sets = asso;
	const int total_lines = 1 << (index_size);			// number of rows
    const int vc_lines = vc_size;

	//set rows
	tag_store = new uint64_t*[total_lines];
	dirty = new int*[total_lines];
	age = new uint64_t*[total_lines];
	p_bit = new int*[total_lines];
	valid = new int*[total_lines];

	//set victom cache rows and inititalize dirty bits to 0
	vc_tag = new uint64_t[vc_lines];
	vc_index = new uint64_t[vc_lines];
	vc_dirty = new int[vc_lines];
	vc_age = new int[vc_lines];
	vc_p_bit = new int[vc_lines];
	vc_valid = new int[vc_lines];
	for(int i=0; i<vc_lines; i++) {
        vc_age[i] = 0;
        vc_tag[i] = 0;
        vc_index[i] = 0;
        vc_dirty[i] = 0;
        vc_p_bit[i] = 0;
        vc_valid[i] = 0;
    }

	//set coloumns
	for (int i=0; i < total_lines; i++) {
		tag_store[i] = new uint64_t[sets];
		dirty[i] = new int[sets];
		age[i] = new uint64_t[sets];
		p_bit[i] = new int[sets];
		valid[i] = new int[sets];

		//initialize dirty bits to 0
		for(int j=0; j<sets; j++) {
            tag_store[i][j] = 0;
            dirty[i][j] = 0;
            age[i][j] = 0;
            p_bit[i][j] = 0;
            valid[i][j] = 0;
		}
	}
}
/**
 * Subroutine that simulates the cache one trace event at a time.
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
        //IT'S A HIT
        //update LRU bits
        age[index][hit_set] = p_stats->accesses;
        if(p_bit[index][hit_set] == 1) {
            p_bit[index][hit_set] = 0;
            p_stats->useful_prefetches++;
        }
        if(rw == WRITE) {
            dirty[index][hit_set] = 1;
            p_stats->writes++;
        }
        else { p_stats->reads++; }
    }
    else {
        //IT'S A MISS!
        p_stats->misses++;
        //Find block to be replaced
        uint64_t lru_age = p_stats->accesses + 10;
        int lru_block;    //min age => least recently used
        for(int i=0; i<asso; i++) {
            if(age[index][i] < lru_age) {
                lru_age = age[index][i];
                lru_block = i;
            }
        }
        if(vc_size != 0) {    //if vc exists
        //Victim cache Exists
            //Check VC for hit
            int vc_block;
            for(int i=0; i<vc_size; i++) {
                if( (index == vc_index[i]) && (tag == vc_tag[i]) && (vc_valid[i] == 1) ) {
                    vc_hit = true;
                    vc_block = i;
                    break;
                }
            }
            if(vc_hit) {
            //If VC hits
                if(vc_p_bit[vc_block]==1) {
                    vc_p_bit[vc_block] = 0;
                    p_stats->useful_prefetches++;
                }
                //swap blocks
                uint64_t temp_tag = tag_store[index][lru_block];
                int temp_dirty = dirty[index][lru_block];
                int temp_p_bit = p_bit[index][lru_block];
                int temp_valid = valid[index][lru_block];

                tag_store[index][lru_block] = vc_tag[vc_block];
                dirty[index][lru_block] = vc_dirty[vc_block];
                p_bit[index][lru_block] = vc_p_bit[vc_block];
                valid[index][lru_block] = vc_valid[vc_block];

                vc_tag[vc_block] = temp_tag;
                vc_dirty[vc_block] = temp_dirty;
                //vc_index[vc_block] = index;
                vc_p_bit[vc_block] = temp_p_bit;
                vc_valid[vc_block] = temp_valid;

                //update block to MRU
                age[index][lru_block] = p_stats->accesses;
                //update vc fifo
                for (int i=0; i<vc_size; i++){
                    vc_age[i]++;
                }
                vc_age[vc_block] = 0;
                //update dirty bit if write access
                if(rw == WRITE){
                    dirty[index][lru_block] = 1;
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
                    if(vc_age[i] >= max_age) {
                        max_age = vc_age[i];    //max_age => first in
                        first_in = i;
                    }
                }
                if(vc_dirty[first_in] == 1) {
                    p_stats->write_backs++;
                    vc_dirty[first_in] = 0;
                    //Write Back to Memory
                }
                //Replace block in victim cache
                vc_tag[first_in] = tag_store[index][lru_block];
                vc_index[first_in] = index;
                vc_dirty[first_in] = dirty[index][lru_block];
                vc_p_bit[first_in] = p_bit[index][lru_block];
                vc_valid[first_in] = valid[index][lru_block];
                //Update vc FIFO
                for(int i = 0; i<vc_size;i++){
                    vc_age[i]++;
                }
                vc_age[first_in] = 0;
                // Put new block in cache
                tag_store[index][lru_block] = tag;
                dirty[index][lru_block] = 0;
                p_bit[index][lru_block] = 0;
                valid[index][lru_block] = 1;
                //Update LRU bits
                age[index][lru_block] = p_stats->accesses;
                if(rw == WRITE) {
                    dirty[index][lru_block] = 1;    //Due to write allocate
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
            if(prefetch!=0) {
            //Is prefetching defined?
                uint64_t now = (address & ~((1<<offset_size) -1));
                int block;
                d = now - LastMissBlockAddress;
                LastMissBlockAddress = now;
                //Now prefetch
                if(d == PendingStride) {
                    for(int i=1; i<=prefetch; i++) {
                        p_stats->prefetched_blocks++;
                        uint64_t p_offset,p_index,p_tag;
                        uint64_t p_addr;
                        bool p_vc_hit = false;
                        bool p_hit = false;
                        p_addr = now + (i*PendingStride);
                        p_offset = (p_addr & ((1 << offset_size) - 1));
                        p_index = (p_addr & (((1 << index_size) - 1) << offset_size)) >> offset_size;
                        p_tag = (p_addr & (((1 << tag_size) - 1)<< (offset_size + index_size))) >> (offset_size+index_size);
                        uint64_t min_acc = p_stats->accesses+10;
                        //Does prefetched block exist in cache?
                        for( int j=0; j<asso; j++) {
                            if( (p_tag == tag_store[p_index][j]) && (valid[p_index][j] == 1) ) {
                                p_hit = true;
                                p_bit[p_index][j] = 1;
                                break;
                            }
                        }
                        if(!p_hit) {
                        //if prefetched block is not in cache
                            //find lru block to be replaced
                            min_acc = p_stats->accesses+10;
                            for(int j=0; j<asso; j++) {
                                if(age[p_index][j]<min_acc) {
                                        min_acc = age[p_index][j];
                                        block = j;
                                }
                            }
                            //Check if prefetched block exists in VC
                            for(int j=0; j<vc_size; j++) {
                                if(p_index == vc_index[j] && p_tag == vc_tag[j] && vc_valid[j] == 1) {
                                    p_vc_hit = true;
                                    vc_block = j;
                                    vc_p_bit[j] = 1;
                                    break;
                                }
                            }
                            if(p_vc_hit) {
                                //swap blocks
                                uint64_t temp_tag = tag_store[p_index][block];
                                int temp_dirty = dirty[p_index][block];
                                int temp_p_bit = dirty[p_index][block];
                                int temp_valid = valid[p_index][block];

                                tag_store[p_index][block] = vc_tag[vc_block];
                                dirty[p_index][block] = vc_dirty[vc_block];
                                //p_bit[p_index][block] = vc_p_bit[vc_block];
                                p_bit[p_index][block] = 1;
                                valid[p_index][block] = vc_valid[vc_block];

                                vc_tag[vc_block] = temp_tag;
                                //vc_index[vc_block] = p_index;
                                vc_dirty[vc_block] = temp_dirty;
                                vc_p_bit[vc_block] = temp_p_bit;
                                vc_valid[vc_block] = temp_valid;

                                //update vc fifo
                                for(int j=0; j<vc_size; j++) {
                                    vc_age[j]++;
                                }
                                vc_age[vc_block] = 0;
                                // No lru update due to prefetch
                                // Block is already LRU
                                age[p_index][block]--;
                                // No rw due to prefetch
                            }
                            else {
                                uint64_t max_age = 0;
                                int first_in;
                                for(int j = 0; j < vc_size; j++) {
                                    if(vc_age[j] >= max_age) {
                                        max_age = vc_age[j];    //max_age => first in
                                        first_in = j;
                                    }
                                }
                                if(vc_dirty[first_in] == 1) {
                                    p_stats->write_backs++;
                                    vc_dirty[first_in] = 0;
                                    //Write Back to Memory
                                }
                                //Replace block in victim cache
                                vc_tag[first_in] = tag_store[p_index][block];
                                vc_index[first_in] = p_index;
                                vc_dirty[first_in] = dirty[p_index][block];
                                vc_p_bit[first_in] = p_bit[p_index][block];
                                //Update vc FIFO
                                for(int j = 0; j < vc_size; j++){
                                    vc_age[j]++;
                                }
                                vc_age[first_in] = 0;
                                // Put new block in cache
                                tag_store[p_index][block] = p_tag;
                                dirty[p_index][block] = 0;
                                p_bit[p_index][block] = 1;
                                valid[p_index][block] = 1;
                                //Update Age to lru
                                if(age[p_index][block] != 0) {
                                    age[p_index][block]--;
                                }
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
            //Check lru block for dirty
            if(dirty[index][lru_block] == 1) {
                p_stats->write_backs++;
                dirty[index][lru_block] = 0;
                //WriteBack
            }
            //Put new block in cache
            tag_store[index][lru_block] = tag;
            dirty[index][lru_block] = 0;
            valid[index][lru_block] = 1;
            //Update LRU
            age[index][lru_block] = p_stats->accesses;

            if(rw == WRITE) {
                dirty[index][lru_block] = 1; //Write Allocate
                p_stats->writes++;
                p_stats->write_misses++;
                p_stats->write_misses_combined++;
            }
            else {
                p_stats->reads++;
                p_stats->read_misses++;
                p_stats->read_misses_combined++;
            }
            if(prefetch != 0) {
            //Is prefetching defined?
                uint64_t now = (address & ~((1<<offset_size) -1));
                int block;
                d = now - LastMissBlockAddress;
                LastMissBlockAddress = now;
                //Now prefetch
                if(d == PendingStride) {
                    for(int i=1; i<=prefetch; i++) {
                        p_stats->prefetched_blocks++;
                        uint64_t p_addr;
                        uint64_t p_offset,p_index,p_tag;
                        bool p_hit = false;
                        p_addr = now + (i*PendingStride);
                        p_offset = (p_addr & ((1 << offset_size) - 1));
                        p_index = (p_addr & (((1 << index_size) - 1) << offset_size)) >> offset_size;
                        p_tag = (p_addr & (((1 << tag_size) - 1)<< (offset_size + index_size))) >> (offset_size+index_size);
                        //Does prefetched block exist in cache?
                        for( int j=0; j<asso; j++) {
                            if( p_tag == tag_store[p_index][j] && valid[p_index][j] == 1) {
                                p_hit = true;
                                break;
                            }
                        }
                        if(!p_hit) {
                        //if prefetched block is not in cache
                            //find lru block to be replaced
                            uint64_t min_acc = p_stats->accesses;
                            for(int j=0; j<asso; j++) {
                                if(age[p_index][j] < min_acc) {
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
                            dirty[p_index][block] = 0;
                            p_bit[p_index][block] = 1;
                            valid[p_index][block] = 1;
                            //Update Age to lru-1
                            if(age[p_index][block] != 0) {
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


/*
 * Subroutine for cleaning up any outstanding memory operations and calculating overall statistics
 * such as miss rate or average access time.
 *
 * @p_stats Pointer to the statistics structure
 */
void complete_cache(cache_stats_t *p_stats) {

    if(vc_size == 0) { p_stats->vc_misses = p_stats->misses;  }
    p_stats->bytes_transferred = (p_stats->read_misses_combined + p_stats->write_misses_combined + p_stats->write_backs + p_stats->prefetched_blocks)*(1<<offset_size);
    p_stats->miss_rate = (p_stats->misses)/(double)p_stats->accesses;
    p_stats->avg_access_time = p_stats->hit_time + (p_stats->miss_rate * p_stats->miss_penalty);

}
