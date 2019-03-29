#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <fstream>
#include <math.h>
#include <bitset>
#include <x86intrin.h>

#include "countingquotientfilter.h"

using namespace std;

double FALSE_POSITIVE_RATE = 0.1;
uint32_t FILTER_RAND_SEED = 1;

#define BITS_PER_UINT       32
#define BITS_PER_CHAR       8
static const uint8_t bit_mask[BITS_PER_CHAR] =
{
    0x01,  //00000001
    0x02,  //00000010
    0x04,  //00000100
    0x08,  //00001000
    0x10,  //00010000
    0x20,  //00100000
    0x40,  //01000000
    0x80   //10000000
};

/**
 *
 * This is the implementation of class CountingQuotientFilter
 *
 */
CountingQuotientFilter::CountingQuotientFilter(uint32_t max_items)
{
    if(0 == max_items)
    {
        cout << "Error: Zero item in the set." << endl;
        exit(1);
    }

    num_items = 0;
    num_slots = (uint32_t) ceil(1.0 * max_items / FILTER_LOAD_FACTOR);
    num_blocks = (uint32_t) ceil(1.0 * num_slots / SLOTS_PER_BLOCK);
    num_slots = num_blocks * SLOTS_PER_BLOCK;
    fingerprint_bits = (uint32_t) ceil(log2(1.0 / FALSE_POSITIVE_RATE));

    offset_array = new uint8_t [num_blocks];
    memset(offset_array, 0, sizeof(uint8_t) * num_blocks);
    occupied_array = new uint64_t [num_blocks];
    memset(occupied_array, 0, sizeof(uint64_t) * num_blocks);
    runend_array = new uint64_t [num_blocks];
    memset(runend_array, 0, sizeof(uint64_t) * num_blocks);
    num_uints = (uint32_t) ceil(1.0 * num_slots * fingerprint_bits / BITS_PER_UINT);
    fingerprint_array = new uint32_t [num_uints];
    memset(fingerprint_array, 0, sizeof(uint32_t) * num_uints);
    filter_size = (num_blocks * (sizeof(uint64_t) + sizeof(uint64_t) + sizeof(uint8_t))) + (sizeof(uint32_t) * num_uints);

    num_hashes = 1;
    hash_seed = 0;
    GenerateHashSeed();

    ResetMemoryAccesses();
}

CountingQuotientFilter::~CountingQuotientFilter()
{
    delete [] offset_array;
    delete [] occupied_array;
    delete [] runend_array;
    delete [] fingerprint_array;
}

// generate a hash seed for the filter.
void CountingQuotientFilter::GenerateHashSeed()
{
    uint32_t salt_array[128] =
    {
        0xAAAAAAAA, 0x55555555, 0x33333333, 0xCCCCCCCC,
        0x66666666, 0x99999999, 0xB5B5B5B5, 0x4B4B4B4B,
        0xAA55AA55, 0x55335533, 0x33CC33CC, 0xCC66CC66,
        0x66996699, 0x99B599B5, 0xB54BB54B, 0x4BAA4BAA,
        0xAA33AA33, 0x55CC55CC, 0x33663366, 0xCC99CC99,
        0x66B566B5, 0x994B994B, 0xB5AAB5AA, 0xAAAAAA33,
        0x555555CC, 0x33333366, 0xCCCCCC99, 0x666666B5,
        0x9999994B, 0xB5B5B5AA, 0xFFFFFFFF, 0xFFFF0000,
        0xB823D5EB, 0xC1191CDF, 0xF623AEB3, 0xDB58499F,
        0xC8D42E70, 0xB173F616, 0xA91A5967, 0xDA427D63,
        0xB1E8A2EA, 0xF6C0D155, 0x4909FEA3, 0xA68CC6A7,
        0xC395E782, 0xA26057EB, 0x0CD5DA28, 0x467C5492,
        0xF15E6982, 0x61C6FAD3, 0x9615E352, 0x6E9E355A,
        0x689B563E, 0x0C9831A8, 0x6753C18B, 0xA622689B,
        0x8CA63C47, 0x42CC2884, 0x8E89919B, 0x6EDBD7D3,
        0x15B6796C, 0x1D6FDFE4, 0x63FF9092, 0xE7401432,
        0xEFFE9412, 0xAEAEDF79, 0x9F245A31, 0x83C136FC,
        0xC3DA4A8C, 0xA5112C8C, 0x5271F491, 0x9A948DAB,
        0xCEE59A8D, 0xB5F525AB, 0x59D13217, 0x24E7C331,
        0x697C2103, 0x84B0A460, 0x86156DA9, 0xAEF2AC68,
        0x23243DA5, 0x3F649643, 0x5FA495A8, 0x67710DF8,
        0x9A6C499E, 0xDCFB0227, 0x46A43433, 0x1832B07A,
        0xC46AFF3C, 0xB9C8FFF0, 0xC9500467, 0x34431BDF,
        0xB652432B, 0xE367F12B, 0x427F4C1B, 0x224C006E,
        0x2E7E5A89, 0x96F99AA5, 0x0BEB452A, 0x2FD87C39,
        0x74B2E1FB, 0x222EFD24, 0xF357F60C, 0x440FCB1E,
        0x8BBE030F, 0x6704DC29, 0x1144D12F, 0x948B1355,
        0x6D8FD7E9, 0x1C11A014, 0xADD1592F, 0xFB3C712E,
        0xFC77642F, 0xF9C4CE8C, 0x31312FB9, 0x08B0DD79,
        0x318FA6E7, 0xC040D23D, 0xC0589AA7, 0x0CA5C075,
        0xF874B172, 0x0CF914D5, 0x784D3280, 0x4E8CFEBC,
        0xC569F575, 0xCDB2A091, 0x2CC016B4, 0x5C5F4421
    };

    // randomly select a hash seed for hash functions.
    srand(FILTER_RAND_SEED);
    hash_seed = salt_array[(rand() % 128)];
}

// compute a hash slot index and a fingerprint key of an item.
void CountingQuotientFilter::ComputeHashSlotIndexAndFingerprintKey(uint64_t item_key, uint32_t & slot_index, uint32_t & fingerprint_key)
{
    uint64_t hash_value = MurmurHash2_64((const char *)&item_key, sizeof(item_key), hash_seed);
    slot_index = (uint32_t) ((hash_value >> 32) % num_slots);
    fingerprint_key = (uint32_t) (hash_value % (1LL << fingerprint_bits));
    if(0 == fingerprint_key)
    {
        fingerprint_key = 1;
    }
}

// compute the occupied rank index.
uint32_t CountingQuotientFilter::ComputeOccupiedRankIndex(uint32_t slot_index)
{
    uint32_t rank_index = 0;
    uint32_t block_index = slot_index / SLOTS_PER_BLOCK;
    uint32_t bit_index = slot_index % SLOTS_PER_BLOCK;

    // compute the pop count of each 64-bit block in the occupied array.
    for(size_t i = 0; i < block_index; i++)
    {
        rank_index += __builtin_popcountll(occupied_array[i]);
    }

    // compute the pop count of the current block.
    uint64_t bit_value = (uint64_t)((((1LL << bit_index) - 1) << 1) + 1);
    rank_index += __builtin_popcountll(occupied_array[block_index] & bit_value);

    return rank_index;
}

// compute the occupied rank index in block.
uint32_t CountingQuotientFilter::ComputeOccupiedRankIndexInBlock(uint32_t slot_index)
{
    uint32_t rank_index = 0;
    uint32_t block_index = slot_index / SLOTS_PER_BLOCK;
    uint32_t bit_index = slot_index % SLOTS_PER_BLOCK;

    // compute the pop count of the current block.
    uint64_t bit_value = (uint64_t)((((1LL << bit_index) - 1) << 1) + 1);
    rank_index += __builtin_popcountll(occupied_array[block_index] & bit_value);

    return rank_index;
}

// compute the occupied next index.
uint32_t CountingQuotientFilter::ComputeOccupiedNextIndex(uint32_t slot_index)
{
    uint32_t occupied_next_index = 0;
    uint32_t current_slot_index = slot_index + 1;

    while(current_slot_index < num_slots)
    {
        uint32_t block_index = current_slot_index / SLOTS_PER_BLOCK;
        uint32_t bit_index = current_slot_index % SLOTS_PER_BLOCK;
        uint64_t bit_value = (uint64_t) (1LL << bit_index);

        if(bit_value == (occupied_array[block_index] & bit_value))  // the occupied bit is 1.
        {
            occupied_next_index = current_slot_index;
            break;
        }
        else  // the occupied bit is 0.
        {
            current_slot_index++;
            occupied_next_index = current_slot_index;
        }
    }

    return occupied_next_index;
}

// compute the runend select index.
uint32_t CountingQuotientFilter::ComputeRunendSelectIndex(uint32_t slot_index, uint32_t rank_index)
{
    uint32_t select_index = 0;
    uint32_t block_index = slot_index / SLOTS_PER_BLOCK;

    if(0 == rank_index)
    {
        select_index = NULL_SLOT_INDEX;
        return select_index;
    }

    // compute the pop count of each block in the runend array.
    uint32_t cover_block_index = 0;
    uint32_t ones_count = 0;
    for(size_t i = 0; i < num_blocks; i++)
    {
        uint32_t pop_count = __builtin_popcountll(runend_array[i]);
        if((ones_count + pop_count) >= rank_index)
        {
            cover_block_index = i;
            break;
        }
        else  // ones_count + pop_count < rank_index
        {
            ones_count += pop_count;
        }
    }

    if(cover_block_index < block_index)
    {
        if(slot_index > 0)
        {
            select_index = slot_index - 1;
        }
        else
        {
            select_index = NULL_SLOT_INDEX;
        }

        return select_index;
    }
    else
    {
        select_index = cover_block_index * SLOTS_PER_BLOCK;
    }

    uint32_t remain_count = rank_index - ones_count;
    uint32_t low_value = (uint32_t)(runend_array[cover_block_index] % (1LL << 32));
    uint32_t pop_count = __builtin_popcountl(low_value);
    if(pop_count >= remain_count)
    {
        if(32 == remain_count)
        {
            select_index += 31;
        }
        else
        {
            uint32_t pdep_value = _pdep_u32((uint32_t)(1LL << (remain_count - 1)), low_value);
            select_index += __builtin_ctzl(pdep_value);
        }
    }
    else
    {
        select_index += 32;
        if(32 == (remain_count - pop_count))
        {
            select_index += 31;
        }
        else
        {
            uint32_t high_value = (uint32_t)(runend_array[cover_block_index] >> 32);
            uint32_t pdep_value = _pdep_u32((uint32_t)(1LL << (remain_count - pop_count - 1)), high_value);
            select_index += __builtin_ctzl(pdep_value);
        }
    }

    return select_index;
}

// compute the runend select index in block
uint32_t CountingQuotientFilter::ComputeRunendSelectIndexInBlock(uint32_t slot_index, uint32_t rank_index)
{
    uint32_t select_index = 0;
    uint32_t block_index = slot_index / SLOTS_PER_BLOCK;

    if(0 == rank_index)
    {
        select_index = NULL_SLOT_INDEX;
        return select_index;
    }

    if(rank_index == 1)
    {
        select_index = slot_index + (uint32_t)(offset_array[block_index]);
        return select_index;
    }

    //find the location of the first 1 in the occupied array.
    uint32_t start_slot_index = block_index * SLOTS_PER_BLOCK;
    uint32_t current_slot_index = start_slot_index;
    uint32_t end_slot_index = SLOTS_PER_BLOCK * (block_index + 1) - 1;
    while(current_slot_index <= end_slot_index)
    {
        uint64_t bit_v = (uint64_t) (1LL << (current_slot_index - start_slot_index));
        if(bit_v == (occupied_array[block_index] & bit_v))
        {
            break;
        }
        current_slot_index ++;
    }

    start_slot_index = current_slot_index + (uint32_t)(offset_array[block_index]);
    block_index = start_slot_index / SLOTS_PER_BLOCK;
    uint32_t bit_index = start_slot_index % SLOTS_PER_BLOCK;
    uint64_t bit_value = (uint64_t)(1LL << 63 );
    bit_value -= 1;
    uint64_t bit_v = (uint64_t)(1LL << 63);
    bit_value += bit_v;
    uint32_t i;
    for(i = 0; i < bit_index; i++)
    {
        bit_v = (uint64_t)(1LL << i);
        bit_value -= bit_v;
    }

    uint32_t current_rank = __builtin_popcountll(runend_array[block_index] & bit_value);
    uint32_t remain_count = rank_index;
    bool in_current_block = true;

    while(current_rank < remain_count)
    {
        block_index ++;
        remain_count -= current_rank;
        current_rank = __builtin_popcountll(runend_array[block_index]);
        in_current_block = false;
    }

    //the end of run is in this block
    if(in_current_block == true)
    {
        current_rank = 1;
        current_slot_index = start_slot_index + 1;
        end_slot_index = SLOTS_PER_BLOCK * (block_index + 1) - 1;
        while(current_slot_index <= end_slot_index)
        {
            uint64_t bit_v = (uint64_t) (1LL << (current_slot_index - block_index * SLOTS_PER_BLOCK));
            if(bit_v == (runend_array[block_index] & bit_v))
            {
                current_rank ++;
            }

            if(current_rank >= remain_count)
                break;
            current_slot_index ++;
        }
        select_index = current_slot_index;
    }
    else //the end of run is in another block
    {
        current_rank = 0;
        uint32_t current_slot_index =  SLOTS_PER_BLOCK * block_index ;
        uint32_t end_slot_index = SLOTS_PER_BLOCK * (block_index + 1) - 1;
        while(current_slot_index <= end_slot_index)
        {
            uint64_t bit_v = (uint64_t) (1LL << (current_slot_index - block_index * SLOTS_PER_BLOCK));
            if(bit_v == (runend_array[block_index] & bit_v))
                current_rank ++;
            if(current_rank >= remain_count)
                break;
            current_slot_index ++;
        }
        select_index = current_slot_index;
    }
    return select_index;
}

// compute the runend select index.
uint32_t CountingQuotientFilter::ComputeRunendSelectIndexForLast(uint32_t slot_index, uint32_t rank_index)
{
    uint32_t select_index = 0;
    uint32_t block_index = slot_index / SLOTS_PER_BLOCK;

    if(0 == rank_index)
    {
        select_index = NULL_SLOT_INDEX;
        return select_index;
    }

    // compute the pop count of each block in the runend array.
    uint32_t cover_block_index = 0;
    uint32_t ones_count = 0;
    for(size_t i = 0; i < num_blocks; i++)
    {
        uint32_t pop_count = __builtin_popcountll(runend_array[i]);
        if((ones_count + pop_count) >= rank_index)
        {
            cover_block_index = i;
            break;
        }
        else  // ones_count + pop_count < rank_index
        {
            ones_count += pop_count;
        }
    }
    if(cover_block_index < block_index)
    {
        if(slot_index > 0)
        {
            select_index = slot_index - 1;
        }
        else
        {
            select_index = NULL_SLOT_INDEX;
        }

        return select_index;
    }
    else
    {
        select_index = cover_block_index * SLOTS_PER_BLOCK;
    }

    uint32_t remain_count = rank_index - ones_count;
    uint32_t low_value = (uint32_t)(runend_array[cover_block_index] % (1LL << 32));
    uint32_t pop_count = __builtin_popcountl(low_value);
    if(pop_count >= remain_count)
    {
        if(32 == remain_count)
        {
            select_index += 31;
        }
        else
        {
            uint32_t pdep_value = _pdep_u32((uint32_t)(1LL << (remain_count - 1)), low_value);
            select_index += __builtin_ctzl(pdep_value);
        }
    }
    else
    {
        select_index += 32;
        if(32 == (remain_count - pop_count))
        {
            select_index += 31;
        }
        else
        {
            uint32_t high_value = (uint32_t)(runend_array[cover_block_index] >> 32);
            uint32_t pdep_value = _pdep_u32((uint32_t)(1LL << (remain_count - pop_count - 1)), high_value);
            select_index += __builtin_ctzl(pdep_value);
        }
    }

    return select_index;
}

// compute the fingerprint cluster index.
uint32_t CountingQuotientFilter::ComputeFingerprintClusterIndex(uint32_t slot_index)
{
    uint32_t cluster_slot_index = slot_index + 1;

    while(cluster_slot_index < num_slots)
    {
        // check each fingerprint key until it is 0.
        uint32_t fingerprint_value = 0;
        FetchFingerprintKey(cluster_slot_index, fingerprint_value);
        if(0 == fingerprint_value)
        {
            cluster_slot_index -= 1;
            break;
        }
        cluster_slot_index++;
    }

    return cluster_slot_index;
}

// store a fingerprint key into a slot.
void CountingQuotientFilter::StoreFingerprintKey(uint32_t slot_index, uint32_t fingerprint_key)
{
    uint32_t start_bit_index = slot_index * fingerprint_bits;
    uint32_t end_bit_index = start_bit_index + fingerprint_bits - 1;
    uint32_t start_uint_index = start_bit_index / BITS_PER_UINT;
    uint32_t start_intrin_index = start_bit_index % BITS_PER_UINT;
    uint32_t end_uint_index = end_bit_index / BITS_PER_UINT;
    uint32_t end_intrin_index = end_bit_index % BITS_PER_UINT;

    if(start_uint_index == end_uint_index)
    {
        uint32_t low_key = _bzhi_u32(fingerprint_array[start_uint_index], start_intrin_index);
        uint32_t high_key = (uint32_t)((((uint64_t)(fingerprint_array[start_uint_index])) >> (end_intrin_index + 1)) << (end_intrin_index + 1));
        fingerprint_array[start_uint_index] = low_key + (fingerprint_key << start_intrin_index) + high_key;
    }
    else
    {
        uint32_t low_key = _bzhi_u32(fingerprint_array[start_uint_index], start_intrin_index);
        uint32_t high_key = _bzhi_u32(fingerprint_key, (BITS_PER_UINT - start_bit_index));
        fingerprint_array[start_uint_index] = low_key + (high_key << start_intrin_index);

        low_key = (uint32_t)(fingerprint_key >> (BITS_PER_UINT - start_intrin_index));
        high_key = (uint32_t)((((uint64_t)(fingerprint_array[end_uint_index])) >> (end_intrin_index + 1)) << (end_intrin_index + 1));
        fingerprint_array[end_uint_index] = low_key + high_key;
    }
}

// fetch a fingerprint key from a slot.
void CountingQuotientFilter::FetchFingerprintKey(uint32_t slot_index, uint32_t & fingerprint_key)
{
    fingerprint_key = 0;

    uint32_t start_bit_index = slot_index * fingerprint_bits;
    uint32_t end_bit_index = start_bit_index + fingerprint_bits - 1;
    uint32_t start_uint_index = start_bit_index / BITS_PER_UINT;
    uint32_t start_intrin_index = start_bit_index % BITS_PER_UINT;
    uint32_t end_uint_index = end_bit_index / BITS_PER_UINT;
    uint32_t end_intrin_index = end_bit_index % BITS_PER_UINT;

    if(start_uint_index == end_uint_index)
    {
        fingerprint_key = _bextr_u32(fingerprint_array[start_uint_index], start_intrin_index, fingerprint_bits);
    }
    else
    {
        uint32_t low_key = _bextr_u32(fingerprint_array[start_uint_index], start_intrin_index, (BITS_PER_UINT - start_intrin_index));
        uint32_t high_key = _bextr_u32(fingerprint_array[end_uint_index], 0, (end_intrin_index + 1));
        fingerprint_key = (high_key << (BITS_PER_UINT - start_intrin_index)) + low_key;
    }
}

// move up one runend bit at the slot.
bool CountingQuotientFilter::MoveUpOneRunendBit(uint32_t slot_index)
{
    if(slot_index > (num_slots - 1))
    {
        return false;
    }

    uint32_t block_index = slot_index / SLOTS_PER_BLOCK;
    uint32_t bit_index = slot_index % SLOTS_PER_BLOCK;
    uint64_t bit_value = (uint64_t) (1LL << bit_index);

    if(bit_value == (runend_array[block_index] & bit_value))  // the current runend bit is 1.
    {
        if(slot_index == (num_slots - 1))  // the runend bit at the last slot is 1.
        {
            return false;
        }

        // reset the current runend bit to 0.
        runend_array[block_index] &= ~bit_value;

        // set the next runend bit to 1.
        block_index = (slot_index + 1) / SLOTS_PER_BLOCK;
        bit_index = (slot_index + 1) % SLOTS_PER_BLOCK;
        bit_value = (uint64_t) (1LL << bit_index);
        runend_array[block_index] |= bit_value;
    }
    else  // the current runend bit is 0.
    {
        if(slot_index == (num_slots - 1))  // the runend bit at the last slot is 0.
        {
            return true;
        }

        // reset the current runend bit to 0.
        runend_array[block_index] &= ~bit_value;

        // reset the next runend bit to 0.
        block_index = (slot_index + 1) / SLOTS_PER_BLOCK;
        bit_index = (slot_index + 1) % SLOTS_PER_BLOCK;
        bit_value = (uint64_t) (1LL << bit_index);
        runend_array[block_index] &= ~bit_value;
    }

    return true;
}

// move down one runend bit at the slot.
bool CountingQuotientFilter::MoveDownOneRunendBit(uint32_t slot_index)
{
    if(slot_index > (num_slots - 1))
    {
        return false;
    }

    if(0 == slot_index)  // reset the runend bit at the start slot.
    {
        // reset the runend bit at the start slot to 0.
        uint32_t block_index = slot_index / SLOTS_PER_BLOCK;
        uint32_t bit_index = slot_index % SLOTS_PER_BLOCK;
        uint64_t bit_value = (uint64_t) (1LL << bit_index);
        runend_array[block_index] &= ~bit_value;
        return true;
    }

    uint32_t block_index = slot_index / SLOTS_PER_BLOCK;
    uint32_t bit_index = slot_index % SLOTS_PER_BLOCK;
    uint64_t bit_value = (uint64_t) (1LL << bit_index);

    if(bit_value == (runend_array[block_index] & bit_value))  // the current runend bit is 1.
    {
        // reset the current runend bit to 0.
        runend_array[block_index] &= ~bit_value;

        // set the prev runend bit to 1.
        block_index = (slot_index - 1) / SLOTS_PER_BLOCK;
        bit_index = (slot_index - 1) % SLOTS_PER_BLOCK;
        bit_value = (uint64_t) (1LL << bit_index);
        runend_array[block_index] |= bit_value;
    }
    else  // the current runend bit is 0.
    {
        // reset the current runend bit to 0.
        runend_array[block_index] &= ~bit_value;

        // set the prev runend bit to 0.
        block_index = (slot_index - 1) / SLOTS_PER_BLOCK;
        bit_index = (slot_index - 1) % SLOTS_PER_BLOCK;
        bit_value = (uint64_t) (1LL << bit_index);
        runend_array[block_index] &= ~bit_value;
    }

    return true;
}

// move up one fingerprint key at the slot.
bool CountingQuotientFilter::MoveUpOneFingerprint(uint32_t slot_index)
{
    if(slot_index > (num_slots - 1))
    {
        return false;
    }

    if(slot_index == (num_slots - 1))  // move up the fingerprint key at the last slot.
    {
        uint32_t fingerprint_value = 0;
        FetchFingerprintKey(slot_index, fingerprint_value);
        if(0 != fingerprint_value)
        {
            return false;
        }
    }
    else  // slot_index < (num_slots - 1).
    {
        // store the current fingerprint key at the next slot.
        uint32_t fingerprint_value = 0;
        FetchFingerprintKey(slot_index, fingerprint_value);
        StoreFingerprintKey(slot_index + 1, fingerprint_value);
        // reset the current fingerprint key to 0.
        StoreFingerprintKey(slot_index, 0);
    }

    return true;
}

// move down one fingerprint key at the slot.
bool CountingQuotientFilter::MoveDownOneFingerprint(uint32_t slot_index)
{
    if(slot_index > (num_slots - 1))
    {
        return false;
    }

    if(0 == slot_index)  // reset the fingerprint key at the start slot.
    {
        // reset the fingerprint key at the start slot to 0.
        StoreFingerprintKey(slot_index, 0);
        return true;
    }

    // store the current fingerprint key at the prev slot.
    uint32_t fingerprint_value = 0;
    FetchFingerprintKey(slot_index, fingerprint_value);
    StoreFingerprintKey(slot_index - 1, fingerprint_value);
    // reset the current fingerprint key to 0.
    StoreFingerprintKey(slot_index, 0);

    return true;
}

//update offset array
bool CountingQuotientFilter::UpdateOffsetArray(uint32_t block_index)
{
    int i;
    uint64_t bit_value;
    bool update_result = false;
    for(i = 0; i < 64; i++)
    {
        bit_value = (uint64_t) (1LL << i);

        // find the first 1 in the block
        if(bit_value == (occupied_array[block_index] & bit_value))
        {
            uint32_t slot_index = block_index * SLOTS_PER_BLOCK + i;
            uint32_t rank_index = ComputeOccupiedRankIndex(slot_index);
            uint32_t end_slot_index = ComputeRunendSelectIndex(slot_index, rank_index);
            uint8_t new_offset = (uint8_t)(end_slot_index - slot_index);
            if(offset_array[block_index] != new_offset)
            {
                update_result = true;
                offset_array[block_index] = new_offset;
            }
            break;
        }
    }
    return update_result;
}

// insert a fingerprint key into a run of slots.
bool CountingQuotientFilter::InsertFingerprintKey(uint32_t start_slot_index, uint32_t end_slot_index, uint32_t fingerprint_key)
{
    bool insert_result = false;
    insert_accesses += 1;
    uint32_t empty_slot_index;

    uint32_t block_index = start_slot_index / SLOTS_PER_BLOCK;
    uint32_t bit_index = start_slot_index % SLOTS_PER_BLOCK;
    uint64_t bit_value = (uint64_t) (1LL << bit_index);

    if(NULL_SLOT_INDEX == end_slot_index)  // the slot index is NULL.
    {
        occupied_array[block_index] |= bit_value;  // set the occupied bit to 1.
        runend_array[block_index] |= bit_value;  // set the runend bit to 1.
        StoreFingerprintKey(start_slot_index, fingerprint_key);  // store the fingerprint key at the start slot.
        insert_result = true;
        UpdateOffsetArray(block_index);
        return insert_result;
    }

    if(bit_value != (occupied_array[block_index] & bit_value))  // the occupied bit is 0.
    {
        if(start_slot_index > end_slot_index)
        {
            occupied_array[block_index] |= bit_value;  // set the occupied bit to 1.
            runend_array[block_index] |= bit_value;  // set the runend bit to 1.
            StoreFingerprintKey(start_slot_index, fingerprint_key);  // store the fingerprint key at the start slot.
            insert_result = true;
            UpdateOffsetArray(block_index);
        }
        else  // start_slot_index <= end_slot_index.
        {
            empty_slot_index = ComputeFingerprintClusterIndex(end_slot_index) + 1;
            if((num_slots - 1) == end_slot_index ||(empty_slot_index >= num_slots))  // the filter is assumed to be full.
            {
                insert_result = false;
            }
            else  // move up runend bits and fingerprint keys.
            {
                // set the occupied bit at the start slot to 1.
                occupied_array[block_index] |= bit_value;

                // compute the fingerprint cluster index.
                uint32_t cluster_slot_index = ComputeFingerprintClusterIndex(end_slot_index);

                // move up runend bits and fingerprint keys from cluster slot index to end+1 slot index.
                insert_result = true;
                uint32_t current_slot_index = cluster_slot_index;
                while(current_slot_index > end_slot_index)
                {
                    if((true == MoveUpOneRunendBit(current_slot_index)) && (true == MoveUpOneFingerprint(current_slot_index)))
                    {
                        insert_result = true;
                        current_slot_index--;
                    }
                    else
                    {
                        insert_result = false;
                        break;
                    }
                }

                // set the runend bit and store the fingerprint key at the end+1 slot.
                if(true == insert_result)
                {
                    block_index = (end_slot_index + 1) / SLOTS_PER_BLOCK;
                    bit_index = (end_slot_index + 1) % SLOTS_PER_BLOCK;
                    bit_value = (uint64_t) (1LL << bit_index);
                    runend_array[block_index] |= bit_value;  // set the runend bit at the end+1 slot to 1.
                    StoreFingerprintKey(end_slot_index + 1, fingerprint_key);  // store the fingerprint key at the end+1 slot.

                    // update the offset array
                    uint32_t block_start_index = start_slot_index / SLOTS_PER_BLOCK;
                    UpdateOffsetArray(block_start_index);

                    uint32_t current_block_index = block_start_index + 1;
                    for(; current_block_index < num_blocks; current_block_index++ )
                    {
                        if(UpdateOffsetArray(current_block_index) == false)
                        {
                            break;
                        }
                    }
                }
            }
        }
    }
    else  // the occupied bit is 1.
    {
        empty_slot_index = ComputeFingerprintClusterIndex(end_slot_index) + 1;
        if((num_slots - 1) == end_slot_index ||(empty_slot_index >= num_slots))  // the filter is assumed to be full.
        {
            insert_result = false;
        }
        else  // we need to move up runend bits and fingerprint keys.
        {
            // compute the fingerprint cluster index.
            uint32_t cluster_slot_index = ComputeFingerprintClusterIndex(end_slot_index);

            // move up runend bits and fingerprint keys from cluster slot index to end+1 slot index.
            insert_result = true;
            uint32_t current_slot_index = cluster_slot_index;
            while(current_slot_index > end_slot_index)
            {
                if((true == MoveUpOneRunendBit(current_slot_index)) && (true == MoveUpOneFingerprint(current_slot_index)))
                {
                    insert_result = true;
                    current_slot_index--;
                }
                else
                {
                    insert_result = false;
                    break;
                }
            }

            // move up the runend bit at the end slot and store the fingerprint key at the end+1 slot.
            if(true == insert_result)
            {
                MoveUpOneRunendBit(end_slot_index);
                StoreFingerprintKey(end_slot_index + 1, fingerprint_key);

                // update the offset array
                uint32_t block_start_index = start_slot_index / SLOTS_PER_BLOCK;
                UpdateOffsetArray(block_start_index);

                uint32_t current_block_index = block_start_index + 1;
                for(; current_block_index < num_blocks; current_block_index++ )
                {
                    if(UpdateOffsetArray(current_block_index) == false)
                    {
                        break;
                    }
                }
            }
        }
    }
    return insert_result;
}

// search a fingerprint key in a run of slots.
bool CountingQuotientFilter::SearchFingerprintKey(uint32_t start_slot_index, uint32_t end_slot_index, uint32_t fingerprint_key)
{
    bool search_result = false;
    lookup_accesses += 1;

    if(NULL_SLOT_INDEX == end_slot_index)  // the slot index is NULL.
    {
        search_result = false;
        return search_result;
    }

    uint32_t block_index = start_slot_index / SLOTS_PER_BLOCK;
    uint32_t bit_index = start_slot_index % SLOTS_PER_BLOCK;
    uint64_t bit_value = (uint64_t) (1LL << bit_index);
    if(bit_value != (occupied_array[block_index] & bit_value))  // the occupied bit is 0.
    {
        search_result = false;
        return search_result;
    }

    // check each fingerprint key in the run of slots.
    uint32_t current_slot_index = end_slot_index;
    while(current_slot_index >= start_slot_index)
    {
        if(current_slot_index == end_slot_index)  // check the fingerprint key at the end slot.
        {
            uint32_t fingerprint_value = 0;
            FetchFingerprintKey(current_slot_index, fingerprint_value);
            if(fingerprint_key == fingerprint_value)
            {
                search_result = true;
                break;
            }
        }
        else  // check the fingerprint key at the middle slot in the run.
        {
            uint32_t block_index = current_slot_index / SLOTS_PER_BLOCK;
            uint32_t bit_index = current_slot_index % SLOTS_PER_BLOCK;
            uint64_t bit_value = (uint64_t) (1LL << bit_index);
            if(bit_value == (runend_array[block_index] & bit_value))  // the runend bit at the middle slot is 1.
            {
                search_result = false;
                break;
            }

            uint32_t fingerprint_value = 0;
            FetchFingerprintKey(current_slot_index, fingerprint_value);
            if(fingerprint_key == fingerprint_value)
            {
                search_result = true;
                break;
            }
        }

        if(current_slot_index > 0)
        {
            current_slot_index--;
        }
        else
        {
            break;
        }
    }

    return search_result;
}

// remove a fingerprint key from a run of slots.
bool CountingQuotientFilter::RemoveFingerprintKey(uint32_t start_slot_index, uint32_t end_slot_index, uint32_t fingerprint_key)
{
    bool remove_result = false;
    uint32_t remove_slot_index = 0;
    uint32_t fingerprint_count = 0;
    delete_accesses += 1;

    if(NULL_SLOT_INDEX == end_slot_index)  // the slot index is NULL.
    {
        remove_result = false;
        return remove_result;
    }

    uint32_t block_index = start_slot_index / SLOTS_PER_BLOCK;
    uint32_t bit_index = start_slot_index % SLOTS_PER_BLOCK;
    uint64_t bit_value = (uint64_t) (1LL << bit_index);
    if(bit_value != (occupied_array[block_index] & bit_value))  // the occupied bit is 0.
    {
        remove_result = false;
        return remove_result;
    }

    // check each fingerprint key in the run of slots to find the slot that matches the fingerprint key.
    uint32_t current_slot_index = end_slot_index;
    while(current_slot_index >= start_slot_index)
    {
        if(current_slot_index == end_slot_index)  // check the fingerprint key at the end slot.
        {
            uint32_t fingerprint_value = 0;
            FetchFingerprintKey(current_slot_index, fingerprint_value);
            if(fingerprint_key == fingerprint_value)
            {
                remove_slot_index = current_slot_index;
                remove_result = true;
            }
            fingerprint_count++;
        }
        else  // check the fingerprint key at the middle slot in the run.
        {
            uint32_t block_index = current_slot_index / SLOTS_PER_BLOCK;
            uint32_t bit_index = current_slot_index % SLOTS_PER_BLOCK;
            uint64_t bit_value = (uint64_t) (1LL << bit_index);
            if(bit_value == (runend_array[block_index] & bit_value))  // the runend bit at the middle slot is 1.
            {
                break;
            }
            else
            {
                fingerprint_count++;
            }

            uint32_t fingerprint_value = 0;
            FetchFingerprintKey(current_slot_index, fingerprint_value);
            if(fingerprint_key == fingerprint_value)
            {
                    remove_slot_index = current_slot_index;
                    remove_result = true;
                    break;
            }
        }

        if(current_slot_index > 0)
        {
            current_slot_index--;
        }
        else
        {
            break;
        }
    }

    // remove the fingerprint key and move down the runend bit at the slot.
    if(true == remove_result)
    {
        if(1 == fingerprint_count)  // there is one fingerprint key to be removed in the run.
        {
            // set the occupied bit at the start slot to 0.
            uint32_t block_index = start_slot_index / SLOTS_PER_BLOCK;
            uint32_t bit_index = start_slot_index % SLOTS_PER_BLOCK;
            uint64_t bit_value = (uint64_t) (1LL << bit_index);
            occupied_array[block_index] &= (~bit_value);

            // set the runend bit at the end slot to 0.
            block_index = end_slot_index / SLOTS_PER_BLOCK;
            bit_index = end_slot_index % SLOTS_PER_BLOCK;
            bit_value = (uint64_t) (1LL << bit_index);
            runend_array[block_index] &= (~bit_value);

            // reset the fingerprint key at the end slot to 0.
            StoreFingerprintKey(end_slot_index, 0);
        }
        else  // there are multiple fingerprint keys to be moved down.
        {
            if(remove_slot_index == end_slot_index)
            {
                // move down one runend bit at the end slot in the run.
                MoveDownOneRunendBit(end_slot_index);

                // reset the fingerprint key at the end slot to 0.
                StoreFingerprintKey(end_slot_index, 0);
            }
            else if (remove_slot_index < end_slot_index)
            {
                uint32_t current_slot_index = remove_slot_index + 1;
                while(current_slot_index <= end_slot_index)
                {
                    // move down one runend bit at the current slot in the run.
                    MoveDownOneRunendBit(current_slot_index);

                    // move down one fingerprint key at the current slot in the run.
                    MoveDownOneFingerprint(current_slot_index);

                    current_slot_index++;
                }
            }
        }

        uint32_t empty_slot_index = ComputeFingerprintClusterIndex(end_slot_index) + 1;
        uint32_t end_value = end_slot_index;
        if(empty_slot_index > num_slots)
        {
            empty_slot_index = num_slots;
        }

        while(end_value + 1 < empty_slot_index)
        {
            uint32_t occupied_next_index = ComputeOccupiedNextIndex(start_slot_index);
            uint32_t rank_index = ComputeOccupiedRankIndex(occupied_next_index);
            uint32_t end_slot_index = ComputeRunendSelectIndex(occupied_next_index, rank_index);
            uint32_t last_end_slot_index = ComputeRunendSelectIndexForLast(start_slot_index, rank_index - 1);
            start_slot_index = occupied_next_index;
            end_value = end_slot_index;

            if((last_end_slot_index != NULL_SLOT_INDEX)&&(last_end_slot_index + 2 > occupied_next_index))
            {
                for(current_slot_index = last_end_slot_index + 2; current_slot_index <= end_slot_index; current_slot_index++ )
                {
                    // move down one runend bit at the current slot in the next run.
                    MoveDownOneFingerprint(current_slot_index);

                    // move down one fingerprint key at the current slot in the next run.
                    MoveDownOneRunendBit(current_slot_index);
                }
            }
            else
            {
                     break;
            }
        }

        // update the offset array
        block_index = start_slot_index / SLOTS_PER_BLOCK;
        uint32_t block_end_index = empty_slot_index / SLOTS_PER_BLOCK;
        uint32_t current_block_index = block_index ;
        while(current_block_index <= block_end_index)
        {
            UpdateOffsetArray(current_block_index);
            current_block_index ++ ;
        }
    }

    return remove_result;
}

// insert an item into the filter.
bool CountingQuotientFilter::Insert(uint64_t item_key)
{
    bool insert_result = false;

    // compute the slot index and the fingerprint key of the item.
    uint32_t slot_index = 0;
    uint32_t fingerprint_key = 0;
    ComputeHashSlotIndexAndFingerprintKey(item_key, slot_index, fingerprint_key);

    // compute the start slot index and the end slot index of the run.
    uint32_t start_slot_index = slot_index;
    uint32_t rank_index = ComputeOccupiedRankIndex(slot_index);
    uint32_t end_slot_index = ComputeRunendSelectIndex(slot_index, rank_index);

    // insert the fingerprint key into the run of slots.
    insert_result = InsertFingerprintKey(start_slot_index, end_slot_index, fingerprint_key);
    if(true == insert_result)
    {
        num_items++;
    }

    return insert_result;
}

// lookup whether an item is in the filter.
bool CountingQuotientFilter::Lookup(uint64_t item_key)
{
    bool lookup_result = false;

    // compute the slot index and the fingerprint key of the item.
    uint32_t slot_index = 0;
    uint32_t fingerprint_key = 0;
    ComputeHashSlotIndexAndFingerprintKey(item_key, slot_index, fingerprint_key);

    uint32_t block_index = slot_index / SLOTS_PER_BLOCK;
    uint32_t bit_index = slot_index % SLOTS_PER_BLOCK;
    uint64_t bit_value = (uint64_t) (1LL << bit_index);

    // check the occupied id at the slot.
    if(bit_value != (occupied_array[block_index] & bit_value))  // the occupied bit is 0.
    {
        lookup_result = false;
        return lookup_result;
    }

    // compute the start slot index and the end slot index.
    uint32_t start_slot_index = slot_index;
    uint32_t rank_index = ComputeOccupiedRankIndexInBlock(slot_index);
    uint32_t end_slot_index = ComputeRunendSelectIndexInBlock(slot_index, rank_index);

    // search the fingerprint key at the run of slots.
    lookup_result = SearchFingerprintKey(start_slot_index, end_slot_index, fingerprint_key);

    return lookup_result;
}

// delete an item from the filter.
bool CountingQuotientFilter::Delete(uint64_t item_key)
{
    bool delete_result = false;

    // compute the slot index and the fingerprint key of the item.
    uint32_t slot_index = 0;
    uint32_t fingerprint_key = 0;
    ComputeHashSlotIndexAndFingerprintKey(item_key, slot_index, fingerprint_key);

    uint32_t block_index = slot_index / SLOTS_PER_BLOCK;
    uint32_t bit_index = slot_index % SLOTS_PER_BLOCK;
    uint64_t bit_value = (uint64_t) (1LL << bit_index);

    // check the occupied id at the slot.
    if(bit_value != (occupied_array[block_index] & bit_value))  // the occupied bit is 0.
    {
        delete_result = false;
        return delete_result;
    }

    // compute the start slot index and the end slot index.
    uint32_t start_slot_index = slot_index;
    uint32_t rank_index = ComputeOccupiedRankIndex(slot_index);
    uint32_t end_slot_index = ComputeRunendSelectIndex(slot_index, rank_index);

    // remove the fingerprint key from the run of slots.
    delete_result = RemoveFingerprintKey(start_slot_index, end_slot_index, fingerprint_key);
    if(true == delete_result)
    {
        num_items--;
    }

    return delete_result;
}

// calculate the theoretical size of the filter.
uint32_t CountingQuotientFilter::CalculateFilterSize()
{
    return filter_size;
}

// calculate the bits per item of the filter.
double CountingQuotientFilter::CalculateBitsPerItem()
{
    return (1.0 * filter_size * BITS_PER_CHAR / num_items);
}

// reset the operation memory accesses.
void CountingQuotientFilter::ResetMemoryAccesses()
{
    insert_accesses = 0;
    lookup_accesses = 0;
    delete_accesses = 0;
}

// calculate the insert memory accesses.
uint64_t CountingQuotientFilter::CalculateInsertAccesses()
{
    return insert_accesses;
}

// calculate the lookup memory accesses.
uint64_t CountingQuotientFilter::CalculateLookupAccesses()
{
    return lookup_accesses;
}

// calculate the delete memory accesses.
uint64_t CountingQuotientFilter::CalculateDeleteAccesses()
{
    return delete_accesses;
}

// write the log of the filter.
void CountingQuotientFilter::WriteFilterLog()
{
    ofstream bloom_result_txt;
    bloom_result_txt.open("results/bloom_result.txt", ios_base::app);

    bloom_result_txt << "false_positive_rate: " << FALSE_POSITIVE_RATE
                     << ", load_factor: " << FILTER_LOAD_FACTOR
                     << ", num_items: " << num_items
                     << ", num_hashes: " << num_hashes
                     << ", num_blocks: " << num_blocks
                     << ", num_slots: " << num_slots
                     << ", fingerprint_bits: " << fingerprint_bits
                     << ", filter_size (bytes): " << CalculateFilterSize()
                     << ", bits_per_item: " << CalculateBitsPerItem() << endl;

    bloom_result_txt.close();
}
