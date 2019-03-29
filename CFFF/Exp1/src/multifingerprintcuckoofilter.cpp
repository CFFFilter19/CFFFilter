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

#include "multifingerprintcuckoofilter.h"

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
 * This is the implementation of class MultiFingerprintCuckooFilter
 *
 */
MultiFingerprintCuckooFilter::MultiFingerprintCuckooFilter(uint32_t max_items)
{
    if(0 == max_items)
    {
        cout << "Error: Zero item in the set." << endl;
        exit(1);
    }

    num_items = 0;
    num_buckets = (uint32_t) ceil(1.0 * max_items / (CUCKOO_LOAD_FACTOR * SLOTS_PER_BUCKET));
    choice_bits = (uint32_t) ceil(log2(CUCKOO_NUM_HASHES));
    fingerprint_bits = (uint32_t) ceil(log2(1.0 * CUCKOO_NUM_HASHES * SLOTS_PER_BUCKET / FALSE_POSITIVE_RATE));

    num_uints = (uint32_t) ceil(1.0 * num_buckets * SLOTS_PER_BUCKET * fingerprint_bits / BITS_PER_UINT);
    slot_array = new uint32_t [num_uints];
    memset(slot_array, 0, sizeof(uint32_t) * num_uints);
    filter_size = sizeof(uint32_t) * num_uints;

    item_hash_seed = 0;
    fingerprint_hash_seed = 0;
    GenerateHashSeed();

    choice_hash_value_array = new uint32_t [CUCKOO_NUM_HASHES];
    memset(choice_hash_value_array, 0, sizeof(uint32_t) * CUCKOO_NUM_HASHES);
    ComputeChoiceHashValueArray();

    ResetMemoryAccesses();
}

MultiFingerprintCuckooFilter::~MultiFingerprintCuckooFilter()
{
    delete [] choice_hash_value_array;
    delete [] slot_array;
}

// generate an item hash seed and a fingerprint hash seed for the filter.
void MultiFingerprintCuckooFilter::GenerateHashSeed()
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

    // select a item seed and a fingerprint seed for hash functions.
    srand(FILTER_RAND_SEED);
    item_hash_seed = salt_array[(rand() % 128)];
    fingerprint_hash_seed = salt_array[(rand() % 128)];
}

// compute a choice hash value array.
void MultiFingerprintCuckooFilter::ComputeChoiceHashValueArray()
{
    // calculate each hash value of choices.
    for(uint32_t i = 0; i < CUCKOO_NUM_HASHES; i++)
    {
        uint32_t choice_value = i + 1;
        choice_hash_value_array[i] = MurmurHash3_32((const char *)&choice_value, sizeof(uint32_t), fingerprint_hash_seed);
    }
}

// compute a beacon bucket index and a fingerprint key of an item.
void MultiFingerprintCuckooFilter::ComputeBeaconBucketIndexFingerprintKey(uint64_t item_key, uint32_t & beacon_bucket_index, uint32_t & fingerprint_key)
{
    uint64_t hash_value = MurmurHash2_64((const char *)&item_key, sizeof(item_key), item_hash_seed);
    beacon_bucket_index = (uint32_t) ((hash_value >> 32) % num_buckets);
    fingerprint_key = (uint32_t) (hash_value % (1LL << fingerprint_bits));

    // calculate the remainder of the fingerprint key.
    uint32_t remainder_bits = fingerprint_bits - choice_bits;
    uint32_t remainder_value = fingerprint_key % (1LL << remainder_bits);
    if(0 == remainder_value)
    {
        remainder_value = 1;
    }

    // calculate the choice value of the fingerprint key.
    uint32_t choice_value = (fingerprint_key >> remainder_bits) % CUCKOO_NUM_HASHES;

    // compose the fingerprint key.
    fingerprint_key = remainder_value + (choice_value << remainder_bits);
}

// fetch the choice value and remainder value of a fingerprint key.
void MultiFingerprintCuckooFilter::FetchChoiceRemainderValue(uint32_t fingerprint_key, uint32_t & choice_value, uint32_t & remainder_value)
{
    uint32_t remainder_bits = fingerprint_bits - choice_bits;
    choice_value = (fingerprint_key >> remainder_bits) % CUCKOO_NUM_HASHES;
    remainder_value = fingerprint_key % (1LL << remainder_bits);
}

// compute the remainder hash value.
uint32_t MultiFingerprintCuckooFilter::ComputeRemainderHashValue(uint32_t remainder_value)
{
    // calculate the hash value of the remainder.
    return (MurmurHash3_32((const char *)&remainder_value, sizeof(uint32_t), fingerprint_hash_seed));
}

// compute the alternative bucket index by Addition.
uint32_t MultiFingerprintCuckooFilter::ComputeAlternativeBucketIndex(uint32_t beacon_bucket_index, uint32_t choice_value, uint32_t remainder_hash_value)
{
    // XOR the remainder hash value and choice hash value.
    uint64_t hash_value = remainder_hash_value ^ choice_hash_value_array[choice_value];

    // calculate the alternative bucket index of the fingerprint key.
    uint32_t alternative_bucket_index = ((uint64_t)beacon_bucket_index + hash_value) % num_buckets;

    return alternative_bucket_index;
}

// compute the beacon bucket index by Subtraction.
uint32_t MultiFingerprintCuckooFilter::ComputeBeaconBucketIndex(uint32_t bucket_index, uint32_t choice_value, uint32_t remainder_hash_value)
{
    // XOR the remainder hash value and choice hash value.
    uint64_t hash_value = remainder_hash_value ^ choice_hash_value_array[choice_value];

    // calculate the beacon bucket index of the fingerprint key.
    uint64_t beacon_bucket_index = (bucket_index + ((uint64_t)num_buckets << 32) - hash_value) % num_buckets;

    return beacon_bucket_index;
}

// compute the alternative fingerprint key.
uint32_t MultiFingerprintCuckooFilter::ComputeAlternativeFingerprintKey(uint32_t fingerprint_key)
{
    // calculate the remainder value of the fingerprint key.
    uint32_t remainder_bits = fingerprint_bits - choice_bits;
    uint32_t remainder_value = fingerprint_key % (1LL << remainder_bits);

    // calculate the choice value of the fingerprint key.
    uint32_t choice_value = (fingerprint_key >> remainder_bits) % CUCKOO_NUM_HASHES;
    choice_value = (choice_value + 1) % CUCKOO_NUM_HASHES;

    // compose the alternative fingerprint key.
    uint32_t alternative_fingerprint_key = remainder_value + (choice_value << remainder_bits);

    return alternative_fingerprint_key;
}

// store a fingerprint key into a bucket.
void MultiFingerprintCuckooFilter::StoreFingerprintKey(uint32_t bucket_index, uint32_t slot_index, uint32_t fingerprint_key)
{
    uint32_t start_bit_index = ((bucket_index * SLOTS_PER_BUCKET) + slot_index) * fingerprint_bits;
    uint32_t end_bit_index = start_bit_index + fingerprint_bits - 1;
    uint32_t start_uint_index = start_bit_index / BITS_PER_UINT;
    uint32_t start_intrin_index = start_bit_index % BITS_PER_UINT;
    uint32_t end_uint_index = end_bit_index / BITS_PER_UINT;
    uint32_t end_intrin_index = end_bit_index % BITS_PER_UINT;

    if(start_uint_index == end_uint_index)
    {
        uint32_t low_key = _bzhi_u32(slot_array[start_uint_index], start_intrin_index);
        uint32_t high_key = (uint32_t)((((uint64_t)(slot_array[start_uint_index])) >> (end_intrin_index + 1)) << (end_intrin_index + 1));
        slot_array[start_uint_index] = low_key + (fingerprint_key << start_intrin_index) + high_key;
    }
    else
    {
        uint32_t low_key = _bzhi_u32(slot_array[start_uint_index], start_intrin_index);
        uint32_t high_key = _bzhi_u32(fingerprint_key, (BITS_PER_UINT - start_bit_index));
        slot_array[start_uint_index] = low_key + (high_key << start_intrin_index);

        low_key = (uint32_t)(fingerprint_key >> (BITS_PER_UINT - start_intrin_index));
        high_key = (uint32_t)((((uint64_t)(slot_array[end_uint_index])) >> (end_intrin_index + 1)) << (end_intrin_index + 1));
        slot_array[end_uint_index] = low_key + high_key;
    }
}

// fetch a fingerprint key from a bucket.
void MultiFingerprintCuckooFilter::FetchFingerprintKey(uint32_t bucket_index, uint32_t slot_index, uint32_t & fingerprint_key)
{
    uint32_t start_bit_index = ((bucket_index * SLOTS_PER_BUCKET) + slot_index) * fingerprint_bits;
    uint32_t end_bit_index = start_bit_index + fingerprint_bits - 1;
    uint32_t start_uint_index = start_bit_index / BITS_PER_UINT;
    uint32_t start_intrin_index = start_bit_index % BITS_PER_UINT;
    uint32_t end_uint_index = end_bit_index / BITS_PER_UINT;
    uint32_t end_intrin_index = end_bit_index % BITS_PER_UINT;

    if(start_uint_index == end_uint_index)
    {
        fingerprint_key = _bextr_u32(slot_array[start_uint_index], start_intrin_index, fingerprint_bits);
    }
    else
    {
        uint32_t low_key = _bextr_u32(slot_array[start_uint_index], start_intrin_index, (BITS_PER_UINT - start_intrin_index));
        uint32_t high_key = _bextr_u32(slot_array[end_uint_index], 0, (end_intrin_index + 1));
        fingerprint_key = (high_key << (BITS_PER_UINT - start_intrin_index)) + low_key;
    }
}

// check whether a bucket has a vacant slot.
bool MultiFingerprintCuckooFilter::IsBucketVacant(uint32_t bucket_index, uint32_t & slot_index)
{
    bool vacant_result = false;
    insert_accesses += 1;

    // check each slot in the bucket.
    slot_index = 0;
    for(size_t i = 0; i < SLOTS_PER_BUCKET; i++)
    {
        uint32_t fingerprint_value = 0;
        FetchFingerprintKey(bucket_index, i, fingerprint_value);
        if(0 == fingerprint_value)
        {
            slot_index = i;
            vacant_result = true;
            break;
        }
    }

    return vacant_result;
}

// search a fingerprint key in the bucket.
bool MultiFingerprintCuckooFilter::SearchFingerprintKey(uint32_t bucket_index, uint32_t fingerprint_key)
{
    bool lookup_result = false;
    lookup_accesses += 1;

    // check each slot in the bucket.
    for(size_t i = 0; i < SLOTS_PER_BUCKET; i++)
    {
        uint32_t fingerprint_value = 0;
        FetchFingerprintKey(bucket_index, i, fingerprint_value);
        if((0 != fingerprint_value) && (fingerprint_key == fingerprint_value))
        {
            lookup_result = true;
            break;
        }
    }

    return lookup_result;
}

// remove a fingerprint key from a bucket.
bool MultiFingerprintCuckooFilter::RemoveFingerprintKey(uint32_t bucket_index, uint32_t fingerprint_key)
{
    bool remove_result = false;
    delete_accesses += 1;

    // check each slot in the bucket.
    for(size_t i = 0; i < SLOTS_PER_BUCKET; i++)
    {
        uint32_t fingerprint_value = 0;
        FetchFingerprintKey(bucket_index, i, fingerprint_value);
        if((0 != fingerprint_value) && (fingerprint_key == fingerprint_value))
        {
            StoreFingerprintKey(bucket_index, i, 0);
            remove_result = true;
            break;
        }
    }

    return remove_result;
}

// write the log of a bucket.
void MultiFingerprintCuckooFilter::WriteBucketLog(uint32_t bucket_index)
{
    ofstream bloom_result_txt;
    bloom_result_txt.open("results/bloom_result.txt", ios_base::app);
    bloom_result_txt << "Bucket [" << bucket_index << "]: " << endl;

    for(size_t i = 0; i < SLOTS_PER_BUCKET; i++)
    {
        uint32_t fingerprint_value = 0;
        FetchFingerprintKey(bucket_index, i, fingerprint_value);
        bitset <16> char_bitset(fingerprint_value);
        bloom_result_txt << "Slot [" << i << "]: " << fingerprint_value << "(" << char_bitset << ")"<< endl;
    }

    bloom_result_txt.close();
}


// move items along the cuckoo path for inserting an item.
bool MultiFingerprintCuckooFilter::CuckooMove(uint32_t bucket_index, uint32_t slot_index, uint32_t fingerprint_key, uint32_t & num_kicks)
{
    bool move_result = false;

    uint32_t * bucket_index_array = new uint32_t [CUCKOO_NUM_HASHES];
    uint32_t * fingerprint_key_array = new uint32_t [CUCKOO_NUM_HASHES];

    uint32_t insert_bucket_index = bucket_index;
    uint32_t insert_slot_index = slot_index;
    uint32_t insert_fingerprint_key = fingerprint_key;
    uint32_t insert_choice_value = 0;
    uint32_t insert_remainder_value = 0;
    FetchChoiceRemainderValue(insert_fingerprint_key, insert_choice_value, insert_remainder_value);
    uint32_t insert_remainder_hash_value = ComputeRemainderHashValue(insert_remainder_value);

    uint32_t victim_bucket_index = 0;
    uint32_t victim_slot_index = 0;
    uint32_t victim_fingerprint_key = 0;

    while((false == move_result) && (num_kicks <= CUCKOO_MAX_KICKS))
    {
        // compute the insert beacon bucket index of the fingerprint.
        uint32_t insert_beacon_bucket_index = ComputeBeaconBucketIndex(insert_bucket_index, insert_choice_value, insert_remainder_hash_value);

        bucket_index_array[0] = insert_bucket_index;
        fingerprint_key_array[0] = insert_fingerprint_key;
        uint32_t current_fingerprint_key = insert_fingerprint_key;

        // store the alternative fingerprint key into a vacant bucket.
        for(size_t i = 1; i < CUCKOO_NUM_HASHES; i++)
        {
            // compute the alternative fingerprint key, the choice value, and the alternative bucket index.
            uint32_t next_fingerprint_key = ComputeAlternativeFingerprintKey(current_fingerprint_key);
            insert_choice_value = (insert_choice_value + 1) % CUCKOO_NUM_HASHES;
            uint32_t next_bucket_index = ComputeAlternativeBucketIndex(insert_beacon_bucket_index, insert_choice_value, insert_remainder_hash_value);

            bucket_index_array[i] = next_bucket_index;
            fingerprint_key_array[i] = next_fingerprint_key;

            // check whether the next bucket has a vacant slot.
            uint32_t next_slot_index = 0;
            if(true == IsBucketVacant(next_bucket_index, next_slot_index))
            {
                // store the alternative fingerprint key into the next bucket.
                StoreFingerprintKey(next_bucket_index, next_slot_index, next_fingerprint_key);
                num_kicks++;
                move_result = true;
                break;
            }

            current_fingerprint_key = next_fingerprint_key;
        }

        // kick out other buckets to store the fingerprint key.
        if(false == move_result)
        {
            // calculate the victim bucket index of the item.
            uint32_t rand_index = (uint32_t)(rand() % CUCKOO_NUM_HASHES);
            victim_bucket_index = bucket_index_array[rand_index];

            // calculate the victim slot index of the item.
            if(insert_bucket_index == victim_bucket_index)
            {
                if(SLOTS_PER_BUCKET > 1)
                {
                    victim_slot_index = (uint32_t)(rand() % SLOTS_PER_BUCKET);
                    if(insert_slot_index == victim_slot_index)
                    {
                        victim_slot_index = (victim_slot_index + 1) % SLOTS_PER_BUCKET;
                    }
                }
                else
                {
                    rand_index = (rand_index + 1) % CUCKOO_NUM_HASHES;
                    victim_bucket_index = bucket_index_array[rand_index];
                    victim_slot_index = (uint32_t)(rand() % SLOTS_PER_BUCKET);
                }
            }
            else
            {
                victim_slot_index = (uint32_t)(rand() % SLOTS_PER_BUCKET);
            }

            // kick out a random fingerprint key from the victim bucket.
            FetchFingerprintKey(victim_bucket_index, victim_slot_index, victim_fingerprint_key);

            // store the fingerprint key into the victim bucket.
            insert_fingerprint_key = fingerprint_key_array[rand_index];
            StoreFingerprintKey(victim_bucket_index, victim_slot_index, insert_fingerprint_key);
            num_kicks++;

            insert_bucket_index = victim_bucket_index;
            insert_slot_index = victim_slot_index;
            insert_fingerprint_key = victim_fingerprint_key;
            FetchChoiceRemainderValue(insert_fingerprint_key, insert_choice_value, insert_remainder_value);
            insert_remainder_hash_value = ComputeRemainderHashValue(insert_remainder_value);
        }
    }

    delete [] bucket_index_array;
    delete [] fingerprint_key_array;
    return move_result;
}

// insert an item into the filter.
bool MultiFingerprintCuckooFilter::Insert(uint64_t item_key)
{
    bool insert_result = false;

    uint32_t * bucket_index_array = new uint32_t [CUCKOO_NUM_HASHES];
    uint32_t * fingerprint_key_array = new uint32_t [CUCKOO_NUM_HASHES];

    // compute the beacon bucket index and the fingerprint key of the item.
    uint32_t beacon_bucket_index = 0;
    uint32_t fingerprint_key = 0;
    ComputeBeaconBucketIndexFingerprintKey(item_key, beacon_bucket_index, fingerprint_key);

    // decompose the fingerprint key and fetch its choice value and remainder value.
    uint32_t choice_value = 0;
    uint32_t remainder_value = 0;
    FetchChoiceRemainderValue(fingerprint_key, choice_value, remainder_value);

    // compute the remainder hash value.
    uint32_t remainder_hash_value = ComputeRemainderHashValue(remainder_value);

    for(size_t i = 0; i < CUCKOO_NUM_HASHES; i++)
    {
        // compute each alternative bucket index of the item.
        uint32_t bucket_index = ComputeAlternativeBucketIndex(beacon_bucket_index, choice_value, remainder_hash_value);

        bucket_index_array[i] = bucket_index;
        fingerprint_key_array[i] = fingerprint_key;

        // check whether the alternative bucket has a vacant slot.
        uint32_t slot_index = 0;
        if(true == IsBucketVacant(bucket_index, slot_index))
        {
            // store the fingerprint key into the alternative bucket.
            StoreFingerprintKey(bucket_index, slot_index, fingerprint_key);
            num_items++;
            insert_result = true;
            break;
        }
        else
        {
            // compute the alternative fingerprint key of the item.
            fingerprint_key = ComputeAlternativeFingerprintKey(fingerprint_key);
            choice_value = (choice_value + 1) % CUCKOO_NUM_HASHES;
            insert_result = false;
        }
    }

    // kick out other items to insert the fingerprint key.
    if(false == insert_result)
    {
        uint32_t num_kicks = 0;
        uint32_t rand_index = (uint32_t)(rand() % CUCKOO_NUM_HASHES);
        uint32_t victim_bucket_index = bucket_index_array[rand_index];
        uint32_t victim_slot_index = (uint32_t)(rand() % SLOTS_PER_BUCKET);
        uint32_t victim_fingerprint_key = 0;
        FetchFingerprintKey(victim_bucket_index, victim_slot_index, victim_fingerprint_key);

        // store the fingerprint key into the victim bucket.
        uint32_t insert_fingerprint_key = fingerprint_key_array[rand_index];
        StoreFingerprintKey(victim_bucket_index, victim_slot_index, insert_fingerprint_key);
        num_kicks++;

        // kick out other fingerprint keys.
        insert_result = CuckooMove(victim_bucket_index, victim_slot_index, victim_fingerprint_key, num_kicks);
        if(true == insert_result)
        {
            num_items++;
        }
    }

    delete [] bucket_index_array;
    delete [] fingerprint_key_array;
    return insert_result;
}

// lookup whether an item is in the filter.
bool MultiFingerprintCuckooFilter::Lookup(uint64_t item_key)
{
    bool lookup_result = false;

    // compute the beacon bucket index and the fingerprint key of the item.
    uint32_t beacon_bucket_index = 0;
    uint32_t fingerprint_key = 0;
    ComputeBeaconBucketIndexFingerprintKey(item_key, beacon_bucket_index, fingerprint_key);

    // decompose the fingerprint key and fetch its choice value and remainder value.
    uint32_t choice_value = 0;
    uint32_t remainder_value = 0;
    FetchChoiceRemainderValue(fingerprint_key, choice_value, remainder_value);

    // compute the remainder hash value.
    uint32_t remainder_hash_value = ComputeRemainderHashValue(remainder_value);

    for(size_t i = 0; i < CUCKOO_NUM_HASHES; i++)
    {
        // compute each alternative bucket index of the item.
        uint32_t bucket_index = ComputeAlternativeBucketIndex(beacon_bucket_index, choice_value, remainder_hash_value);

        // search the fingerprint key in the alternative bucket.
        if(true == SearchFingerprintKey(bucket_index, fingerprint_key))
        {
            lookup_result = true;
            break;
        }
        else
        {
            // compute the alternative fingerprint key of the item.
            fingerprint_key = ComputeAlternativeFingerprintKey(fingerprint_key);
            choice_value = (choice_value + 1) % CUCKOO_NUM_HASHES;
            lookup_result = false;
        }
    }

    return lookup_result;
}

// delete an item from the filter.
bool MultiFingerprintCuckooFilter::Delete(uint64_t item_key)
{
    bool delete_result = false;

    // compute the beacon bucket index and the fingerprint key of the item.
    uint32_t beacon_bucket_index = 0;
    uint32_t fingerprint_key = 0;
    ComputeBeaconBucketIndexFingerprintKey(item_key, beacon_bucket_index, fingerprint_key);

    // decompose the fingerprint key and fetch its choice value and remainder value.
    uint32_t choice_value = 0;
    uint32_t remainder_value = 0;
    FetchChoiceRemainderValue(fingerprint_key, choice_value, remainder_value);

    // compute the remainder hash value.
    uint32_t remainder_hash_value = ComputeRemainderHashValue(remainder_value);

    for(size_t i = 0; i < CUCKOO_NUM_HASHES; i++)
    {
        // compute the alternative bucket index of the item.
        uint32_t bucket_index = ComputeAlternativeBucketIndex(beacon_bucket_index, choice_value, remainder_hash_value);

        // delete the fingerprint key from the alternative bucket.
        if(true == RemoveFingerprintKey(bucket_index, fingerprint_key))
        {
            num_items--;
            delete_result = true;
            break;
        }
        else
        {
            // calculate the alternative fingerprint key of the item.
            fingerprint_key = ComputeAlternativeFingerprintKey(fingerprint_key);
            choice_value = (choice_value + 1) % CUCKOO_NUM_HASHES;
            delete_result = false;
        }
    }

    return delete_result;
}

// write the log of a cuckoo move path.
void MultiFingerprintCuckooFilter::WriteCuckooMoveLog(uint32_t bucket_index, uint64_t slot_index, uint32_t fingerprint_key, uint32_t num_kicks)
{
    ofstream cuckoo_log_txt;
    cuckoo_log_txt.open("results/bloom_result.txt", ios_base::app);

    bitset <32> char_bitset(fingerprint_key);
    cuckoo_log_txt << "Move " << num_kicks << ": [" << bucket_index << ", " << slot_index << ", "
                    << fingerprint_key << "(" << char_bitset << ")]" << endl;
    cuckoo_log_txt.close();
}

// calculate the theoretical size of the filter.
uint32_t MultiFingerprintCuckooFilter::CalculateFilterSize()
{
    return filter_size;
}

// calculate the bits per item of the filter.
double MultiFingerprintCuckooFilter::CalculateBitsPerItem()
{
    return (1.0 * filter_size * BITS_PER_CHAR / num_items);
}

// reset the operation memory accesses.
void MultiFingerprintCuckooFilter::ResetMemoryAccesses()
{
    insert_accesses = 0;
    lookup_accesses = 0;
    delete_accesses = 0;
}

// calculate the insert memory accesses.
uint64_t MultiFingerprintCuckooFilter::CalculateInsertAccesses()
{
    return insert_accesses;
}

// calculate the lookup memory accesses.
uint64_t MultiFingerprintCuckooFilter::CalculateLookupAccesses()
{
    return lookup_accesses;
}

// calculate the delete memory accesses.
uint64_t MultiFingerprintCuckooFilter::CalculateDeleteAccesses()
{
    return delete_accesses;
}

// write the log of the filter.
void MultiFingerprintCuckooFilter::WriteFilterLog()
{
    ofstream bloom_result_txt;
    bloom_result_txt.open("results/bloom_result.txt", ios_base::app);

    bloom_result_txt << "false_positive_rate: " << FALSE_POSITIVE_RATE
                     << ", load_factor: " << CUCKOO_LOAD_FACTOR
                     << ", num_items: " << num_items
                     << ", num_hashes: " << CUCKOO_NUM_HASHES
                     << ", slots_per_bucket: " << SLOTS_PER_BUCKET
                     << ", num_buckets: " << num_buckets
                     << ", choice_bits: " << choice_bits
                     << ", fingerprint_bits: " << fingerprint_bits
                     << ", filter_size (bytes): " << CalculateFilterSize()
                     << ", bits_per_item: " << CalculateBitsPerItem() << endl;

    bloom_result_txt.close();
}
