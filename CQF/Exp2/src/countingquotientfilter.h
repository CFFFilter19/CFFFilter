#ifndef COUNTINGQUOTIENTFILTER_H_INCLUDED
#define COUNTINGQUOTIENTFILTER_H_INCLUDED

#include "murmurhash3.h"

extern uint32_t FILTER_BIT_SIZE;
extern double   FALSE_POSITIVE_RATE;
extern uint32_t BITS_PER_FINGERPRINT;
extern uint32_t FILTER_MAX_ITEMS;
extern uint32_t FILTER_RAND_SEED;

#define FILTER_LOAD_FACTOR      0.75
#define SLOTS_PER_BLOCK         64
#define NULL_SLOT_INDEX         (uint32_t) ((1LL << 32) - 1)

class CountingQuotientFilter
{
public:
    uint32_t num_items;
    uint32_t num_hashes;
    uint32_t hash_seed;
    uint32_t num_blocks;
    uint32_t num_slots;
    uint32_t fingerprint_bits;
    uint8_t  * offset_array;
    uint64_t * occupied_array;
    uint64_t * runend_array;
    uint16_t * fingerprint_array;
    uint32_t filter_size;
    uint64_t insert_accesses;
    uint64_t lookup_accesses;
    uint64_t delete_accesses;

private:
    void GenerateHashSeed();
    void ComputeHashSlotIndexAndFingerprintKey(uint64_t item_key, uint32_t & slot_index, uint32_t & fingerprint_key);
    uint32_t ComputeOccupiedRankIndex(uint32_t slot_index);
    uint32_t ComputeOccupiedRankIndexInBlock(uint32_t slot_index);
    uint32_t ComputeOccupiedNextIndex(uint32_t slot_index);
    uint32_t ComputeRunendSelectIndex(uint32_t slot_index, uint32_t rank_index);
    uint32_t ComputeRunendSelectIndexInBlock(uint32_t slot_index, uint32_t rank_index);
    uint32_t ComputeRunendSelectIndexForLast(uint32_t slot_index, uint32_t rank_index);
    uint32_t ComputeFingerprintClusterIndex(uint32_t slot_index);
    void StoreFingerprintKey(uint32_t slot_index, uint32_t fingerprint_key);
    void FetchFingerprintKey(uint32_t slot_index, uint32_t & fingerprint_key);
    bool MoveUpOneRunendBit(uint32_t slot_index);
    bool MoveDownOneRunendBit(uint32_t slot_index);
    bool MoveUpOneFingerprint(uint32_t slot_index);
    bool MoveDownOneFingerprint(uint32_t slot_index);
    bool UpdateOffsetArray(uint32_t block_index);
    bool InsertFingerprintKey(uint32_t start_slot_index, uint32_t end_slot_index, uint32_t fingerprint_key);
    bool SearchFingerprintKey(uint32_t start_slot_index, uint32_t end_slot_index, uint32_t fingerprint_key);
    bool RemoveFingerprintKey(uint32_t start_slot_index, uint32_t end_slot_index, uint32_t fingerprint_key);

public:
    CountingQuotientFilter(uint32_t max_items);
    ~CountingQuotientFilter();
    bool Insert(uint64_t item_key);
    bool Lookup(uint64_t item_key);
    bool Delete(uint64_t item_key);
    uint32_t CalculateFilterSize();
    double CalculateBitsPerItem();
    void ResetMemoryAccesses();
    uint64_t CalculateInsertAccesses();
    uint64_t CalculateLookupAccesses();
    uint64_t CalculateDeleteAccesses();
    void WriteFilterLog();
};

#endif // COUNTINGQUOTIENTFILTER_H_INCLUDED
