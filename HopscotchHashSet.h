#include <iostream>
#include <vector>
#include <cstdint>
#include <string>
#include <random>
#include <stdexcept>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <windows.h>
#include <excpt.h>
#include <ctime>
#include <unordered_set>
#include <chrono>
#include <numeric>
//#include <boost/functional/hash.hpp>

using std::vector;
using std::cin;
using std::cout;
using std::endl;
using std::pair;
using std::string;
using std::to_string;
using std::uniform_int_distribution;
using std::random_device;
using std::mt19937;
using std::swap;
using std::unordered_set;
using std::unordered_multiset;


// helper functions
inline uint32_t bit_set(uint32_t number, uint32_t n) {
    return number | ((uint32_t)1 << n);
}

inline uint32_t bit_clear(uint32_t number, uint32_t n) {
    return number & ~((uint32_t)1 << n);
}

inline uint32_t bit_toggle(uint32_t number, uint32_t n) {
    return number ^ ((uint32_t)1 << n);
}

inline bool bit_check(uint32_t number, uint32_t n) {
    return (number >> n) & (uint32_t)1;
}

inline uint32_t bit_set_to(uint32_t number, uint32_t n, bool x) {
    return (number & ~((uint32_t)1 << n)) | ((uint32_t)x << n);
}

// following functions change the number instead of returning it
inline void bit_set_change(uint32_t& number, uint32_t n) {
    number |= ((uint32_t)1 << n);
}

inline void bit_clear_change(uint32_t& number, uint32_t n) {
    number &= ~((uint32_t)1 << n);
}

inline void bit_toggle_change(uint32_t& number, uint32_t n) {
    number ^= ((uint32_t)1 << n);
}

inline void bit_set_to_change(uint32_t& number, uint32_t n, bool x) {
    number = (number & ~((uint32_t)1 << n)) | ((uint32_t)x << n);
}

uint32_t minbit(uint32_t x) { // returns minimal bit that equals 1
    return static_cast<uint32_t>(log2(x & ~(x-1)));
}

int myPow(int x, uint32_t p)
{
    if (p == 0) return 1;
    if (p == 1) return x;

    int tmp = myPow(x, p/2);
    if (p % 2 == 0) return tmp * tmp;
    else return x * tmp * tmp;
}

uint32_t myMod(int x, uint32_t p) {
    // Returns MATHEMATICALLY x % p (the number from 0 to p-1)
    return ((x % p) + p) % p;
}

const uint32_t Prime = 0x01000193; //   16777619
const uint32_t recommendedSeed  = 0x811C9DC5; // 2166136261

inline uint32_t fnv1a(unsigned char oneByte, uint32_t hash /*= Seed*/)
{
    return (oneByte ^ hash) * Prime;
}

uint32_t fnv1a(const void* data, size_t numBytes, uint32_t hash /*= Seed*/)
{
    const unsigned char* ptr = (const unsigned char*)data;
    while (numBytes--)
        hash = fnv1a(*ptr++, hash);
    return hash;
}

uint32_t myhash(int key, uint32_t seed) {
    return fnv1a(&key, sizeof(key), seed);
}

uint32_t generateSeed() {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<uint32_t> distrib(0, 4294967295); // from 0 to max uint32_t
    uint32_t newSeed = distrib(gen);
    return newSeed;
}

/*struct Key {
    int value;
};

bool operator==(const Key& a, const Key& b) {
    return (a.value == b.value);
}

uint32_t myhash(Key key) {
    return fnv1a(&key.value, sizeof(key.value));
}*/

class HopscotchHashSet {
private:
    static const uint32_t HOP_RANGE = 32; // Can't be increased without changing bitmap info
    static const uint32_t ADD_RANGE = 256;
    static const uint32_t MAX_SEGMENTS = 1024;
    static const uint32_t MAX_TRIES = 5;
    static const uint32_t Prime = 0x01000193; //   16777619
    static const uint32_t modPrime = 100000007; // hash will be mod this, then mod table size
    uint32_t Seed;

    // We'll need these 2 for every call of ADD
    uint32_t bad_bucket_ind; // index of bucket with -1's
    uint32_t bad_bucket_bitmap; // bitmap of that bucket

    // initially filled with key=-1 and bitmap=0
    vector<pair<int, uint32_t>> values; // key + bitmap that contains info about ith bucket


    void resize(); // double the size, rehash
    bool tryadd(int key); // add element without resize

public:
    void init(uint32_t size, uint32_t seed); // init table of this size and this seed
    bool contains(int key);
    bool add (int key); // add with resize if needed
    bool remove(int key); // TODO optional: add runtime errors if not exists
    void print(); // prints table
};

void HopscotchHashSet::print() {
    cout << "Table: ";
    for (auto v: values) {
        cout << v.first << " ";
    }
    cout << "\nBitmaps: ";
    for (auto v: values) {
        cout << v.second << " ";
    }
    cout << endl;
    cout << "Bad_bucket_ind: " << bad_bucket_ind << endl;
    cout << "Bad_bucket_bitmap: " << bad_bucket_bitmap << endl;
}

void HopscotchHashSet::init(uint32_t size = 1024, uint32_t seed = recommendedSeed) {
    vector<pair<int, uint32_t>> temp;
    if (size > myPow(2, 30)) // just above 1 billion
        throw std::runtime_error("Size is too big!");
    temp.resize(size, pair(-1, 0));
    values = temp;
    Seed = seed;
    bad_bucket_bitmap = 0;
    bad_bucket_ind = (myhash(-1, Seed) % modPrime) % values.size();
}

void HopscotchHashSet::resize() {
    //cout << "Resize happening" << endl;
    for (uint32_t iteration = 1; iteration <= MAX_TRIES; ++iteration) {
        //cout << "Iteration " << iteration << " of resize happening" << endl;
        uint32_t seed = generateSeed(); // generate new seed, since 2xing the size won't resolve hash collision -- just make them 2x better!

        // since operating big sized tables is just painful, increase the size LINEARLY
        // create new table with size 2^i * prev_size and previous seed, then add elements 1 by 1
        HopscotchHashSet newSet;
        //newSet.init(myPow(2, iteration) * values.size(), seed);
        newSet.init((2 + iteration) * values.size(), seed);
        //newSet.init(myPow(2, iteration) * values.size(), Seed);

        bool flag = true; // is resize successful
        for (uint32_t i = 0; i < values.size(); ++i) {
            if (values[i].first != -1 or (bad_bucket_ind <= i // add if elem != -1 or is a real -1
                                          and static_cast<int>(bad_bucket_ind) >= static_cast<int>(i) - static_cast<int>(HOP_RANGE)
                                          and bit_check(bad_bucket_bitmap, i - bad_bucket_ind))) {
                bool isAddSuccessful = newSet.tryadd(values[i].first);
                if (!isAddSuccessful) {
                    flag = false;
                    break;
                }
            }
        }
        if (flag) {
            values = newSet.values;
            Seed = newSet.Seed;
            bad_bucket_bitmap = newSet.bad_bucket_bitmap;
            bad_bucket_ind = newSet.bad_bucket_ind;
            return;
        }
    }
    // Failed to resize & rehash table for some reason
    string exceptionmsg = "Can not resize table";
    throw std::runtime_error("Error: " + exceptionmsg);
}

/*void HopscotchHashSet::shrink() {
    // create table with 1/2 size
    HopscotchHashSet newSet;
    newSet.init(values.size() / 2, Seed);

    bool flag = true; // is shrink successful
    for (auto & value : values) {
        bool isAddSuccessful = newSet.tryadd(value.first);
        if (!isAddSuccessful) {
            flag = false;
            break;
        }
    }
    if (flag) {
        values = newSet.values;
        Seed = newSet.Seed;
        return;
    }
    // otherwise we failed -- do nothing
}*/

bool HopscotchHashSet::contains(int key) {
    uint32_t bucket_ind = (myhash(key, Seed) % modPrime) % values.size();
    uint32_t bucket_bitmap = values[bucket_ind].second; // get bitmap
    for (uint32_t i = 0; i < HOP_RANGE; ++i) {
        // check if (bucket_ind + i)'th cell contains value corresponding to this bucket
        if (bit_check(bucket_bitmap, i) && values[bucket_ind + i].first == key) {
            return true;
        }
    }
    return false;
}

bool HopscotchHashSet::remove(int key) { // returns true if element was in set and false if wasn't
    uint32_t bucket_ind = (myhash(key, Seed) % modPrime) % values.size();
    uint32_t bucket_bitmap = values[bucket_ind].second; // get bitmap
    for (uint32_t i = 0; i < HOP_RANGE; ++i) {
        if (bit_check(bucket_bitmap, i) && values[bucket_ind + i].first == key) {
            values[bucket_ind + i].first = -1;
            bit_clear_change(values[bucket_ind].second, i);
            if (key == -1)
                bit_clear_change(bad_bucket_bitmap, i);
            return true;
        }
    }
    return false;
}

bool HopscotchHashSet::tryadd(int key) { // true if successful, false if failed
    //cout << "h";
    uint32_t bucket_ind = (myhash(key, Seed) % modPrime) % values.size();
    // uint32_t bucket_bitmap = values[bucket_ind].second; // get bitmap
    bool found = false;
    uint32_t freeaddind;

    for (uint32_t addind = 0; addind < std::min(static_cast<int>(ADD_RANGE),
                                                static_cast<int>(values.size()) - static_cast<int>(bucket_ind)); ++addind) {
        //cout << "i";
        if (values[bucket_ind + addind].first == -1) {
            // current index = bucket_ind + ind
            // equals bad_bucket_ind + ...
            //cout << "j";
            if (bad_bucket_ind > bucket_ind + addind or static_cast<int>(bad_bucket_ind) <
                                                        static_cast<int>(bucket_ind) + static_cast<int>(addind) - static_cast<int>(HOP_RANGE)) {
                //cout << "k";
                // bad_bucket is out of range for current cell, thus -1 is fake
                //cout << "Bad_bucket out of range. addind = " << addind  << ", bucket_ind = " << bucket_ind
                //<< ", bad_bucket_ind = " << bad_bucket_ind << endl;
                //cout << "1st bound == " << bucket_ind + addind << ", 2nd bound == "
                //<< static_cast<int>(bucket_ind) + static_cast<int>(addind) - static_cast<int>(HOP_RANGE) << endl;
                //cout << (bad_bucket_ind > bucket_ind + addind) << " " << (static_cast<int>(bad_bucket_ind) <
                //    static_cast<int>(bucket_ind) + static_cast<int>(addind) - static_cast<int>(HOP_RANGE)) << endl;
                found = true;
                freeaddind = addind;
                break;
            } else {
                //cout << "l";
                if (!bit_check(bad_bucket_bitmap, bucket_ind + addind - bad_bucket_ind)) {
                    //cout << "m";
                    //cout << "Checked bad_bucket_bitmap. Bit = " << bucket_ind + addind - bad_bucket_ind << endl;
                    // found 0 on bad_bucket_bitmap ==> -1 is fake
                    found = true;
                    freeaddind = addind;
                    break;
                }
                // else -1 is real, just do nothing
            }
        }
    }
    //cout << "n";
    if (!found)
        return false;
    // otherwise found free space at bucket_ind + found_ind
    //cout << "o";
    if (freeaddind < HOP_RANGE) {
        //cout << "p";
        // we can insert without moving
        values[bucket_ind + freeaddind].first = key;
        bit_set_change(values[bucket_ind].second, freeaddind);
        //cout << "q";
        if (key == -1)
            bad_bucket_bitmap = values[bad_bucket_ind].second;
        //cout << "r";
        return true;
    }

    // else we have to move elements
    //cout << "s";
    while (freeaddind >= HOP_RANGE) {
        //cout << "t";
        bool is_moved = false;
        // check from left to right to see if we can swap some element with free cell
        for (uint32_t i = HOP_RANGE - 1; i > 0; --i) {
            //cout << "u";
            // try to move free cell i places to the left
            assert(0 < bucket_ind + freeaddind - i);
            //assert(bucket_ind + freeaddind - i < values.size());
            uint32_t minind = minbit(values[bucket_ind + freeaddind - i].second); // index of minimal bit in bucket bucket_ind+freeind-i
            //cout << "v";
            if (values[bucket_ind + freeaddind - i].second && minind < i) {
                //cout << "w";
                // if there is a key in bucket freeind-i below index freeind
                // then we can move
                bit_clear_change(values[bucket_ind + freeaddind - i].second, minind);
                bit_set_change(values[bucket_ind + freeaddind - i].second, i);
                swap(values[bucket_ind + freeaddind - i + minind].first, values[bucket_ind + freeaddind].first);
                freeaddind = freeaddind - i + minind;
                is_moved = true;
                break;
            }
        }
        //cout << "x";
        if (!is_moved)
            // couldn't move free space
            return false;
    }
    //cout << "y";
    // now that we moved elements, free space is in range HOP_RANGE
    values[bucket_ind + freeaddind].first = key;
    bit_set_change(values[bucket_ind].second, freeaddind);
    //cout << "z";
    if (key == -1)
        bad_bucket_bitmap = values[bad_bucket_ind].second;
    //cout << "A";
    return true;
}

bool HopscotchHashSet::add(int key) { // true if no resize happened, false if resize
    //cout << "e";
    bool is_successful = tryadd(key);
    //cout << "g";
    if (is_successful) {
        //cout << "Added without resize" << endl;
        return true;
    }
    uint32_t iter_count = 0;
    //cout << "f";
    while (iter_count < MAX_TRIES) {
        resize();
        //cout << values.size() << endl;
        bool is_successful_now = tryadd(key);
        if (is_successful_now) {
            //cout << "Successful adding with resize" << endl;
            return false;
        }
        ++iter_count;
    }
    // Failed to add element
    uint32_t bucket_ind = (myhash(key, Seed) % modPrime) % values.size();
    if (values[bucket_ind].second + 1 == 0)
        throw std::runtime_error("33 elems with same hash");
    string exceptionmsg = "Can not add element";
    throw std::runtime_error("Error: " + exceptionmsg);
}


#endif //HOPSCOTCH_HOPSCOTCHHASHTABLE_H
