#include <algorithm>
#include <cassert>
#include <chrono>
#include <climits>
#include <intrin.h>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#pragma intrinsic(_BitScanForward)

using std::cout;
using std::endl;
using std::pair;
using std::vector;
using std::random_device;
using std::string;
using std::swap;
using std::uniform_int_distribution;
using std::unordered_multiset;
using std::unordered_set;
using std::mt19937;

// helper functions TODO add intrinsics?
bool bit_check(uint32_t number, uint32_t n) {
    return (number >> n) & (uint32_t)1;
}

// following functions change the number instead of returning it
void bit_set_change(uint32_t& number, uint32_t n) {
    number |= ((uint32_t)1 << n);
}

void bit_clear_change(uint32_t& number, uint32_t n) {
    number &= ~((uint32_t)1 << n);
}

uint32_t minbit(uint32_t x) {
    unsigned long res;
    unsigned char isNonzero = _BitScanForward(&res, x);
    return res * isNonzero;
}

uint32_t math_mod(int x, int p) {
    // Returns MATHEMATICALLY x % p (the number from 0 to p-1)
    // p must be positive
    assert(p > 0);
    return ((x % p) + p) % p;
}

const uint32_t Prime = 0x01000193; //   16777619
const uint32_t default_seed = 0x811C9DC5; // 2166136261

uint32_t fnv1a(unsigned char oneByte, uint32_t hash /*= Seed*/)
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

/* WORKS:
 * HopscotchHashSet<int>
 * HopscotchHashSet<double>
 *
 * struct Point { int x, y, z; };
 * HopscotchHashSet<Point>
 *
 * and anything else that has correct sizeof
 *
 *
 * DOESN'T WORK:
 * HopscotchHashSet<std::string>
 * HopscotchHashSet<std::vector>
 *
 * and anything else that doesn't really support sizeof
 * */

template<typename T> uint32_t myhash(T key, uint32_t seed) {
    return fnv1a(&key, sizeof(key), seed);
}

uint32_t generate_seed() { // TODO move to xorshift
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<uint32_t> distrib(0, 4294967295); // from 0 to max uint32_t
    uint32_t newSeed = distrib(gen);
    return newSeed;
}

// https://en.wikipedia.org/wiki/Hopscotch_hashing
template<typename T> class HopscotchHashSet {
private:
    const T default_value{};

    // hopscotch parameters, initialized to defaults
    uint32_t HOP_RANGE = 32; // must be <= 32, default==32
    uint32_t ADD_RANGE = 256; // should be >= HOP_RANGE, default==256
    uint32_t MAX_TRIES = 5; // must be > 0, default==5
    uint32_t Seed = 0x811C9DC5;

    uint32_t bad_bucket_ind = 0; // index of bucket with ind==hash(default_value)
    uint32_t bad_bucket_bitmap = 0; // bitmap of that bucket
    uint32_t num_elements = 0; // number of elements

    // initially filled with key=default_key and bitmap=0
    vector<pair<T, uint32_t>> values; // key + bitmap that contains info about ith bucket

    bool is_resize_allowed = true; // if false, table will just die instead of resizing -- for testing purposes

    void resize(); // double the size, rehash
    bool tryadd(T key); // add element without resize

public:
    // create with default parameters
    HopscotchHashSet() {
        // use defaults
    }
    // create with specific parameters TODO validate params
    HopscotchHashSet(uint32_t init_HOP_RANGE, uint32_t init_ADD_RANGE, uint32_t init_MAX_TRIES, uint32_t init_Seed) {
        HOP_RANGE = init_HOP_RANGE;
        ADD_RANGE = init_ADD_RANGE;
        MAX_TRIES = init_MAX_TRIES;
        Seed = init_Seed;
    }
    void init(uint32_t size = 1024, uint32_t seed = default_seed); // init table of this size and this seed
    bool contains(T key) const;
    void add (T key); // add with resize if needed
    void remove(T key); // throws exception if no element found
    void print() const; // prints table
    void allow_resize(bool allow); // toggle is_resize_allowed

    [[nodiscard]] double load_factor() const; // get load factor of a table, 0 <= load_factor <= 1

    // for debugging purposes
    [[nodiscard]] vector<T> get_values() const; // returns values as vector
    [[nodiscard]] vector<uint32_t> get_bitmaps() const; // returns bitmaps
    [[nodiscard]] int get_num_elements() const; // return number of elements
};

template<typename T> int HopscotchHashSet<T>::get_num_elements() const {
    return num_elements;
}

template<typename T> vector<uint32_t> HopscotchHashSet<T>::get_bitmaps() const {
    vector<uint32_t> res;
    res.reserve(values.size());
    for (auto v : values) {
        res.push_back(v.second);
    }
    return res;
}

template<typename T> vector<T> HopscotchHashSet<T>::get_values() const {
    vector<T> res;
    res.reserve(values.size());
    for (auto v : values) {
        res.push_back(v.first);
    }
    return res;
}

template<typename T> double HopscotchHashSet<T>::load_factor() const {
    if (values.empty())
        return 0.0; // for empty table I think it makes sense
    return static_cast<double>(num_elements) / static_cast<double>(values.size());
}

template<typename T> void HopscotchHashSet<T>::allow_resize(bool allow) {
    is_resize_allowed = allow;
}

template<typename T> void HopscotchHashSet<T>::print() const {
    // T must be cout-able
    cout << "Table: ";
    for (const auto& v: values) {
        cout << v.first << " ";
    }
    cout << "\nBitmaps: ";
    for (const auto& v: values) {
        cout << v.second << " ";
    }
    cout << endl;
    cout << "Bad_bucket_ind: " << bad_bucket_ind << endl;
    cout << "Bad_bucket_bitmap: " << bad_bucket_bitmap << endl;
    cout << "Elements: " << num_elements << endl;
    cout << "Size: " << values.size() << endl;
}

template<typename T> void HopscotchHashSet<T>::init(uint32_t size, uint32_t seed) {
    vector<pair<T, uint32_t>> temp;
    temp.resize(size, pair(default_value, 0));
    swap(values, temp);
    Seed = seed;
    bad_bucket_bitmap = 0;
    bad_bucket_ind = size ? myhash(default_value, Seed) % size : 0; // init bad_bucket_ind with something if we init table of size 0 here
    num_elements = 0;
    is_resize_allowed = true;
}

template<typename T> void HopscotchHashSet<T>::resize() {
    if (!is_resize_allowed)
        throw std::runtime_error("Resize is not allowed!");
    for (uint32_t iteration = 1; iteration <= MAX_TRIES; ++iteration) {
        uint32_t seed = generate_seed(); // generate new seed,
        // because 2xing the size won't resolve hash collision -- just make 2x fewer collisions in any bucket

        // since operating big sized tables is just painful, increase the size linearly
        // create new table with size (2 + i) * prev_size and previous seed, then add elements 1 by 1
        HopscotchHashSet<T> newSet;
        newSet.init((2 + iteration) * values.size(), seed);

        bool flag = true; // is resize successful
        for (uint32_t i = 0; i < values.size(); ++i) {
            if (values[i].first != default_value || (math_mod(i - bad_bucket_ind, values.size()) < HOP_RANGE // add if elem != default_value or is a stored default_value
                    && bit_check(bad_bucket_bitmap, i - bad_bucket_ind))) {
                bool isAddSuccessful = newSet.tryadd(values[i].first);
                if (!isAddSuccessful) {
                    flag = false;
                    break;
                }
            }
        }
        if (flag) {
            swap(values, newSet.values);
            Seed = newSet.Seed;
            bad_bucket_bitmap = newSet.bad_bucket_bitmap;
            bad_bucket_ind = newSet.bad_bucket_ind;
            num_elements = newSet.num_elements;
            return;
        }
    }
    // Failed to resize & rehash table for some reason
    throw std::runtime_error("Error: Can not resize table");
}

template<typename T> bool HopscotchHashSet<T>::contains(T key) const {
    int size = static_cast<int>(values.size());
    uint32_t bucket_ind = myhash(key, Seed) % size;
    uint32_t bucket_bitmap = values[bucket_ind].second; // get bitmap
    // iterate through 1s in bucket_bitmap, check values inside
    while (bucket_bitmap) {
        uint32_t ind = minbit(bucket_bitmap);
        if (values[(bucket_ind + ind) % size].first == key) {
            return true;
        }
        bit_clear_change(bucket_bitmap, ind);
    }
    return false;
}

template<typename T> void HopscotchHashSet<T>::remove(T key) {
    int size = static_cast<int>(values.size());
    uint32_t bucket_ind = myhash(key, Seed) % values.size();
    uint32_t bucket_bitmap = values[bucket_ind].second; // get bitmap
    for (uint32_t i = 0; i < HOP_RANGE; ++i) { // minbit optimization slows things down here -- probably because overhead is too big
        if (bit_check(bucket_bitmap, i) && values[(bucket_ind + i) % values.size()].first == key) {
            values[(bucket_ind + i) % values.size()].first = default_value;
            bit_clear_change(values[bucket_ind].second, i);
            if (bucket_ind == bad_bucket_ind)
                bit_clear_change(bad_bucket_bitmap, i);
            //bit_set_to_change(bad_bucket_bitmap, i, bit_check(bad_bucket_bitmap, i) * (bucket_ind == bad_bucket_ind));
            // the line above works slower than if -- tested
            --num_elements;
            return;
        }
    }

    // key not found
    throw std::runtime_error("Tried to remove non-existent element");
}

template<typename T> bool HopscotchHashSet<T>::tryadd(T key) { // true if successful, false if failed
    if (values.empty()) {
        pair<T, uint32_t> temp = pair(key, 1);
        values.push_back(temp);
        bad_bucket_bitmap = 1;
        num_elements = 1;
        return true;
    }
    int size = static_cast<int>(values.size());
    uint32_t bucket_ind = myhash(key, Seed) % size;
    bool found = false;
    uint32_t freeaddind;

    for (uint32_t addind = 0; addind < static_cast<int>(ADD_RANGE); ++addind) {
        // current index = bucket_ind + addind
        // equals bad_bucket_ind + ..., 0 <= ... < HOP_RANGE
        if (values[(bucket_ind + addind) % size].first == default_value && (math_mod(static_cast<int>(bucket_ind + addind - bad_bucket_ind), size) >= HOP_RANGE
                || !bit_check(bad_bucket_bitmap,math_mod(static_cast<int>(bucket_ind + addind - bad_bucket_ind), size)))) {
            // if default_value is not stored, then we found free space
            found = true;
            freeaddind = addind;
            break;
        }
        // else default_value is stored, just do nothing
    }

    if (!found)
        return false;

    // otherwise found free space at bucket_ind + found_ind
    if (freeaddind < HOP_RANGE) {
        // we can insert without moving
        values[(bucket_ind + freeaddind) % size].first = key;
        bit_set_change(values[bucket_ind].second, freeaddind);
        if (bucket_ind == bad_bucket_ind)
            bad_bucket_bitmap = values[bad_bucket_ind].second;
        ++num_elements;
        return true;
    }

    // else we have to move elements=
    while (freeaddind >= HOP_RANGE) {
        bool is_moved = false;
        // check from left to right to see if we can swap some element with free cell
        for (uint32_t i = HOP_RANGE - 1; i > 0; --i) {
            // look at bucket that is i places to the left, check if we can move something from this bucket to the right
            uint32_t check_ind = math_mod(static_cast<int>(bucket_ind + freeaddind - i), size);
            uint32_t minind = minbit(values[check_ind].second); // index of minimal bit in bucket check_ind
            if (values[check_ind].second && minind < i) {
                // if there is a key in bucket bucket_ind + freeaddind - i below index bucket_ind + freeaddind
                // then we can move
                bit_clear_change(values[check_ind].second, minind);
                bit_set_change(values[check_ind].second, i);
                swap(values[(check_ind + minind) % size].first,
                     values[(check_ind + i) % size].first);
                freeaddind = freeaddind - i + minind;
                is_moved = true;
                break;
            }
        }
        if (!is_moved)
            // couldn't move free space
            return false;
    }
    // now that we moved elements, free space is in range HOP_RANGE
    values[(bucket_ind + freeaddind) % size].first = key;
    bit_set_change(values[bucket_ind].second, freeaddind);
    if (bucket_ind == bad_bucket_ind)
        bad_bucket_bitmap = values[bad_bucket_ind].second;
    ++num_elements;
    return true;
}

template<typename T> void HopscotchHashSet<T>::add(T key) { // true if no resize happened, false if resize
    bool is_successful = tryadd(key);
    if (is_successful)
        return;
    uint32_t iter_count = 0;
    while (iter_count < MAX_TRIES) {
        resize();
        bool is_successful_now = tryadd(key);
        if (is_successful_now)
            return;
        ++iter_count;
    }
    // Failed to add element
    uint32_t bucket_ind = myhash(key, Seed) % values.size();
    T elem = values[bucket_ind].first;
    bool flag = true; // check if we have HOP_RANGE + 1 equal elems
    for (int i = 1; i < HOP_RANGE; ++i) {
        uint32_t check_ind = math_mod(i + bucket_ind, values.size());
        if (values[check_ind].first != elem) {
            flag = false;
            break;
        }
    }
    if (flag)
        throw std::runtime_error("Error: You can't store more than HOP_RANGE equal elements");
    throw std::runtime_error("Error: Can not add element");
}

void testIntAdd(int numTries, int numElems) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // times
    vector<double> hoptimes;
    vector<double> stltimes;

    // load factors
    vector<double> hoploads;
    vector<double> stlloads;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distrib(INT_MIN, INT_MAX);

    for (int i = 0; i < numTries; ++i) {
        // generate elems
        vector<int> keys;
        keys.reserve(numElems);
        for (int j = 0; j < numElems; ++j) {
            keys.push_back(distrib(gen));
        }

        // do testing
        std::chrono::steady_clock::time_point hopstart = std::chrono::steady_clock::now();
        HopscotchHashSet<int> hoptable;
        for (auto key : keys) {
            hoptable.add(key);
        }
        std::chrono::steady_clock::time_point hopend = std::chrono::steady_clock::now();
        unordered_multiset<int> stlset;
        for (auto key : keys) {
            stlset.insert(key);
        }
        std::chrono::steady_clock::time_point stlend = std::chrono::steady_clock::now();

        std::chrono::duration<double> hoptime = hopend - hopstart;
        std::chrono::duration<double> stltime = stlend - hopend;

        hoptimes.push_back(hoptime.count());
        stltimes.push_back(stltime.count());

        hoploads.push_back(hoptable.load_factor());
        stlloads.push_back(stlset.load_factor());
    }

    double hoptimesum = std::accumulate(hoptimes.begin(), hoptimes.end(), 0.0);
    double hoptimemean = hoptimesum / static_cast<double>(hoptimes.size());
    double hoptimesq_sum = std::inner_product(hoptimes.begin(), hoptimes.end(), hoptimes.begin(), 0.0);
    double hoptimestdev = std::sqrt(hoptimesq_sum / static_cast<double>(hoptimes.size()) - hoptimemean * hoptimemean);

    double stltimesum = std::accumulate(stltimes.begin(), stltimes.end(), 0.0);
    double stltimemean = stltimesum / static_cast<double>(stltimes.size());
    double stltimesq_sum = std::inner_product(stltimes.begin(), stltimes.end(), stltimes.begin(), 0.0);
    double stltimestdev = std::sqrt(stltimesq_sum / static_cast<double>(stltimes.size()) - stltimemean * stltimemean);

    double hoploadsum = std::accumulate(hoploads.begin(), hoploads.end(), 0.0);
    double hoploadmean = hoploadsum / static_cast<double>(hoploads.size());
    double hoploadsq_sum = std::inner_product(hoploads.begin(), hoploads.end(), hoploads.begin(), 0.0);
    double hoploadstdev = std::sqrt(hoploadsq_sum / static_cast<double>(hoploads.size()) - hoploadmean * hoploadmean);

    double stlloadsum = std::accumulate(stlloads.begin(), stlloads.end(), 0.0);
    double stlloadmean = stlloadsum / static_cast<double>(stlloads.size());
    double stlloadsq_sum = std::inner_product(stlloads.begin(), stlloads.end(), stlloads.begin(), 0.0);
    double stlloadstdev = std::sqrt(stlloadsq_sum / static_cast<double>(stlloads.size()) - stlloadmean * stlloadmean);

    cout << "Testing adds:" << endl;
    cout << "Average Hopscotch time: " << hoptimemean << endl;
    cout << "Hopscotch time STD: " << hoptimestdev << endl;
    cout << "Average STL time: " << stltimemean << endl;
    cout << "STL time STD: " << stltimestdev << endl;
    cout << "Average Hopscotch load: " << hoploadmean << endl;
    cout << "Hopscotch load STD: " << hoploadstdev << endl;
    cout << "Average STL load: " << stlloadmean << endl;
    cout << "STL load STD: " << stlloadstdev << endl;
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
    cout << "Time elapsed: " << time_elapsed.count() << "s" << endl;
    cout << "Ratio: " << stltimemean / hoptimemean << endl;
    cout << endl;
}

void testIntRemove(int numTries, int numElems, int numChecks) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // sanity check
    assert(numElems >= numChecks);

    // times
    vector<double> hoptimes;
    vector<double> stltimes;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distrib(INT_MIN, INT_MAX);

    for (int i = 0; i < numTries; ++i) {
        // generate elems
        vector<int> keys;
        keys.reserve(numElems);
        for (int j = 0; j < numElems; ++j) {
            keys.push_back(distrib(gen));
        }

        // get test
        vector<int> testkeys;
        testkeys = keys;
        auto rng = std::default_random_engine {rd()};
        std::shuffle(std::begin(testkeys), std::end(testkeys), rng);
        HopscotchHashSet<int> hoptable;
        for (auto key : keys) {
            hoptable.add(key);
        }

        // do testing
        std::chrono::steady_clock::time_point hopstart = std::chrono::steady_clock::now();
        for (int j = 0; j < numChecks; ++j) {
            hoptable.remove(testkeys[j]);
        }
        std::chrono::steady_clock::time_point hopend = std::chrono::steady_clock::now();
        unordered_multiset<int> stlset;
        for (auto key : keys) {
            stlset.insert(key);
        }
        std::chrono::steady_clock::time_point stlstart = std::chrono::steady_clock::now();
        for (int j = 0; j < numChecks; ++j) {
            stlset.erase(stlset.find(testkeys[j]));
        }
        std::chrono::steady_clock::time_point stlend = std::chrono::steady_clock::now();

        std::chrono::duration<double> hoptime = hopend - hopstart;
        std::chrono::duration<double> stltime = stlend - stlstart;

        hoptimes.push_back(hoptime.count());
        stltimes.push_back(stltime.count());
    }

    double hopsum = std::accumulate(hoptimes.begin(), hoptimes.end(), 0.0);
    double hopmean = hopsum / static_cast<double>(hoptimes.size());
    double hopsq_sum = std::inner_product(hoptimes.begin(), hoptimes.end(), hoptimes.begin(), 0.0);
    double hopstdev = std::sqrt(hopsq_sum / static_cast<double>(hoptimes.size()) - hopmean * hopmean);

    double stlsum = std::accumulate(stltimes.begin(), stltimes.end(), 0.0);
    double stlmean = stlsum / static_cast<double>(stltimes.size());
    double stlsq_sum = std::inner_product(stltimes.begin(), stltimes.end(), stltimes.begin(), 0.0);
    double stlstdev = std::sqrt(stlsq_sum / static_cast<double>(stltimes.size()) - stlmean * stlmean);

    cout << "Testing removes:" << endl;
    cout << "Average Hopscotch time: " << hopmean << endl;
    cout << "Hopscotch STD: " << hopstdev << endl;
    cout << "Average STL time: " << stlmean << endl;
    cout << "STL STD: " << stlstdev << endl;
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
    cout << "Time elapsed: " << time_elapsed.count() << "s" << endl;
    cout << "Ratio: " << stlmean / hopmean << endl;
    cout << endl;
}

void testIntTrueContains(int numTries, int numElems, int numChecks) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // times
    vector<double> hoptimes;
    vector<double> stltimes;

    // loads
    vector<double> hoploads;
    vector<double> stlloads;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distrib(INT_MIN, INT_MAX);

    bool temp;
    int counter = 0; // to force compiler to do something on release

    for (int i = 0; i < numTries; ++i) {
        // generate elems
        vector<int> keys;
        keys.reserve(numElems);
        for (int j = 0; j < numElems; ++j) {
            keys.push_back(distrib(gen));
        }

        // get test
        vector<int> testkeys;
        testkeys = keys;
        auto rng = std::default_random_engine {rd()};
        std::shuffle(std::begin(testkeys), std::end(testkeys), rng);

        // do testing
        HopscotchHashSet<int> hoptable;
        unordered_multiset<int> stlset;
        for (auto key : keys) {
            hoptable.add(key);
        }
        for (auto key : keys) {
            stlset.insert(key);
        }

        std::chrono::steady_clock::time_point hopstart = std::chrono::steady_clock::now();
        for (int j = 0; j < numChecks; ++j) {
            temp = hoptable.contains(testkeys[j % numElems]);
            if (temp)
                ++counter;
        }
        std::chrono::steady_clock::time_point hopend = std::chrono::steady_clock::now();
        std::chrono::steady_clock::time_point stlstart = std::chrono::steady_clock::now();
        for (int j = 0; j < numChecks; ++j) {
            temp = stlset.contains(testkeys[j % numElems]);
            if (temp)
                ++counter;
        }
        std::chrono::steady_clock::time_point stlend = std::chrono::steady_clock::now();

        std::chrono::duration<double> hoptime = hopend - hopstart;
        std::chrono::duration<double> stltime = stlend - stlstart;

        hoptimes.push_back(hoptime.count());
        stltimes.push_back(stltime.count());

        hoploads.push_back(hoptable.load_factor());
        stlloads.push_back(stlset.load_factor());
    }

    double hoptimesum = std::accumulate(hoptimes.begin(), hoptimes.end(), 0.0);
    double hoptimemean = hoptimesum / static_cast<double>(hoptimes.size());
    double hoptimesq_sum = std::inner_product(hoptimes.begin(), hoptimes.end(), hoptimes.begin(), 0.0);
    double hoptimestdev = std::sqrt(hoptimesq_sum / static_cast<double>(hoptimes.size()) - hoptimemean * hoptimemean);

    double stltimesum = std::accumulate(stltimes.begin(), stltimes.end(), 0.0);
    double stltimemean = stltimesum / static_cast<double>(stltimes.size());
    double stltimesq_sum = std::inner_product(stltimes.begin(), stltimes.end(), stltimes.begin(), 0.0);
    double stltimestdev = std::sqrt(stltimesq_sum / static_cast<double>(stltimes.size()) - stltimemean * stltimemean);

    double hoploadsum = std::accumulate(hoploads.begin(), hoploads.end(), 0.0);
    double hoploadmean = hoploadsum / static_cast<double>(hoploads.size());
    double hoploadsq_sum = std::inner_product(hoploads.begin(), hoploads.end(), hoploads.begin(), 0.0);
    double hoploadstdev = std::sqrt(hoploadsq_sum / static_cast<double>(hoploads.size()) - hoploadmean * hoploadmean);

    double stlloadsum = std::accumulate(stlloads.begin(), stlloads.end(), 0.0);
    double stlloadmean = stlloadsum / static_cast<double>(stlloads.size());
    double stlloadsq_sum = std::inner_product(stlloads.begin(), stlloads.end(), stlloads.begin(), 0.0);
    double stlloadstdev = std::sqrt(stlloadsq_sum / static_cast<double>(stlloads.size()) - stlloadmean * stlloadmean);

    cout << "Testing true contains:" << endl;
    cout << "Average Hopscotch time: " << hoptimemean << endl;
    cout << "Hopscotch time STD: " << hoptimestdev << endl;
    cout << "Average STL time: " << stltimemean << endl;
    cout << "STL time STD: " << stltimestdev << endl;
    /*cout << "Average Hopscotch load: " << hoploadmean << endl;
    cout << "Hopscotch load STD: " << hoploadstdev << endl;
    cout << "Average STL load: " << stlloadmean << endl;
    cout << "STL load STD: " << stlloadstdev << endl;*/
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
    cout << "Time elapsed: " << time_elapsed.count() << "s" << endl;
    cout << "Counter: " << counter << endl;
    cout << "Ratio: " << stltimemean / hoptimemean << endl;
    cout << endl;
}

void testIntFalseContains(int numTries, int numElems, int numChecks) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // times
    vector<double> hoptimes;
    vector<double> stltimes;

    // loads
    vector<double> hoploads;
    vector<double> stlloads;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distrib(INT_MIN, INT_MAX);

    bool temp;
    int counter = 0; // to force compiler to do something on release

    for (int i = 0; i < numTries; ++i) {
        // generate elems
        vector<int> keys;
        keys.reserve(numElems);
        for (int j = 0; j < numElems; ++j) {
            keys.push_back(distrib(gen));
        }

        // fill stlset first! required for testing
        unordered_multiset<int> stlset;
        for (auto key : keys) {
            stlset.insert(key);
        }

        // generate test
        vector<int> testkeys;
        testkeys.reserve(numChecks);
        while (static_cast<int>(testkeys.size()) < numChecks) {
            int testcase = distrib(gen);
            if (!stlset.contains(testcase))
                testkeys.push_back(testcase);
        }

        // do testing
        HopscotchHashSet<int> hoptable;
        for (auto key : keys) {
            hoptable.add(key);
        }
        std::chrono::steady_clock::time_point hopstart = std::chrono::steady_clock::now();
        for (auto key: testkeys) {
            temp = hoptable.contains(key);
            if (!temp)
                ++counter;
        }
        std::chrono::steady_clock::time_point hopend = std::chrono::steady_clock::now();

        std::chrono::steady_clock::time_point stlstart = std::chrono::steady_clock::now();
        for (auto& key : testkeys) {
            temp = stlset.contains(key);
            if (!temp)
                ++counter;
        }
        std::chrono::steady_clock::time_point stlend = std::chrono::steady_clock::now();

        std::chrono::duration<double> hoptime = hopend - hopstart;
        std::chrono::duration<double> stltime = stlend - stlstart;

        hoptimes.push_back(hoptime.count());
        stltimes.push_back(stltime.count());

        hoploads.push_back(hoptable.load_factor());
        stlloads.push_back(stlset.load_factor());
    }

    double hoptimesum = std::accumulate(hoptimes.begin(), hoptimes.end(), 0.0);
    double hoptimemean = hoptimesum / static_cast<double>(hoptimes.size());
    double hoptimesq_sum = std::inner_product(hoptimes.begin(), hoptimes.end(), hoptimes.begin(), 0.0);
    double hoptimestdev = std::sqrt(hoptimesq_sum / static_cast<double>(hoptimes.size()) - hoptimemean * hoptimemean);

    double stltimesum = std::accumulate(stltimes.begin(), stltimes.end(), 0.0);
    double stltimemean = stltimesum / static_cast<double>(stltimes.size());
    double stltimesq_sum = std::inner_product(stltimes.begin(), stltimes.end(), stltimes.begin(), 0.0);
    double stltimestdev = std::sqrt(stltimesq_sum / static_cast<double>(stltimes.size()) - stltimemean * stltimemean);

    double hoploadsum = std::accumulate(hoploads.begin(), hoploads.end(), 0.0);
    double hoploadmean = hoploadsum / static_cast<double>(hoploads.size());
    double hoploadsq_sum = std::inner_product(hoploads.begin(), hoploads.end(), hoploads.begin(), 0.0);
    double hoploadstdev = std::sqrt(hoploadsq_sum / static_cast<double>(hoploads.size()) - hoploadmean * hoploadmean);

    double stlloadsum = std::accumulate(stlloads.begin(), stlloads.end(), 0.0);
    double stlloadmean = stlloadsum / static_cast<double>(stlloads.size());
    double stlloadsq_sum = std::inner_product(stlloads.begin(), stlloads.end(), stlloads.begin(), 0.0);
    double stlloadstdev = std::sqrt(stlloadsq_sum / static_cast<double>(stlloads.size()) - stlloadmean * stlloadmean);

    cout << "Testing false contains:" << endl;
    cout << "Average Hopscotch time: " << hoptimemean << endl;
    cout << "Hopscotch time STD: " << hoptimestdev << endl;
    cout << "Average STL time: " << stltimemean << endl;
    cout << "STL time STD: " << stltimestdev << endl;
    /*cout << "Average Hopscotch load: " << hoploadmean << endl;
    cout << "Hopscotch load STD: " << hoploadstdev << endl;
    cout << "Average STL load: " << stlloadmean << endl;
    cout << "STL load STD: " << stlloadstdev << endl;*/
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
    cout << "Time elapsed: " << time_elapsed.count() << "s" << endl;
    cout << "Counter: " << counter << endl;
    cout << "Ratio: " << stltimemean / hoptimemean << endl;
    cout << endl;
}

int main() {
    // TODO split in several files (HopscotchHashSet into .h, testing functions into another one)
    cout << "Testing...\n\n";
    testIntAdd(100, 1'000'000);
    testIntRemove(100, 1'000'000, 1'000'000);
    testIntTrueContains(100, 1'000'000, 1'000'001);
    testIntFalseContains(100, 1'000'000, 1'000'001);
}
