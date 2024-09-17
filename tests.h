#ifndef TESTS_H
#define TESTS_H

#include "HopscotchHashSet.h"
#include <chrono>

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

    double hoptimesum = std::accumulate(hoptimes.begin(), hoptimes.end(), 0.0); // TODO change to reduce?
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
        unordered_set<int> stlset;
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

#endif /* TESTS_H */
