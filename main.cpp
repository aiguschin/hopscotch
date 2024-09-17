#include <iostream>
#include "tests.h"

using std::cout;
using std::endl;

int main() {
    cout << "Testing...\n\n";
    cout << "---------------1k entries:-------------------" << endl;
    testIntAdd(50000, 1'000);
    testIntRemove(500, 1'000, 1'000);
    testIntTrueContains(500, 1'000, 2'000);
    testIntFalseContains(500, 1'000, 2'000);
    cout << "---------------10k entries:-------------------" << endl;
    testIntAdd(5000, 10'000);
    testIntRemove(500, 10'000, 10'000);
    testIntTrueContains(50, 10'000, 20'000);
    testIntFalseContains(500, 10'000, 20'000);
    cout << "---------------100k entries:-------------------" << endl;
    testIntAdd(500, 100'000);
    testIntRemove(50, 100'000, 100'000);
    testIntTrueContains(50, 100'000, 200'000);
    testIntFalseContains(50, 100'000, 200'000);
    cout << "---------------1mil entries:-------------------" << endl;
    testIntAdd(50, 1'000'000);
    testIntRemove(50, 1'000'000, 1'000'000);
    testIntTrueContains(50, 1'000'000, 2'000'000);
    testIntFalseContains(50, 1'000'000, 2'000'000);
    cout << "---------------10mil entries:-------------------" << endl;
    testIntAdd(5, 10'000'000);
    testIntRemove(5, 10'000'000, 10'000'000);
    testIntTrueContains(5, 10'000'000, 20'000'000);
    testIntFalseContains(5, 10'000'000, 20'000'000);
    cout << "---------------100mil entries:-------------------" << endl;
    testIntAdd(1, 100'000'000);
    testIntRemove(1, 100'000'000, 100'000'000);
    testIntTrueContains(1, 100'000'000, 200'000'000);
    testIntFalseContains(1, 100'000'000, 200'000'000);
}
