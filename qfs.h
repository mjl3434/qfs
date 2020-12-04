/*
Copyright (c) 2020 Marcus J. Larwill

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#ifndef QFS_H
#define QFS_H

#include <stdint.h>
#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include <fstream>
#include <list>

using std::vector;
using std::list;
using std::fstream;
using std::ifstream;
using std::ofstream;

typedef struct
{
    uint64_t number;
    uint8_t flags;
} prime_canidate_t;

class relations
{
public:
	inline relations(mpz_class f_xi_val, mpz_class xi_val) :
			f_xi(f_xi_val),
			xi(xi_val),
			remaining(f_xi_val) {}
    mpz_class f_xi;
    mpz_class xi;
    mpz_class remaining;
};

namespace qfs
{

class qfs
{

public:
    qfs(char* numberToFactor, char* bound, char* numberOfThreads);
    ~qfs();

    bool factor(void);
    void printResult(void);

private:

    void getPrimesUpToBound(void);
    void fillFactorBase(void);
    void sieveFactorBase(void);
    int legendreSymbol(mpz_class& n, uint64_t p);
    void sieveOfEratosthenes(void);

    fstream primesCacheFile;
    mpz_class numberToFactor;
    uint64_t bound;
    vector<uint64_t> primes;
    vector<relations> relationsVector;
    list<uint64_t> factorBase;
    unsigned int numberOfThreads;

    // Static Data
    static const char cacheFilePath[];
    static const unsigned int MAX_THREADS;
};

} // namespace qfs

#endif // QFS_H
