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
#include "qfs.h"

#include <stdexcept>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iostream> // FIXME: Remove
#include <unistd.h>

using std::string;
using std::ofstream;
using std::cout;

namespace qfs
{


enum PrimeCanidateFlags
{
    PRIME             = 0x01,
    COMPOSITE         = 0x02,
    QUADRATIC_RESIDUE = 0x04
};

const char qfs::cacheFilePath[] = "listOfPrimes.txt";
const unsigned int qfs::MAX_THREADS = 8;

qfs::qfs(char* theNumberToFactor, char* theBound, char* theNumberOfThreads) :
		numberToFactor(theNumberToFactor)
{
	bound = strtol(theBound, NULL, 10);
	numberOfThreads = strtol(theNumberOfThreads, NULL, 10);

	if (numberOfThreads == 0)
		numberOfThreads = 1;
	if (numberOfThreads > MAX_THREADS)
		numberOfThreads = MAX_THREADS;
}

qfs::~qfs()
{
}


// The QFS algorithm requires a list of prime numbers as an input. There is no need to calculate them more than once,
// hence it is a good idea to keep a cache of them saved to a file. Therefore this algorithm first checks to see if a
// cache file exists and if so loads the primes from the cache file. If there were enough primes in the cache file to
// fulfill the bound requested then the function will fill the primes vector up to the required limit and exit. If not
// it will perform a Sieve of Eratosthenes (276 BC - 194 BC) to calculate the rest of the primes. The format of the
// cache file is just one prime on each line starting with 2 on the first line.
void qfs::getPrimesUpToBound(void)
{
    string currentLine;
    bool fileExists = false;

    if (access(qfs::cacheFilePath, R_OK|W_OK) != -1) {
    	fileExists = true;
    }

    if (fileExists) {
    	primesCacheFile.open(qfs::cacheFilePath, std::fstream::in | std::fstream::out | std::fstream::app);
    	if (primesCacheFile.is_open() == false) {
    		throw std::runtime_error("Could not open prime number cache file for reading.\n");
    	}
    }

    // Read in primes from our cache file
    while (getline(primesCacheFile, currentLine)) {

        uint64_t currentPrime = strtoull(currentLine.c_str(), NULL, 10);

        if (currentPrime <= bound) {

        	// And we should add it
            if (primes.size() == 0 ||                                   // We have no primes, just add it
               (primes.size() > 0 && currentPrime > primes.back())) {   // Only add the prime if new prime is larger
                primes.push_back(currentPrime);
            }
        }
        else {
        	primesCacheFile.close();
            return;
        }
    }
    primesCacheFile.clear();

    // If we get here then there weren't enough primes in the cache file to meet our demands. Perform the modified
    // Sieve Of Eratosthenes as described above, to fill in the rest of the primes.
    sieveOfEratosthenes();

    /*
    std::cout << "List of primes:\n";
    for (auto it = primes.begin(); it != primes.end(); it++) {
        std::cout << *it << "\n";
    }
    */
}

void qfs::sieveFactorBase(void)
{
	int exact = 0;
	mpz_class a;
	mpz_class aMin;
	mpz_class aMax;
	mpz_class twoN;
	mpz_class fOfX;

    // Now that we have the prime factor candidates that make up the factor base, compute the polynomials f(x) and
	// begin sieving by dividing out the prime factor candidates until we find enough smooth numbers. If there are N
	// numbers in the factor base then we need N+1 smooth numbers to guarantee a unique non-zero solution.

	// let:
	//     n = numberToFactor
	//     a = ceil(sqrt(n))		aMin
    //     a' = ceil(sqrt(2n))		aMax
	//     f(x) = x^2 - n
	//
	// for ai = a, a + 1, a + 2, ... , a'
	// calculate f(ai) = ai^2 - n
	//
	// Keep track of pairs: f(ai), ai, and reminder ri
	//
	// f(a1), f(a2), f(a3), f(a4), ... , f(a')
	//   a1 ,   a2 ,   a3 ,   a4 , ... ,   a'
	//   r1 ,   r2 ,   r3 ,   r4 , ... ,   r'
	//
	// Divide f(ai) by each prime p in the factor base, until ri is 1
	//

	// Calculate aMin = ceil(sqrt(n))
	exact =  mpz_root(aMin.get_mpz_t(), numberToFactor.get_mpz_t(), 2);	// Square root
	if (exact == 0) {
		// Because mpz_root truncates, we add 1 to achieve ceil(result)
		aMin += 1;
	}

	// Calculate aMax = ceil(sqrt(2*n))
	//mpz_mul_ui(twoN, numberToFactor, 2);		// Multiply by 2
	twoN = numberToFactor * 2;
	exact = mpz_root(aMax.get_mpz_t(), twoN.get_mpz_t(), 2);
	if (exact == 0) {
		// Because mpz_root truncates, we add 1 to achieve ceil(result)
		aMax += 1;
	}

	// Calculate relations over the bound a to a' and fill the relationsVector
	for (a = aMin; a < aMax; a++) {
		// f(x) = a^2 - n
		mpz_pow_ui(fOfX.get_mpz_t(), a.get_mpz_t(), 2); 	// FIXME: Large numbers better use exponentiaion by squaring
		fOfX -= numberToFactor;
		relationsVector.push_back(relations(fOfX, a));
	}

	cout << "  xi = ";
	int j = 0;
	for (auto i = relationsVector.begin(); i != relationsVector.end(); ++i) {
		cout.width(5);
		cout << i->xi << " ";
		j++;
	}
	cout << "\n";

	cout << "f(xi)= ";
	j = 0;
	for (auto i = relationsVector.begin(); i != relationsVector.end(); ++i) {
		cout.width(5);
		cout << i->f_xi << " ";
		j++;
	}
	cout << "\n";


	// Divide up work into threads


	for (auto p = factorBase.begin(); p != factorBase.end(); ++p) {

		// Print out the current prime
        cout.width(6);
        cout << *p;

        // For each relation we want to see if the factors of f(x) are all small numbers in our factor base
        // i.e. we are searching for B-smooth numbers.
        for (auto rel = relationsVector.begin(); rel != relationsVector.end(); ++rel) {

			mpz_class remainder;
			mpz_class quotient;
			uint64_t exponentFactor = 0;
			bool factorsLeft = true;

			do
			{
				// quotient = remaining / p
				// remainder = remaining % p
				mpz_tdiv_qr_ui(quotient.get_mpz_t(), remainder.get_mpz_t(), rel->remaining.get_mpz_t(), *p);

				// If the remainder is zero
				if (mpz_cmp_ui(remainder.get_mpz_t(), 0) == 0) {

					// We found a factor so count it
					exponentFactor++;

					// Update remaining to what's left after the prime factor is divided out
					rel->remaining = quotient;

					// We should start collecting our data here. What format do we need?
					// 		x
					// 		f(x)
					// 		f(x)_remaining
					// p1	0
					// p2	1
					// p3	0
					// p4	2


					// What is the format here that the matrix library uses?

				}
				else {
					factorsLeft = false;
				}

			} while (factorsLeft == true);

			// Print out the exponent factor
			cout.width(6);
			cout << exponentFactor;
		}
		cout << "\n";
	}

}

void qfs::sieveOfEratosthenes(void)
{
    vector<prime_canidate_t> primeCanidates;
    uint64_t sievingStartPoint = 2;

    // Populate the field of prime candidates that we will sieve from
    for (uint64_t number = 0; number <= bound; number++) {
        prime_canidate_t canidate;
        canidate.number = number;
        canidate.flags = 0;
        primeCanidates.push_back(canidate);
    }

    // Now check if we already have some primes passed in or not, if so we can save some work
    if (primes.size() > 0) {

        // For each prime passed in
        for (auto it = primes.begin(); it != primes.end(); it++) {

            // Mark it as prime
            primeCanidates[*it].flags |= PRIME;

            // Mark it's multiples as composite
            for (uint64_t index = 2 * (*it); index <= bound; index += *it) {
                primeCanidates[index].flags |= COMPOSITE;
            }

            // Update the point where we can start sieving (everything below this has been visited and discovered as
            // either prime or composite.
            if (*it > sievingStartPoint) {
                sievingStartPoint = *it;
            }
        }
    }

    // Set output stream position to end of primes cache file
    //primesCacheFile.seekp(0, std::ios_base::end);

    // Now we begin sieving for primes, but first we have to determine a starting point. If we had no primes in the
    // cache file we just start at 2 like normal, otherwise if we did have primes, we start at the next largest odd
    // number (the prime + 2).
    uint64_t i = (sievingStartPoint > 2) ? sievingStartPoint + 2 : 2;
    for (/* initalized above */; i <= bound; i++) {

        // You can think of primeCanidates[i] as "current", the current number being evaluated.

        // If the number is marked as composite we have already visited it and we can skip it
        if ((primeCanidates[i].flags & COMPOSITE) == COMPOSITE) {
            continue;
        }

        // This number has not been visited so it must be prime
        primeCanidates[i].flags |= PRIME;
        primes.push_back(primeCanidates[i].number);
        primesCacheFile << std::to_string(primeCanidates[i].number) << "\n";


        // Mark it's multiples as composite
        for (uint64_t j = 2 * primeCanidates[i].number; j <= bound; j += primeCanidates[i].number) {
            primeCanidates[j].flags |= COMPOSITE;
        }
    }

    primesCacheFile.close();
}

int qfs::legendreSymbol(mpz_class& n, uint64_t p)
{
	int retval = 0;

	mpz_class result;
	mpz_class base(n); 	// Set the base
	mpz_class mod(p); 	// Set the modulus
	mpz_class exp;

    // Set the exponent
    uint64_t power = (p-1)/2;
    exp = power;

    //
    // LegendreSymbol(N|p) = (N|p) = N^((p-1)/2) (mod p)
    //
    mpz_powm(result.get_mpz_t(), base.get_mpz_t(), exp.get_mpz_t(), mod.get_mpz_t()); // result = (base2^exp) % mod

    //cout << "Legendre(N|p)=" << result << "\n";

    if (result == 0) {
    	retval = 0;
    }
    else if(result == 1) {   // FIXME: does this need to be congruent to 1 mod p ??
    	retval = 1;
    }
    else {
        retval = -1;
    }

    return retval;
}

void qfs::fillFactorBase(void)
{
	for (auto primeItr = primes.begin(); primeItr != primes.end(); ++primeItr) {

        // The Legendre Symbol is undefined for "2", so always add 2 as a factor
        if (*primeItr == 2) {
            factorBase.push_back(*primeItr);
            printf("Added %lu to list of prime factor canidates.\n", *primeItr);
            continue;
        }

        if (legendreSymbol(numberToFactor, *primeItr) == 1) {
            factorBase.push_back(*primeItr);
            printf("Added %lu to list of prime factor canidates.\n", *primeItr);
        }
    }
}

bool qfs::factor(void)
{
    bool factorResult = false;

    getPrimesUpToBound();

    fillFactorBase();

    sieveFactorBase();

    return factorResult;
}

void qfs::printResult(void)
{

}

} // namespace qfs
