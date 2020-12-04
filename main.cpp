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
#include <iostream>
#include <ctype.h>
#include "qfs.h"

using std::cout;

void usage(char*);
bool isInteger(const char*);

int main(int argc, char* argv[])
{
	if (argc != 3) {
		usage(argv[0]);
		exit(1);
	}

	if (isInteger(argv[1]) == false ||
		isInteger(argv[2]) == false ||
		isInteger(argv[3]) == false) {
		usage(argv[0]);
		exit(1);
	}

    qfs::qfs quadraticFieldSieve(argv[1], argv[2], argv[3]);
    quadraticFieldSieve.factor();
    quadraticFieldSieve.printResult();

	return 0;
}

void usage(char* appName)
{
	cout << appName << " [Number To Factor] [Bound] [Number of Threads]\n";
}

// Can't use strtol since we're expecting numbers > 64 bits
bool isInteger(const char* input)
{
	for (const char* c = &input[0]; *c != '\0'; c++) {
		if (isdigit(*c) != true)
			return false;
	}

	return true;
}
