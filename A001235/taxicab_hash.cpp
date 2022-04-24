#include<iostream>
#include<fstream>
#include<cstdlib>
#include<unordered_set>
#include<algorithm>
#include<cmath>

using namespace std;

/**
 * Enumerate all (a,b) -> a^3 + b^3 with b <= a
 * Store in a hash set (unordered_set) occasionally removing
 * when the hashset's size grows large
 */
void findHardyNumbers(long maxA)
{
	double oneThird = ((double) 1) / ((double) 3);
	long maxA3B3 = maxA * maxA * maxA;
	unordered_set<long> taxiCabNumbers;
	unordered_set<long> a3b3;
	for (long a = 0; ; a ++)
	{
		long cubeA = a * a * a;
		long nextCubeA = (a+1) * (a+1) * (a+1);

		if (cubeA > maxA3B3)
			break;

        // This should maybe grow if percent removed isn't large enough.
		if (a3b3.size() > 5'000'000)
		{
			cout << "a3b3 has to many items, at a = " << a << endl;
            // Remove any items below a^3
			int count = 0;
			unordered_set<long>::iterator iter;
			for (iter = a3b3.begin(); iter != a3b3.end(); ) {
				if ((*iter) < cubeA) {
					iter = a3b3.erase(iter);
					count += 1;
				} else {
                    iter++;
                }
            }
            printf("\tremoved %d elements, new size: %lu\t removed: %.1f\n",
                count, a3b3.size(), (float) count / (count + a3b3.size()));
		}

		long maxB = ((cubeA + cubeA) < maxA3B3) ? a : (long) pow((double) maxA3B3 - cubeA, oneThird);

		for (long b = 0; b <= maxB; b ++)
		{
			long temp = cubeA + b * b * b;
			if (a3b3.find(temp) != a3b3.end())
				taxiCabNumbers.insert(temp);
			else
				if (temp > nextCubeA)
					a3b3.insert(temp);
		}
	}

    cout << "\tmaxA: " << maxA << "\ta3b3 elements: " << a3b3.size() << endl;

	a3b3.clear();

	long taxiCabArray[taxiCabNumbers.size()];
	int length = 0;

	unordered_set<long>::iterator iter;
	for (iter = taxiCabNumbers.begin(); iter != taxiCabNumbers.end(); iter ++)
		taxiCabArray[length++] = (*iter);

	taxiCabNumbers.clear();

	sort(taxiCabArray, taxiCabArray + length);

	cout << endl << endl << "SORTED all " << length << " taxicab numbers." << endl;
}

long nthHardyNumber(long N)
{
	long maxA3B3 = 10000LL * 10000LL * 10000LL;

	unordered_set<long> taxiCabNumbers;
	unordered_set<long> a3b3;

	for (long a = 0; ; a++)
	{
		long cubeA = a * a * a;
		long nextCubeA = (a+1) * (a+1) * (a+1);

		if (cubeA > maxA3B3)
			break;

		if (a3b3.size() > 8'000'000)
		{
			cout << "a3b3 has to many items, at a = " << a << endl;
			int count = 0;
			unordered_set<long>::iterator iter;
			for (iter = a3b3.begin(); iter != a3b3.end(); ) {
				if ((*iter) < cubeA) {
					iter = a3b3.erase(iter);
					count += 1;
				} else {
                    iter++;
                }
            }
            printf("\tremoved %d elements, new size: %lu\t removed: %.1f\n",
                count, a3b3.size(), (float) count / (count + a3b3.size()));
		}

		for (long b = 0; b <= a; b ++)
		{
			long cubeB = b * b * b;

			long temp = cubeA + cubeB;

			if (temp > maxA3B3)
				break;

			if (a3b3.find(temp) != a3b3.end())
			{
				taxiCabNumbers.insert(temp);

				int found = taxiCabNumbers.size();
				if (found == N)
				{
					unordered_set<long>::iterator iter;

					maxA3B3 = 0;
					for (iter = taxiCabNumbers.begin(); iter != taxiCabNumbers.end(); iter ++)
						if ((*iter) > maxA3B3)
							maxA3B3 = (*iter);
					cout << endl << "FOUND N MANY, max is " << maxA3B3 << endl << endl;
				}
			}
			else
				if (temp > nextCubeA)
					a3b3.insert(temp);
		}
	}

	a3b3.clear();

	long taxiCabArray[taxiCabNumbers.size()];
	int length = 0;

	unordered_set<long>::iterator iter;
	for (iter = taxiCabNumbers.begin(); iter != taxiCabNumbers.end(); iter ++)
		taxiCabArray[length++] = (*iter);

	taxiCabNumbers.clear();

	sort(taxiCabArray, taxiCabArray + length);

	cout << endl << endl << "SORTED all " << length << " taxicab numbers." << endl;

	return taxiCabArray[N-1];
}


int main(int argc, char** argv)
{
    int maxA = 10000;
	findHardyNumbers(maxA);

    // Slower? but doesn't require guestimate on A.
	int n = 40000;
	long nth = nthHardyNumber(n);
    cout << "n" << n << ": " << nth << endl;

    /*
	if (argc != 2)
	{
		cout << "takes one paramater (output file)" << endl;
		return -1;
	}
    */

//	ofstream outFile;
//	outFile.open(fileName);
//	for (int index = 0; index < length; index ++)
//	{
//		cout << (index + 1) << ": " << taxiCabArray[index] << endl;
//		outFile << (index + 1) << ": " << taxiCabArray[index] << endl;
//	}
//	outFile.close();

	return 0;
}
