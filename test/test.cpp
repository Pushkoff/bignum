#include "BigNum.h"
#include <assert.h>
#include <string>
#include <algorithm>
#include <chrono>
#include <time.h>

template<int N>
std::string toString(const BigNum::Num<N>& num)
{
	std::string ret;
	BigNum::Num<N> Q = num, R;

	while (Q > BigNum::Num<N>(0))
	{
		R = Q % BigNum::Num<N>(10);
		Q = Q / BigNum::Num<N>(10);
		ret.push_back(R[0] + '0');
	}
	std::reverse(ret.begin(), ret.end());
	return ret;
}

template<int N>
const BigNum::Num<N> fromString(const char* str)
{
	BigNum::Num<N> ret;
	while (*str && isdigit(*str))
	{
		ret = ret * 10;
		ret = ret + BigNum::Num<N>(*str - '0');
		str++;
	}
	return ret;
}

template<int N>
const BigNum::Num<N> rand()
{
	BigNum::Num<N> ret;
	for (int i = 0; i < BigNum::Num<N>::Size; ++i)
	{
		ret[i] = std::rand() % 256;
	}
	return ret;
}

// http://en.literateprograms.org/index.php?title=Special:DownloadCode/Miller-Rabin_primality_test_(C)&oldid=18973

template<int N>
bool millerRabinPass(const BigNum::Num<N>& num, const BigNum::Num<N>& a)
{
	BigNum::Num<N> d = num - 1;
	
	int s = 0;
	while (d.bit(s) == false)
		s++;

	d = d >> s;
	BigNum::Num<N> a_to_power = BigNum::modExp(a, d, num);

	if (a_to_power == 1)
		return true;

	for (int i = 0; i < s - 1; i++)
	{
		if (a_to_power == num - 1)
			return true;
		
		BigNum::Num<N * 2> temp = BigNum::mul2N(a_to_power,a_to_power);

		BigNum::Num<N * 2> q;
		BigNum::div(temp, num, q, a_to_power);
	}

	if (a_to_power == num - 1)
		return true;

	return false;
}

template<int N>
bool millerRabinTest(const BigNum::Num<N>& num)
{
	for (int i = 0; i < 10; i++)
	{
		BigNum::Num<N> a;
		do {
			a = rand<N>() % num;
		} while (a == 0);

		if (millerRabinPass(num, a) == false)
			return false;
	}
	return true;
}



template<int N>
const BigNum::Num<N> findPrime(const BigNum::Num<N>& from)
{
	BigNum::Num<N> prime = from | 1;
	
	const unsigned int primes[] = { 3, 5, 7, 11, 13,  17, 19, 23 , 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 , 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251,
};
		//257,    263,    269,    271,    277,    281,
		//283,    293,    307,    311,    313,    317,    331,    337,    347,    349,
		//353,    359,    367,    373,    379,    383,    389,    397,    401,    409,
		//419,    421,    431,    433,    439,    443,    449,    457,    461,    463,
		//467,    479,    487,    491,    499,    503,    509,    521,    523,    541,
		//547,    557,    563,    569,    571,    577,    587,    593,    599,    601,
		//607,    613,    617,    619,    631,    641,    643,    647,    653,    659,
		//661,    673,    677,    683,    691,    701,    709,    719,    727,    733,
		//739,    743,    751,    757,    761,    769,    773,    787,    797,    809,
		//811,    821,    823,    827,    829,    839,    853,    857,    859,    863,
		//877,    881,    883,    887,    907,    911,    919,    929,    937,    941,
		//947,    953,    967,    971,    977,    983,    991,    997, };
	unsigned int rests[sizeof(primes) / sizeof(primes[0])] = { 0 };

	bool simpleTest = true;
	for (int i = 0; i < sizeof(rests) / sizeof(rests[0]); i++)
	{
		rests[i] = prime % primes[i];
		simpleTest &= (rests[i] != 0);
	}

	while (!(simpleTest && millerRabinTest(prime)))
	{
		//unsigned char addToNum = 0;
		//do{
		//	addToNum += 2;
			simpleTest = true;
			for (int i = 0; i < sizeof(rests) / sizeof(rests[0]); i++)
			{
				rests[i] += 2;
				if (rests[i] >= primes[i])
				{
					rests[i] -= primes[i];
					simpleTest &= (rests[i] != 0);
				}
			}
		//} while (simpleTest == false && addToNum < 200);
				
		//prime = prime + addToNum;
		prime = prime + 2;
	}

	return prime;
}

template<int N>
const BigNum::Num<N> randPrime()
{
	BigNum::Num<N> num = rand<N>() | (BigNum::Num<N>(1) | (BigNum::Num<N>(1) << (N - 1)));
	return findPrime(num);
}

#define PROFILING 0

int main()
{
#if PROFILING
	auto prime = findPrime<1024>(BigNum::Num<1024>(1) << 767);
#else
	//srand(time(nullptr));
	assert(BigNum::Num<256>(fromString<128>("118802731")) == fromString<256>("118802731"));

	assert(BigNum::Num<64>(0) == BigNum::Num<64>(0));
	assert(BigNum::Num<64>(0) != BigNum::Num<64>(1));

	assert(BigNum::Num<64>(1) + BigNum::Num<64>(1) == BigNum::Num<64>(2));
	assert(BigNum::Num<64>(2) - BigNum::Num<64>(1) == BigNum::Num<64>(1));

	assert(BigNum::Num<64>(256) - BigNum::Num<64>(1) == BigNum::Num<64>(255));
	assert(BigNum::Num<16>(255<<8) - BigNum::Num<16>(1) == BigNum::Num<16>((255 << 8) -1));

	assert(BigNum::Num<64>(256) - 1 == BigNum::Num<64>(255));
	assert(BigNum::Num<16>(255 << 8) - 1 == BigNum::Num<16>((255 << 8) - 1));

	assert(BigNum::Num<64>(20) * BigNum::Num<64>(45) == BigNum::Num<64>(20 * 45));

	BigNum::Num<1024> fact(1);
	for (unsigned char i = 2; i <= 100; i++)
	{
		fact = fact * i;
	}
	assert(fact == fromString<1024>("93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000"));
	printf("100!=%s\n", toString(fact).c_str());


	assert((BigNum::Num<64>(1) > BigNum::Num<64>(0)));
	assert((BigNum::Num<64>(1) < BigNum::Num<64>(0)) == false);

	assert((BigNum::Num<64>(0) < BigNum::Num<64>(1)));
	assert((BigNum::Num<64>(1) < BigNum::Num<64>(0)) == false);

	assert((BigNum::Num<64>(0) >= BigNum::Num<64>(0)));
	assert((BigNum::Num<64>(1) >= BigNum::Num<64>(0)));
	assert((BigNum::Num<64>(0) >= BigNum::Num<64>(1)) == false);

	assert((BigNum::Num<64>(0) <= BigNum::Num<64>(0)));
	assert((BigNum::Num<64>(0) <= BigNum::Num<64>(1)));
	assert((BigNum::Num<64>(1) <= BigNum::Num<64>(0)) == false);

	assert((BigNum::Num<64>(1 << 8) == BigNum::shift_left(BigNum::Num<64>(1), 8)));
	assert((BigNum::Num<64>(1 << 7) == BigNum::shift_left(BigNum::Num<64>(1), 7)));
	assert((BigNum::Num<64>(135 << 12) == BigNum::shift_left(BigNum::Num<64>(135), 12)));
	assert((BigNum::Num<64>(135056 << 12) == BigNum::shift_left(BigNum::Num<64>(135056), 12)));
	assert((BigNum::Num<64>(135056 << 7) == BigNum::shift_left(BigNum::Num<64>(135056), 7)));
	assert((BigNum::Num<64>(13056 << 15) == BigNum::shift_left(BigNum::Num<64>(13056), 15)));

	assert((BigNum::Num<64>(1) == BigNum::shift_right(BigNum::Num<64>(1<<8), 8)));
	assert((BigNum::Num<64>(135056 >> 1) == BigNum::shift_right(BigNum::Num<64>(135056), 1)));
	assert((BigNum::Num<64>(135056 >> 4) == BigNum::shift_right(BigNum::Num<64>(135056), 4)));
	assert((BigNum::Num<64>(135056 >> 12) == BigNum::shift_right(BigNum::Num<64>(135056), 12)));
	assert((BigNum::Num<64>(13505675 >> 15) == BigNum::shift_right(BigNum::Num<64>(13505675), 15)));
	assert((BigNum::Num<64>(13505675 >> 22) == BigNum::shift_right(BigNum::Num<64>(13505675), 22)));
	assert(BigNum::Num<1024>(1) << 1023 == fromString<1024>("89884656743115795386465259539451236680898848947115328636715040578866337902750481566354238661203768010560056939935696678829394884407208311246423715319737062188883946712432742638151109800623047059726541476042502884419075341171231440736956555270413618581675255342293149119973622969239858152417678164812112068608"));

	//{
	//	unsigned int num = 58369;
	//	BigNum::Num<16> bn(num);
	//	for (int i = 0; i < 16; i++)
	//	{
	//		bn = bn >> 1;
	//		num = num >> 1;
	//		assert(bn == num);
	//	}
	//}

	assert(((BigNum::Num<64>(1) << 63) + 1) / (BigNum::Num<64>(1) << 63) == 1);
	assert(((BigNum::Num<64>(1) << 63) + 1) % (BigNum::Num<64>(1) << 63) == 1);
	assert(((BigNum::Num<64>(1) << 63) - 1) / (BigNum::Num<64>(1) << 63) == 0);
	assert(((BigNum::Num<64>(1) << 63) - 1) % (BigNum::Num<64>(1) << 63) == ((BigNum::Num<64>(1) << 63) - 1));

	assert(BigNum::Num<64>(524568123) / BigNum::Num<64>(45654) == BigNum::Num<64>(524568123 / 45654));
	assert(BigNum::Num<64>(524568123) % BigNum::Num<64>(45654) == BigNum::Num<64>(524568123 % 45654));

	assert(BigNum::Num<64>(524568123) / static_cast<unsigned char>(10) == BigNum::Num<64>(524568123 / 10));
	assert(BigNum::Num<64>(524568123) / static_cast<unsigned short>(10) == BigNum::Num<64>(524568123 / 10));
	assert(BigNum::Num<64>(524568123) / 10u == BigNum::Num<64>(524568123 / 10));
	assert(BigNum::Num<64>(524568123) % static_cast<unsigned char>(10) == 524568123 % 10);
	assert(BigNum::Num<64>(524568123) % static_cast<unsigned short>(10) == 524568123 % 10);
	assert(BigNum::Num<64>(524568123) % 10u == 524568123 % 10);

	printf("BigNum::Num<128>(1025445861) = %s\n", toString(BigNum::Num<128>(1025445861ull)).c_str());

	printf("BigNum::Num<128>(125632587586954854585354458654754812345678908765678765467543456765678) = %s\n", toString(fromString<1024>("125632587586954854585354458654754812345678908765678765467543456765678")).c_str());

	assert(BigNum::modExp<32>(BigNum::Num<32>(4), BigNum::Num<32>(13), BigNum::Num<32>(497)) == BigNum::Num<32>(445));

	millerRabinPass(fromString<128>("118802731"), fromString<128>("74812"));

	srand(0);
		
	assert(millerRabinTest(fromString<128>("961748941")) == true);
	assert(millerRabinTest(fromString<128>("982451653")) == true);
	assert(millerRabinTest(fromString<128>("89685068870671644428994813094690803553")) == true);
	assert(millerRabinTest(fromString<128>("982451655")) == false);


	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = findPrime<64>(BigNum::Num<64>(1) << 63);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime (64)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = findPrime<128>(BigNum::Num<128>(1) << 63);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime (128)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = findPrime<256>(BigNum::Num<256>(1) << 63);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime (256)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = findPrime<512>(BigNum::Num<512>(1) << 63);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime (512)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
		{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = findPrime<1024>(BigNum::Num<1024>(1) << 63);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime (1024)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = findPrime<128>(BigNum::Num<128>(1) << 127);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime (128)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = findPrime<256>(BigNum::Num<256>(1) << 255);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime (256)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = findPrime<512>(BigNum::Num<512>(1) << 511);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime (512)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = findPrime<1024>(BigNum::Num<1024>(1) << 1023);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime (1024)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = findPrime<2048>(BigNum::Num<2048>(1) << 2047);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime (2048)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	printf("Press any key...");
	getchar();
#endif
	return 0;
}