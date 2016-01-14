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
		
		BigNum::Num<N * 2 + 1> temp = BigNum::Num<N*2 + 1>(a_to_power) * BigNum::Num<N*2 + 1>(a_to_power);

		a_to_power = BigNum::Num<N>(temp % BigNum::Num<N * 2 + 1>(num));
	}

	if (a_to_power == num - 1)
		return true;

	return false;
}

template<int N>
bool millerRabinTest(const BigNum::Num<N>& num)
{
	for (int i = 0; i < 20; i++)
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
	

	const unsigned char primes[] = { 3, 5, 7, 11, 13,  17, 19, 23 , 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251 };
	unsigned short rests[sizeof(primes) / sizeof(primes[0])] = { 0 };

	for (int i = 0; i < sizeof(rests) / sizeof(rests[0]);i++)
		rests[i] = prime % primes[i];

	bool found = false;
	while (found == false)
	{
		bool simpleTest = false;
		for (int i = 0; i < sizeof(rests) / sizeof(rests[0]);i++)
		{
			if (rests[i] == 0)
				simpleTest = true;
			
			rests[i] += 2;
			if (rests[i] >= primes[i])
				rests[i] -= primes[i];
		}

		if (simpleTest == false && millerRabinTest(prime))
			break;
				
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


int main()
{
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

	assert(BigNum::Num<64>(524568123) / BigNum::Num<64>(45654) == BigNum::Num<64>(524568123 / 45654));
	assert(BigNum::Num<64>(524568123) % BigNum::Num<64>(45654) == BigNum::Num<64>(524568123 % 45654));

	assert(BigNum::Num<64>(524568123) / 10 == BigNum::Num<64>(524568123 / 10));
	assert(BigNum::Num<64>(524568123) % 10 == 524568123 % 10);

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
		printf("Prime = %s\n   duration - %lld ms\n", toString(prime).c_str(), std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = findPrime<128>(BigNum::Num<128>(1) << 127);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime = %s\n   duration - %lld ms\n", toString(prime).c_str(), std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = findPrime<256>(BigNum::Num<256>(1) << 255);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime = %s\n   duration - %lld ms\n", toString(prime).c_str(), std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = findPrime<512>(BigNum::Num<512>(1) << 511);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime = %s\n   duration - %lld ms\n", toString(prime).c_str(), std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	printf("Press any key...");
	getchar();
	return 0;
}