#include <assert.h>
#include <string>
#include <string.h>
#include <algorithm>
#include <chrono>
#include <vector>
#include <time.h>
#include <stdint.h>
#include "BigNum.h"

template<size_t N>
std::string toString(BigNum::Num<N> num)
{
	std::string ret;

	while (num > 0)
	{
		auto R = num % 10;
		assert(R < 10);
		num = num / 10;
		ret.push_back(char(R + '0'));
	}
	std::reverse(ret.begin(), ret.end());
	if (ret.empty())
		ret = "0";
	return ret;
}

template<size_t N>
std::string toHex(BigNum::Num<N> num)
{
	const char Nums[] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' };
	std::string ret;
	while (num > 0)
	{
		BigNum::Word R = num[0] & 0xF;
		ret.push_back(Nums[R]);
		num = num >> 4;
	}
	ret.push_back('x');
	ret.push_back('0');
	std::reverse(ret.begin(), ret.end());
	return ret;
}

template<size_t N>
std::string toRaw(BigNum::Num<N> num)
{
	std::string ret = "{ ";
	for (size_t i = 0; i < BigNum::Num<N>::Size; ++i)
	{
		if (i != 0)
			ret += ", ";

		char buf[32] = {0};
		snprintf(buf, sizeof(buf), "%lluull", num[i]);
		ret += buf;
	}
	ret += " }";
	return ret;
}

template<size_t N>
BigNum::Num<N> fromString(const char* str)
{
	BigNum::Num<N> ret;
	while (*str && isdigit(*str))
	{
		ret = ret * 10u;
		ret = ret + (BigNum::Word(*str) - '0');
		str++;
	}
	return ret;
}

template<size_t N>
BigNum::Num<N> fromHex(const char* str)
{
	auto hexToNum = [=](const char a)->int
	{
		if (a >= '0' && a <= '9')
			return (a - '0');
		if (a >= 'A' && a <= 'F')
			return a - 'A' + 10;
		if (a >= 'a' && a <= 'f')
			return a - 'a' + 10;
		return -1;
	};

	BigNum::Num<N> ret;
	while (*str && hexToNum(*str) != -1)
	{
		ret = ret << 4;
		ret[0] = ret[0] + hexToNum(*str);
		str++;
	}
	return ret;
}

int TestsPass = 0;
int TestsFailed = 0;

template<typename Fn>
void doTest(Fn fn, const char* testName, int fileLine = 0)
{
	srand(0);
	printf("%3d  %s\n", fileLine, testName);
	auto start = std::chrono::high_resolution_clock::now();
	bool ret = fn();
	auto stop = std::chrono::high_resolution_clock::now();
	auto elapsed = stop - start;
	(ret ? TestsPass : TestsFailed)++;
	printf("   %s - duration - %lld ms\n", ret ? "Ok" : "Fail", (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());

	if (!ret) throw(1);
}

template<size_t N>
BigNum::Num<N> biggestPrimeProduct() noexcept
{
	BigNum::Num<N + BigNum::kWordSizeBits> mask = 1;
	mask = mask << N;
	BigNum::Num<N + BigNum::kWordSizeBits> test = 1;
	BigNum::Num<N> ret = 0;

	for (size_t i = 0; i < BigNum::primesCount; ++i)
	{
		test = test * BigNum::primes[i];
		if (test < mask)
			ret = test;
		else
			break;
	}
	return ret;
}


#define test(a) do{  doTest([&](){ return (a); }, #a, __LINE__); }while(false);

#define PROFILING 0

int main()
{
#if PROFILING
	test(BigNum::nextPrimeOpt357<2048>(1_bn2048 << 2047) > 0);
	test(BigNum::nextPrime<2048>(1_bn2048 << 2047) > 0);
	test(BigNum::nextPrimeOpt357<2048>(1_bn2048 << 2047) > 0);
#else
	//printf("biggestPrimeProduct<64>() = %s\n", toRaw(biggestPrimeProduct<64>()).c_str());
	//printf("biggestPrimeProduct<128>() = %s\n", toRaw(biggestPrimeProduct<128>()).c_str());
	//printf("biggestPrimeProduct<256>() = %s\n", toRaw(biggestPrimeProduct<256>()).c_str());
	//printf("biggestPrimeProduct<512>() = %s\n", toRaw(biggestPrimeProduct<512>()).c_str());
	//printf("biggestPrimeProduct<1024>() = %s\n", toRaw(biggestPrimeProduct<1024>()).c_str());
	//printf("biggestPrimeProduct<2048>() = %s\n", toRaw(biggestPrimeProduct<2048>()).c_str());

	//srand(time(nullptr));
	try
	{
		test((BigNum::Num<128>(1) << 1) == BigNum::Num<128>(2));

		{
			bool passed = true;
			BigNum::Num<128> bn(1);
			bn <<= 64;
			for (size_t j = 0; j < BigNum::Num<128>::Size; j++)
			{
				if (j == 64 / BigNum::kWordSizeBits)
					passed = passed && (bn[j] == 1);
				else
					passed = passed && (bn[j] == 0);
			}
			test(passed == true && "Left shifts");
		}

		{
			bool passed = true;
			BigNum::Num<128> bn(BigNum::Word(1) << (BigNum::kWordSizeBits - 1));
			bn <<= 1;
			for (size_t j = 0; j < BigNum::Num<128>::Size; j++)
			{
				if (j == 64 / BigNum::kWordSizeBits)
					passed = passed && (bn[j] == 1);
				else
					passed = passed && (bn[j] == 0);
			}
			test(passed == true && "Left shifts");
		}

		{
			BigNum::Num<256> bn = 1_bn256 << 192;
			bn <<= 1;
			test(bn == (1_bn256 << 193));
		}

		{
			BigNum::Num<256> bn = 1_bn256 << 64;
			bn >>= 1;
			test(bn == (1_bn256 << 63));
		}

		{
			BigNum::Num<256> bn = 1_bn256 << 63;
			bn >>= 1;
			test(bn == (1_bn256 << 62));
		}

		{
			bool passed = true;
			BigNum::Num<128> bn(1);
			for (size_t i = 0u; i < 128u; i++)
			{
				BigNum::Num<128> testbn = bn << i;
				for (size_t j = 0; j < BigNum::Num<128>::Size; j++)
				{
					if (j == i / BigNum::kWordSizeBits)
						passed = passed && (testbn[j] == (1ull << (i%BigNum::kWordSizeBits)));
					else
						passed = passed && (testbn[j] == 0);
				}
				if (!passed)
				{
					printf("%zd\n", i);
					break;
				}
			}
			test(passed == true && "Left shifts");
		}

		{
			bool passed = true;
			BigNum::Num<128> bn(0);
			bn[BigNum::Num<128>::Size - 1] = (BigNum::Word(1) << (BigNum::kWordSizeBits - 1));

			for (size_t i = 0u; i < 128u; i++)
			{
				BigNum::Num<128> testbn = bn;
				testbn >>= i;
				for (size_t j = BigNum::Num<128>::Size; j --> 0; )
				{
					if (j == (BigNum::Num<128>::Size - 1 - i / BigNum::kWordSizeBits))
						passed = passed && (testbn[j] == (1ull << (BigNum::kWordSizeBits - 1 - i % BigNum::kWordSizeBits)));
					else
						passed = passed && (testbn[j] == 0);
				}
				if (!passed)
				{
					printf("%zd\n", i);
					break;
				}
			}
			test(passed == true && "Right shifts");
		}

		{
			bool passed = true;
			BigNum::Num<64> bn = BigNum::Num<64>(1) << 63u;
			for (unsigned int i = 0u; i < 64u; i++)
			{
				BigNum::Num<64> testbn = bn >> i;
				passed = passed && (testbn[(63 - i) / BigNum::kWordSizeBits] == (1ull << ((BigNum::kWordSizeBits - 1) - i%BigNum::kWordSizeBits)));
				assert(passed);
			}
			test(passed == true && "Rignt shift");
		}

		{
			bool passed = true;
			BigNum::Num<256> bn(1);
			for (unsigned int i = 1u; i < 256; i++)
			{
				bn <<= 1;
				passed = passed && (bn[(i) / BigNum::kWordSizeBits] == (1ull << (i % BigNum::kWordSizeBits)));
				assert(passed);
			}
			test(passed == true && "Left shift");
		}

		{
			bool passed = true;
			BigNum::Num<256> bn = BigNum::Num<256>(1) << 255u;
			for (unsigned int i = 1u; i < 256; i++)
			{
				bn >>= 1;
				passed = passed && (bn[(255 - i) / BigNum::kWordSizeBits] == (1ull << ((BigNum::kWordSizeBits - 1) - i % BigNum::kWordSizeBits)));
				assert(passed);
			}
			test(passed == true && "Rignt shift");
		}

		{
			bool passed = true;
			BigNum::Num<64> bn1 = BigNum::Num<64>(0) - 1u;
			for (unsigned int i = 0; i < 64u; i++)
			{
				BigNum::Num<64> testbn1 = bn1 << i;
				passed = passed && (testbn1 == (BigNum::Num<64>(0) - BigNum::Num<64>(1) << i));
				assert(passed);
			}
			test(passed == true && "Left shifts");
		}

		{
			bool passed = true;
			BigNum::Num<64> bn1 = BigNum::Num<64>(255);
			BigNum::Num<64> bn2 = BigNum::Num<64>(255) << 56;
			for (unsigned int i = 0; i < 56; i++)
			{
				BigNum::Num<64> testbn1 = bn1 << i;
				BigNum::Num<64> testbn2 = bn2 >> (56 - i);
				passed = passed && (testbn1 == testbn2);
				assert(passed);
			}
			test(passed == true && "Both shifts");
		}

		{
			bool passed = true;
			BigNum::Num<64> bn(3);
			for (unsigned int i = 0; i < 64; i++)
			{
				BigNum::Num<64> testbn = BigNum::Num<64>(3) << i;
				passed = passed && (testbn == bn);

				bn = bn * 2u;
			}
			test(passed == true && "Left shifts");
		}

		test(BigNum::Num<256>(fromString<128>("118802731")) == fromString<256>("118802731"));
		test(BigNum::Num<256>(fromString<128>("118802731")) == 118802731_bn2048);
		test(BigNum::Num<256>(fromString<128>("118802731")) == 118802731_bn1024);
		test(fromString<128>(toString(fromString<128>("118802731")).c_str()) == fromString<128>("118802731"));

		test((1_bn128 << 127) / 10_bn128 == 17014118346046923173168730371588410572_bn128);
		test((BigNum::Num<128>(1) << 127) % BigNum::Num<128>(10) == BigNum::Num<128>(8));


		test(toString(fromString<128>("118802731")) == "118802731");
		test(toString(BigNum::Num<128>(1) << 127) == "170141183460469231731687303715884105728");
		test(toString(BigNum::Num<256>(1) << 255) == "57896044618658097711785492504343953926634992332820282019728792003956564819968");
		test(toString(1_bn1024 << 1023) == "89884656743115795386465259539451236680898848947115328636715040578866337902750481566354238661203768010560056939935696678829394884407208311246423715319737062188883946712432742638151109800623047059726541476042502884419075341171231440736956555270413618581675255342293149119973622969239858152417678164812112068608");

		test(BigNum::Num<64>(0) == BigNum::Num<64>(0));
		test(BigNum::Num<64>(0) != BigNum::Num<64>(1));

		test(BigNum::Num<64>(1) + BigNum::Num<64>(1) == BigNum::Num<64>(2));
		test(BigNum::Num<64>(2) - BigNum::Num<64>(1) == BigNum::Num<64>(1));

		{
			BigNum::Num<128> x(1);
			x += BigNum::Num<128>(1);
			test(x == BigNum::Num<128>(2));
			x -= BigNum::Num<128>(1);
			test(x == BigNum::Num<128>(1));
			x += 1;
			test(x == BigNum::Num<128>(2));
			x -= 1;
			test(x == BigNum::Num<128>(1));
		}

		{
			BigNum::Num<128> x(1);
			x = (x << 64) - 1;
			x += BigNum::Num<128>(1);
			test(x == (BigNum::Num<128>(1) << 64));
			x -= BigNum::Num<128>(1);
			test(x == ((BigNum::Num<128>(1) << 64) - 1));
			x += 1;
			test(x == (BigNum::Num<128>(1) << 64));
			x -= 1;
			test(x == ((BigNum::Num<128>(1) << 64) - 1));
		}

		test(BigNum::Num<64>(256) - BigNum::Num<64>(1) == BigNum::Num<64>(255));
		test(BigNum::Num<16>(255 << 8) - BigNum::Num<16>(1) == BigNum::Num<16>((255 << 8) - 1));

		test(BigNum::Num<64>(256) - 1 == BigNum::Num<64>(255));
		test(BigNum::Num<16>(255 << 8) - 1 == BigNum::Num<16>((255 << 8) - 1));

		test(BigNum::Num<64>(20) * BigNum::Num<64>(45) == BigNum::Num<64>(20 * 45));
		test(BigNum::Num<32>(0xFFFFFFFFu) * BigNum::Num<32>(0xFFFFFFFFu) == BigNum::Num<64>(0xFFFFFFFE00000001ull));
		test(BigNum::Num<64>(0xFFFFFFFFFFFFFFFFull) * BigNum::Num<64>(0xFFFFFFFFFFFFFFFFull) == fromHex<128>("FFFFFFFFFFFFFFFE0000000000000001"));

		BigNum::Num<1024> fact(1);
		for (BigNum::Word i = 2; i <= 100; i++)
		{
			fact = fact * i;
		}
		test(toString(fact) == "93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000");
		test(toHex(fact) == "0x1B30964EC395DC24069528D54BBDA40D16E966EF9A70EB21B5B2943A321CDF10391745570CCA9420C6ECB3B72ED2EE8B02EA2735C61A000000000000000000000000");
		test(fact == 93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000_bn1024);
		test(toString(fromString<1024>("93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000")) == "93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000");


		test((BigNum::Num<64>(1) > BigNum::Num<64>(0)));
		test((BigNum::Num<64>(1) < BigNum::Num<64>(0)) == false);

		test((BigNum::Num<64>(0) < BigNum::Num<64>(1)));
		test((BigNum::Num<64>(1) < BigNum::Num<64>(0)) == false);

		test((BigNum::Num<64>(0) >= BigNum::Num<64>(0)));
		test((BigNum::Num<64>(1) >= BigNum::Num<64>(0)));
		test((BigNum::Num<64>(0) >= BigNum::Num<64>(1)) == false);

		test((BigNum::Num<64>(0) <= BigNum::Num<64>(0)));
		test((BigNum::Num<64>(0) <= BigNum::Num<64>(1)));
		test((BigNum::Num<64>(1) <= BigNum::Num<64>(0)) == false);

		test((BigNum::Num<64>(1) == BigNum::Num<64>(1) << 0));
		test((BigNum::Num<64>(1 << 8) == BigNum::Num<64>(1) << 8));
		test((BigNum::Num<64>(1 << 7) == BigNum::Num<64>(1) << 7));
		test((BigNum::Num<64>(135 << 12) == BigNum::Num<64>(135) << 12));
		test((BigNum::Num<64>(135056 << 12) == BigNum::Num<64>(135056) << 12));
		test((BigNum::Num<64>(135056 << 7) == BigNum::Num<64>(135056) << 7));
		test((BigNum::Num<64>(13056 << 15) == BigNum::Num<64>(13056) << 15));

		test((BigNum::Num<64>(1) == BigNum::Num<64>(1) >> 0));
		test((BigNum::Num<64>(1) == BigNum::Num<64>(1 << 8) >> 8));
		test((BigNum::Num<64>(135056 >> 1) == BigNum::Num<64>(135056) >> 1));
		test((BigNum::Num<64>(135056 >> 4) == BigNum::Num<64>(135056) >> 4));
		test((BigNum::Num<64>(135056 >> 12) == BigNum::Num<64>(135056) >> 12));
		test((BigNum::Num<64>(13505675 >> 15) == BigNum::Num<64>(13505675) >> 15));
		test((BigNum::Num<64>(13505675 >> 22) == BigNum::Num<64>(13505675) >> 22));
		test(toString(89884656743115795386465259539451236680898848947115328636715040578866337902750481566354238661203768010560056939935696678829394884407208311246423715319737062188883946712432742638151109800623047059726541476042502884419075341171231440736956555270413618581675255342293149119973622969239858152417678164812112068608_bn1024) == "89884656743115795386465259539451236680898848947115328636715040578866337902750481566354238661203768010560056939935696678829394884407208311246423715319737062188883946712432742638151109800623047059726541476042502884419075341171231440736956555270413618581675255342293149119973622969239858152417678164812112068608")
		test((BigNum::Num<1024>(1) << 1023) == 89884656743115795386465259539451236680898848947115328636715040578866337902750481566354238661203768010560056939935696678829394884407208311246423715319737062188883946712432742638151109800623047059726541476042502884419075341171231440736956555270413618581675255342293149119973622969239858152417678164812112068608_bn1024);

		test(((BigNum::Num<64>(1) << 63) + 1) / (BigNum::Num<64>(1) << 63) == 1);
		test(((BigNum::Num<64>(1) << 63) + 1) % (BigNum::Num<64>(1) << 63) == 1);
		printf("%s\n", toHex(((BigNum::Num<64>(1) << 63) - 1)).c_str());
		test(((BigNum::Num<64>(1) << 63) - 1) / (BigNum::Num<64>(1) << 63) == 0);
		test(((BigNum::Num<64>(1) << 63) - 1) % (BigNum::Num<64>(1) << 63) == ((BigNum::Num<64>(1) << 63) - 1));

		test(BigNum::Num<64>(524568123) / BigNum::Num<64>(45654) == BigNum::Num<64>(524568123 / 45654));
		test(BigNum::Num<64>(524568123) % BigNum::Num<64>(45654) == BigNum::Num<64>(524568123 % 45654));

		test(BigNum::Num<64>(524568123) / static_cast<unsigned char>(10) == BigNum::Num<64>(524568123 / 10));
		test(BigNum::Num<64>(524568123) / static_cast<unsigned short>(10) == BigNum::Num<64>(524568123 / 10));
		test(BigNum::Num<64>(524568123) / 10u == BigNum::Num<64>(524568123 / 10));
		test(BigNum::Num<64>(524568123) % static_cast<unsigned char>(10) == 524568123 % 10);
		test(BigNum::Num<64>(524568123) % static_cast<unsigned short>(10) == 524568123 % 10);
		test(BigNum::Num<64>(524568123) % 10u == 524568123 % 10);

		printf("BigNum::Num<128>(1025445861) = %s\n", toString(BigNum::Num<128>(1025445861ull)).c_str());

		printf("BigNum::Num<128>(125632587586954854585354458654754812345678908765678765467543456765678) = %s\n", toString(fromString<1024>("125632587586954854585354458654754812345678908765678765467543456765678")).c_str());

		test(BigNum::modExp<32>(BigNum::Num<32>(4), BigNum::Num<32>(13), BigNum::Num<32>(497)) == BigNum::Num<32>(445));
		test(BigNum::modExp<32>(BigNum::Num<32>(4), BigNum::Num<32>(0), BigNum::Num<32>(497)) == BigNum::Num<32>(1));
		test(BigNum::modExp<32>(BigNum::Num<32>(4), BigNum::Num<32>(1), BigNum::Num<32>(497)) == BigNum::Num<32>(4));

		test(BigNum::gcd(fromString<128>("12"), fromString<128>("9")) == fromString<128>("3"));
		test(BigNum::gcd(fromString<128>("12"), fromString<128>("12")) == fromString<128>("12"));
		test(BigNum::gcd(fromString<128>("12"), fromString<128>("24")) == fromString<128>("12"));
		test(BigNum::gcd(fromString<128>("961748941"), fromString<128>("982451653")) == fromString<128>("1"));

		test(BigNum::lcm(fromString<128>("3"), fromString<128>("4")) == fromString<256>("12"));
		test(BigNum::lcm(fromString<128>("6"), fromString<128>("4")) == fromString<64>("12"));

		test((BigNum::modInv(BigNum::Num<128>(7), BigNum::Num<128>(11)) * BigNum::Num<128>(7)) % BigNum::Num<128>(11) == BigNum::Num<128>(1));
		test((BigNum::modInv(BigNum::Num<128>(961748941), BigNum::Num<128>(982451653)) * BigNum::Num<128>(961748941)) % BigNum::Num<128>(982451653) == 1);

		//BigNum::MonMul<16> mod11(BigNum::Num<16>(121));
		//test(mod11(BigNum::Num<16>(100), BigNum::Num<16>(100)) == BigNum::Num<16>(78));
		test(BigNum::monModMul<64>(BigNum::Num<64>(43), BigNum::Num<64>(56), BigNum::Num<64>(97)) == BigNum::Num<64>(80));
		//test(BigNum::monModExp<32>(BigNum::Num<32>(4), BigNum::Num<32>(13), BigNum::Num<32>(497)) == BigNum::Num<32>(445));
		//test(BigNum::monModExp<32>(BigNum::Num<32>(4), BigNum::Num<32>(0), BigNum::Num<32>(497)) == BigNum::Num<32>(1));
		//test(BigNum::monModExp<32>(BigNum::Num<32>(4), BigNum::Num<32>(1), BigNum::Num<32>(497)) == BigNum::Num<32>(4));
		test(BigNum::monModExp<64>(BigNum::Num<64>(4), BigNum::Num<64>(0), BigNum::Num<64>(497)) == BigNum::Num<64>(1));
		test(BigNum::monModExp<64>(BigNum::Num<64>(4), BigNum::Num<64>(13), BigNum::Num<64>(497)) == BigNum::Num<64>(445));
		test(BigNum::monModExp<64>(BigNum::Num<64>(4), BigNum::Num<64>(1), BigNum::Num<64>(497)) == BigNum::Num<64>(4));

		test(BigNum::modExp<64>(fromString<64>("8814054284918744181"), fromString<64>("1152921504606846977"), fromString<64>("9223372036854775817")) == fromString<64>("4724919961687496308"));
		test(BigNum::monModExp<64>(fromString<64>("8814054284918744181"), fromString<64>("1152921504606846977"), fromString<64>("9223372036854775817")) == fromString<64>("4724919961687496308"));
		test(BigNum::monModExp2ary<64>(fromString<64>("8814054284918744181"), fromString<64>("1152921504606846977"), fromString<64>("9223372036854775817")) == fromString<64>("4724919961687496308"));

		//millerRabinPass(fromString<128>("118802731"), fromString<128>("74812"));

		test(millerRabinTest(fromString<32>("17")) == true);

		test(millerRabinTest(fromString<128>("961748941")) == true);
		test(millerRabinTest(fromString<128>("982451653")) == true);
		test(millerRabinTest(fromString<256>("89685068870671644428994813094690803553")) == true);
		test(millerRabinTest(fromString<128>("982451655")) == false);

		test(BigNum::nextPrime<64>(BigNum::Num<64>(1) << 63) == 9223372036854775837_bn64);
		test(BigNum::nextPrimeOpt35711<64>(BigNum::Num<64>(1) << 63) == 9223372036854775837_bn64);

		test(BigNum::nextPrime<128>(BigNum::Num<64>(1) << 63) == 9223372036854775837_bn64);
		test(BigNum::nextPrimeOpt35711<128>(BigNum::Num<64>(1) << 63) == 9223372036854775837_bn64);

		test(BigNum::nextPrime<256>(BigNum::Num<64>(1) << 63) == 9223372036854775837_bn64);
		test(BigNum::nextPrimeOpt35711<256>(BigNum::Num<64>(1) << 63) == 9223372036854775837_bn64);

		test(BigNum::nextPrime<512>(BigNum::Num<64>(1) << 63) == 9223372036854775837_bn64);
		test(BigNum::nextPrimeOpt35711<512>(BigNum::Num<64>(1) << 63) == 9223372036854775837_bn64);

		test(BigNum::nextPrime<1024>(BigNum::Num<64>(1) << 63) == 9223372036854775837_bn64);
		test(BigNum::nextPrimeOpt35711<1024>(BigNum::Num<64>(1) << 63) == 9223372036854775837_bn64);

		test(BigNum::nextPrime<128>(BigNum::Num<128>(1) << 127) == 170141183460469231731687303715884105757_bn128);
		test(BigNum::nextPrimeOpt357<128>(BigNum::Num<128>(1) << 127) == 170141183460469231731687303715884105757_bn128);
		test(BigNum::nextPrimeOpt35711<128>(BigNum::Num<128>(1) << 127) == 170141183460469231731687303715884105757_bn128);

		test(BigNum::nextPrime<256>(BigNum::Num<256>(1) << 255) == 57896044618658097711785492504343953926634992332820282019728792003956564820063_bn256);
		test(BigNum::nextPrimeOpt357<256>(BigNum::Num<256>(1) << 255) == 57896044618658097711785492504343953926634992332820282019728792003956564820063_bn256);
		test(BigNum::nextPrimeOpt35711<256>(BigNum::Num<256>(1) << 255) == 57896044618658097711785492504343953926634992332820282019728792003956564820063_bn256);

		test(BigNum::nextPrime<512>(1_bn512 << 511) == 6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042159_bn512);
		test(BigNum::nextPrimeOpt3<512>(1_bn512 << 511) == 6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042159_bn512);
		test(BigNum::nextPrimeOpt35<512>(1_bn512 << 511) == 6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042159_bn512);
		test(BigNum::nextPrimeOpt357<512>(1_bn512 << 511) == 6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042159_bn512);
		test(BigNum::nextPrimeOpt35711<512>(1_bn512 << 511) == 6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042159_bn512);

		test(BigNum::nextPrime<1024>(1_bn1024 << 1023) == 89884656743115795386465259539451236680898848947115328636715040578866337902750481566354238661203768010560056939935696678829394884407208311246423715319737062188883946712432742638151109800623047059726541476042502884419075341171231440736956555270413618581675255342293149119973622969239858152417678164812112069763_bn1024);
		test(BigNum::nextPrimeOpt357<1024>(1_bn1024 << 1023) == 89884656743115795386465259539451236680898848947115328636715040578866337902750481566354238661203768010560056939935696678829394884407208311246423715319737062188883946712432742638151109800623047059726541476042502884419075341171231440736956555270413618581675255342293149119973622969239858152417678164812112069763_bn1024);
		test(BigNum::nextPrimeOpt35711<1024>(1_bn1024 << 1023) == 89884656743115795386465259539451236680898848947115328636715040578866337902750481566354238661203768010560056939935696678829394884407208311246423715319737062188883946712432742638151109800623047059726541476042502884419075341171231440736956555270413618581675255342293149119973622969239858152417678164812112069763_bn1024);

		test(BigNum::nextPrime<2048>(1_bn2048 << 2047) > 0); //== 16158503035655503650357438344334975980222051334857742016065172713762327569433945446598600705761456731844358980460949009747059779575245460547544076193224141560315438683650498045875098875194826053398028819192033784138396109321309878080919047169238085235290822926018152521443787945770532904303776199561965192760957166694834171210342487393282284747428088017663161029038902829665513096354230157075129296432088558362971801859230928678799175576150822952201848806616643615613562842355410104862578550863465661734839271290328348967522998634176499319107762583194718667771801067716614802322659239302476074096777926805529798117247_bn2048);
		test(BigNum::nextPrimeOpt35711<2048>(1_bn2048 << 2047) > 0);

		printf("biggestPrimeProduct<128>() = %s\n", toString(biggestPrimeProduct<128>()).c_str());
		printf("biggestPrimeProduct<256>() = %s\n", toString(biggestPrimeProduct<256>()).c_str());

		srand(0);
		{
			auto start = std::chrono::high_resolution_clock::now();
			BigNum::randPrime<512>(BigNum::RandType::Simple);
			BigNum::randPrime<512>(BigNum::RandType::Simple);
			auto stop = std::chrono::high_resolution_clock::now();
			auto elapsed = stop - start;
			printf("2x Primes (512) - Simple \n   duration - %lld ms\n", (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
		}

		srand(0);
		{
			auto start = std::chrono::high_resolution_clock::now();
			BigNum::randPrime<512>(BigNum::RandType::FindNext);
			BigNum::randPrime<512>(BigNum::RandType::FindNext);
			auto stop = std::chrono::high_resolution_clock::now();
			auto elapsed = stop - start;
			printf("2x Primes (512) - FindNext Simple\n   duration - %lld ms\n", (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
		}

		srand(0);
		{
			auto start = std::chrono::high_resolution_clock::now();
			BigNum::randPrime<512>(BigNum::RandType::FindNextOpt);
			BigNum::randPrime<512>(BigNum::RandType::FindNextOpt);
			auto stop = std::chrono::high_resolution_clock::now();
			auto elapsed = stop - start;
			printf("2x Primes (512) - FindNext Optimized\n   duration - %lld ms\n", (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
		}

		{
			unsigned char data[128] = { 0 };
			std::generate(&data[0], &data[128], []() { return (unsigned char)(rand() % 256u); });
			BigNum::Num<2048> packed = BigNum::d2i<2048>(&data[0], &data[128]);
			unsigned char unpacked[128] = { 0 };
			BigNum::i2d(packed, &unpacked[0], &unpacked[128]);
			test(memcmp(data, unpacked, 128) == 0);
		}

		{
			printf("Generate RSA keys.");
			constexpr int Module = 2048;

			auto start = std::chrono::high_resolution_clock::now();

			BigNum::Num<Module / 2> p = BigNum::randPrime<Module / 2>();
			BigNum::Num<Module / 2> q = BigNum::randPrime<Module / 2>();

			BigNum::Num<Module> N = p * q;
			BigNum::Num<Module> t = N - p - q + 1; // (p - 1)(q - 1) = pq - p - q + 1 = N - p - q + 1
			BigNum::Num<Module> e = BigNum::Num<Module>(65537);
			BigNum::Num<Module> d = BigNum::modInv(e, t);

			BigNum::Num<Module / 2> qinvp, dp, dq;
			BigNum::CRT_init(p, q, d, dp, dq, qinvp);

			auto stop = std::chrono::high_resolution_clock::now();
			auto elapsed = stop - start;
			printf(" duration - %lld ms\n", (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());

			unsigned char testdata[256];
			std::generate(std::begin(testdata), std::end(testdata), []() { return (unsigned char)(rand() % 256); });

			std::vector<unsigned char> decrypted;

			std::vector<BigNum::Num<Module>> cipperdata;

			constexpr unsigned int blockSize = (Module / 2) / 8;
			static_assert(blockSize % sizeof(BigNum::Word) == 0, "Module have to be power of 2");
			printf("Encrypt");
			start = std::chrono::high_resolution_clock::now();
			for (std::size_t block = 0; block < (sizeof(testdata)/sizeof(*testdata)); block += blockSize)
			{
				BigNum::Num<Module> data = BigNum::d2i<Module>(&testdata[block], &testdata[block + blockSize]);
				BigNum::Num<Module> cipper = BigNum::monModExp(data, e, N);

				cipperdata.push_back(cipper);
			}
			stop = std::chrono::high_resolution_clock::now();
			elapsed = stop - start;
			printf(" duration - %lld ms\n", (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());

			printf("Decrypt");
			start = std::chrono::high_resolution_clock::now();

			for (std::size_t block = 0; block < cipperdata.size(); block++)
			{
				BigNum::Num<Module> data = BigNum::CRT(cipperdata[block], p, q, dp, dq, qinvp);
				unsigned char decrypted_data[blockSize];
				BigNum::i2d(data, std::begin(decrypted_data), std::end(decrypted_data));
				decrypted.insert(decrypted.end(), &decrypted_data[0], &decrypted_data[blockSize]);
			}

			stop = std::chrono::high_resolution_clock::now();
			elapsed = stop - start;
			printf(" duration - %lld ms\n", (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());

			bool ret = decrypted == std::vector<unsigned char>(std::begin(testdata), std::end(testdata));
			printf("%s\n", (ret) ? "Ok" : "Fail");
			if (!ret)
			{
				for (auto i : testdata)
					printf("%02X ", i);

				printf("\n\n");

				for (auto i : decrypted)
					printf("%02X ", i);

				printf("\n");
			}
			if (!ret) throw(1);
		}
		printf("Passed %d tests\nFailed %d tests\n", TestsPass, TestsFailed);
	}
	catch (...)
	{

	}
#endif
	return 0;
}