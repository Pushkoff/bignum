#include <assert.h>
#include <string>
#include <string.h>
#include <algorithm>
#include <chrono>
#include <vector>
#include <time.h>
#include "BigNum.h"

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

#define test(a) do{ bool ret = (a); printf("%s\n - %s\n\n", #a, (ret) ? "Ok" : "Fail" ); assert(ret); }while(false);

#define PROFILING 0

int main()
{
#if PROFILING
	auto prime = findPrime<2048>(BigNum::Num<2048>(1) << 2047);
#else
	//srand(time(nullptr));
	
	test(BigNum::Num<256>(fromString<128>("118802731")) == fromString<256>("118802731"));
	test(fromString<128>(toString(fromString<128>("118802731")).c_str()) == fromString<128>("118802731"));

	printf("%s\n", toString(fromString<128>("118802731")).c_str());

	test(BigNum::Num<64>(0) == BigNum::Num<64>(0));
	test(BigNum::Num<64>(0) != BigNum::Num<64>(1));

	test(BigNum::Num<64>(1) + BigNum::Num<64>(1) == BigNum::Num<64>(2));
	test(BigNum::Num<64>(2) - BigNum::Num<64>(1) == BigNum::Num<64>(1));

	test(BigNum::Num<64>(256) - BigNum::Num<64>(1) == BigNum::Num<64>(255));
	test(BigNum::Num<16>(255<<8) - BigNum::Num<16>(1) == BigNum::Num<16>((255 << 8) -1));

	test(BigNum::Num<64>(256) - 1 == BigNum::Num<64>(255));
	test(BigNum::Num<16>(255 << 8) - 1 == BigNum::Num<16>((255 << 8) - 1));

	test(BigNum::Num<64>(20) * BigNum::Num<64>(45) == BigNum::Num<64>(20 * 45));

	BigNum::Num<1024> fact(1);
	for (unsigned char i = 2; i <= 100; i++)
	{
		fact = fact * i;
	}
	test(fact == fromString<1024>("93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000"));
	printf("100!=%s\n", toString(fact).c_str());


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

	test((BigNum::Num<64>(1 << 8) == BigNum::shift_left(BigNum::Num<64>(1), 8)));
	test((BigNum::Num<64>(1 << 7) == BigNum::shift_left(BigNum::Num<64>(1), 7)));
	test((BigNum::Num<64>(135 << 12) == BigNum::shift_left(BigNum::Num<64>(135), 12)));
	test((BigNum::Num<64>(135056 << 12) == BigNum::shift_left(BigNum::Num<64>(135056), 12)));
	test((BigNum::Num<64>(135056 << 7) == BigNum::shift_left(BigNum::Num<64>(135056), 7)));
	test((BigNum::Num<64>(13056 << 15) == BigNum::shift_left(BigNum::Num<64>(13056), 15)));

	test((BigNum::Num<64>(1) == BigNum::shift_right(BigNum::Num<64>(1<<8), 8)));
	test((BigNum::Num<64>(135056 >> 1) == BigNum::shift_right(BigNum::Num<64>(135056), 1)));
	test((BigNum::Num<64>(135056 >> 4) == BigNum::shift_right(BigNum::Num<64>(135056), 4)));
	test((BigNum::Num<64>(135056 >> 12) == BigNum::shift_right(BigNum::Num<64>(135056), 12)));
	test((BigNum::Num<64>(13505675 >> 15) == BigNum::shift_right(BigNum::Num<64>(13505675), 15)));
	test((BigNum::Num<64>(13505675 >> 22) == BigNum::shift_right(BigNum::Num<64>(13505675), 22)));
	test(BigNum::Num<1024>(1) << 1023 == fromString<1024>("89884656743115795386465259539451236680898848947115328636715040578866337902750481566354238661203768010560056939935696678829394884407208311246423715319737062188883946712432742638151109800623047059726541476042502884419075341171231440736956555270413618581675255342293149119973622969239858152417678164812112068608"));

	{
		bool passed = true;
		BigNum::Num<64> bn(1);
		for (int i = 0; i < 64; i++)
		{
			BigNum::Num<64> testbn = bn << i;
			passed = passed && (testbn[i/8] == (1 << i%8));
		}
		test(passed == true);
	}
	
	{
		bool passed = true;
		BigNum::Num<64> bn = BigNum::Num<64>(1) << 63;
		for (int i = 0; i < 64; i++)
		{
			BigNum::Num<64> testbn = bn >> i;
			passed = passed && (testbn[(63 - i)/8] == (1 << (7 - i%8)));
		}
		test(passed == true && "Rignt shift");
	}
	
	{
		bool passed = true;
		BigNum::Num<64> bn1 = BigNum::Num<64>(255);
		BigNum::Num<64> bn2 = BigNum::Num<64>(255) << 56;
		for (int i = 0; i < 56; i++)
		{
			BigNum::Num<64> testbn1 = bn1 << i;
			BigNum::Num<64> testbn2 = bn2 >> (56-i);
			passed = passed && (testbn1 == testbn2);
		}
		test(passed == true && "shifts");
	}

	test(((BigNum::Num<64>(1) << 63) + 1) / (BigNum::Num<64>(1) << 63) == 1);
	test(((BigNum::Num<64>(1) << 63) + 1) % (BigNum::Num<64>(1) << 63) == 1);
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
	test((BigNum::modInv(BigNum::Num<128>(961748941), BigNum::Num<128>(982451653)) * BigNum::Num<128>(961748941)) % BigNum::Num<128>(982451653) == BigNum::Num<128>(1));

	//BigNum::MonMul<16> mod11(BigNum::Num<16>(121));
	//test(mod11(BigNum::Num<16>(100), BigNum::Num<16>(100)) == BigNum::Num<16>(78));
	test(BigNum::monModMul<32>(BigNum::Num<32>(43), BigNum::Num<32>(56), BigNum::Num<32>(97)) == BigNum::Num<32>(80));
	test(BigNum::monModExp<32>(BigNum::Num<32>(4), BigNum::Num<32>(13), BigNum::Num<32>(497)) == BigNum::Num<32>(445));
	test(BigNum::monModExp<32>(BigNum::Num<32>(4), BigNum::Num<32>(0), BigNum::Num<32>(497)) == BigNum::Num<32>(1));
	test(BigNum::monModExp<32>(BigNum::Num<32>(4), BigNum::Num<32>(1), BigNum::Num<32>(497)) == BigNum::Num<32>(4));
	test(BigNum::monModExp<64>(fromString<64>("8814054284918744181"), fromString<64>("1152921504606846977"), fromString<64>("9223372036854775817")) == fromString<64>("4724919961687496308"));

	//millerRabinPass(fromString<128>("118802731"), fromString<128>("74812"));

	srand(0);
	
	test(millerRabinTest(fromString<32>("17")) == true);

	test(millerRabinTest(fromString<128>("961748941")) == true);
	test(millerRabinTest(fromString<128>("982451653")) == true);
	test(millerRabinTest(fromString<256>("89685068870671644428994813094690803553")) == true);
	test(millerRabinTest(fromString<128>("982451655")) == false);

	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = BigNum::findPrime<64>(BigNum::Num<64>(1) << 63);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime (64)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = BigNum::findPrime<128>(BigNum::Num<128>(1) << 63);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime (128)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = BigNum::findPrime<256>(BigNum::Num<256>(1) << 63);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime (256)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = BigNum::findPrime<512>(BigNum::Num<512>(1) << 63);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime (512)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
		{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = BigNum::findPrime<1024>(BigNum::Num<1024>(1) << 63);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime (1024)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = BigNum::findPrime<128>(BigNum::Num<128>(1) << 127);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime (128)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = BigNum::findPrime<256>(BigNum::Num<256>(1) << 255);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime (256)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = BigNum::findPrime<512>(BigNum::Num<512>(1) << 511);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime (512)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		auto prime = BigNum::findPrime<1024>(BigNum::Num<1024>(1) << 1023);
		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf("Prime (1024)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	}
	//{
	//	auto start = std::chrono::high_resolution_clock::now();
	//	auto prime = findPrime<2048>(BigNum::Num<2048>(1) << 2047);
	//	auto stop = std::chrono::high_resolution_clock::now();
	//	auto elapsed = stop - start;
	//	printf("Prime (2048)= %s\n   duration - %lld ms\n", toString(prime).c_str(), (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
	//}

	{
		
		printf("Generate RSA keys.");
		auto start = std::chrono::high_resolution_clock::now();

		auto p = BigNum::findPrime<512>(BigNum::Num<512>(1) << 511);
		auto q = BigNum::findPrime<512>((BigNum::Num<512>(1) << 511) + (BigNum::Num<512>(1) << 32));

		auto N = BigNum::mul2N(p, q);
		auto t = BigNum::mul2N(p - 1, q - 1);
		auto e = BigNum::Num<1024>(65537);
		auto d = BigNum::modInv(e, t);

		auto qinvp = BigNum::modInv(q, p);
		auto dp = d % (p - 1);
		auto dq = d % (q - 1);

		auto stop = std::chrono::high_resolution_clock::now();
		auto elapsed = stop - start;
		printf(" duration - %lld ms\n", (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());

		std::vector<char> testdata(10240);
		std::generate(std::begin(testdata), std::end(testdata), []() {return rand() % 256; });

		std::vector<char> decrypted;

		std::vector<char> cipperdata;

		printf("Encrypt");
		start = std::chrono::high_resolution_clock::now();
		for (std::size_t block = 0; block < (testdata.size() / (BigNum::Num<1024>::Size / 2)); block++)
		{
			BigNum::Num<1024> data(0);
			for (int i = 0; i < BigNum::Num<1024>::Size/2; ++i)
				data[i] = testdata[block * BigNum::Num<1024>::Size / 2 + i];

			BigNum::Num<1024> cipper = BigNum::monModExp(data, e, N);

			for (int i = 0; i < BigNum::Num<1024>::Size; ++i)
				cipperdata.push_back(cipper[i]);
		}
		stop = std::chrono::high_resolution_clock::now();
		elapsed = stop - start;
		printf(" duration - %lld ms\n", (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());	

		printf("Decrypt");
		start = std::chrono::high_resolution_clock::now();

		for (std::size_t block = 0; block < (cipperdata.size() / BigNum::Num<1024>::Size); block++)
		{
			BigNum::Num<1024> cipper(0);
			for (int i = 0; i < BigNum::Num<1024>::Size; ++i)
				cipper[i] = cipperdata[block * BigNum::Num<1024>::Size + i];

			//BigNum::Num<1024> data = BigNum::monModExp(cipper, d, N);

			BigNum::Num<512> datap = BigNum::monModExp(cipper, dp, p);
			BigNum::Num<512> dataq = BigNum::monModExp(cipper, dq, q);

			BigNum::Num<1024> data = BigNum::Num<1024>(dataq) + BigNum::mul2N(BigNum::mul2N(datap - dataq, qinvp) % p, q);
			for (int i = 0; i < BigNum::Num<1024>::Size/2; ++i)
				decrypted.push_back(data[i]);
		}

		stop = std::chrono::high_resolution_clock::now();
		elapsed = stop - start;
		printf(" duration - %lld ms\n", (long long)std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());

		printf("%s\n", (decrypted == testdata) ? "Ok" : "Fail");
	}
	printf("Press any key...");
	int ret = getchar();
	(void)ret;
#endif
	return 0;
}