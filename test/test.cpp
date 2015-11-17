#include "BigNum.h"
#include <assert.h>
#include <string>
#include <algorithm>

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

int main()
{
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
	assert((BigNum::Num<64>(135056 >> 4) == BigNum::shift_right(BigNum::Num<64>(135056), 4)));
	assert((BigNum::Num<64>(135056 >> 12) == BigNum::shift_right(BigNum::Num<64>(135056), 12)));
	assert((BigNum::Num<64>(13505675 >> 15) == BigNum::shift_right(BigNum::Num<64>(13505675), 15)));
	assert((BigNum::Num<64>(13505675 >> 22) == BigNum::shift_right(BigNum::Num<64>(13505675), 22)));

	assert(BigNum::Num<64>(524568123) / BigNum::Num<64>(45654) == BigNum::Num<64>(524568123 / 45654));
	assert(BigNum::Num<64>(524568123) % BigNum::Num<64>(45654) == BigNum::Num<64>(524568123 % 45654));

	printf("BigNum::Num<128>(1025445861) = %s\n", toString(BigNum::Num<128>(1025445861ull)).c_str());

	printf("BigNum::Num<128>(125632587586954854585354458654754812345678908765678765467543456765678) = %s\n", toString(fromString<1024>("125632587586954854585354458654754812345678908765678765467543456765678")).c_str());

	return 0;
}