#pragma once

#include <algorithm>

namespace BigNum
{
	namespace Core
	{
		template<typename It1, typename It2>
		int cmp(It1 v1begin, It1 v1end, It2 v2begin, It2 v2end) noexcept
		{
			const std::ptrdiff_t v1len = v1end - v1begin;
			const std::ptrdiff_t v2len = v2end - v2begin;
			const std::ptrdiff_t vmax = (v1len > v2len)? v1len : v2len;

			for(std::ptrdiff_t i = vmax; i-->0;)
			{
				const auto v1 = (i < v1len) ? v1begin[i] : 0;
				const auto v2 = (i < v2len) ? v2begin[i] : 0;
				if (v1 != v2)
					return v1 > v2 ? 1 : -1;
			}
			return 0;
		}

		bool adc(unsigned char& ret, const unsigned char v1, const unsigned char v2, bool carry) noexcept
		{
			unsigned short sum = v1 + v2 + (carry ? 1 : 0);
			ret = sum & 0xFF;
			return (sum >> 8) > 0;
		}

		bool sbc(unsigned char& ret, const unsigned char v1, const unsigned char v2, bool carry) noexcept
		{
			unsigned short sum = v1 - (v2 + (carry ? 1 : 0));
			ret = sum & 0xFF;
			return (sum >> 8) > 0;
		}

		template<typename ItRez, typename It1, typename It2>
		int add(ItRez rezbegin, ItRez rezend, It1 v1begin, It1 v1end, It2 v2begin, It2 v2end) noexcept
		{
			const std::ptrdiff_t rezlen = rezend - rezbegin;
			const std::ptrdiff_t v1len = v1end - v1begin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			//assert(rezlen >= v1len && rezlen >= v2len);

			static_assert(sizeof(*rezbegin) == sizeof(*v1begin), "type have to be the same");
			static_assert(sizeof(*rezbegin) == sizeof(*v2begin), "type have to be the same");

			static_assert(sizeof(int) > sizeof(*rezbegin), "cant detect overflow");

			bool carry = false;
			for (std::ptrdiff_t i = 0; i < rezlen; i++)
			{
				const auto v1 = (i < v1len) ? v1begin[i] : 0;
				const auto v2 = (i < v2len) ? v2begin[i] : 0;
				carry = adc(rezbegin[i], v1, v2, carry);
			}
			return carry;
		}

		template<typename ItRez, typename It1, typename It2>
		int sub(ItRez rezbegin, ItRez rezend, It1 v1begin, It1 v1end, It2 v2begin, It2 v2end) noexcept
		{
			const std::ptrdiff_t rezlen = rezend - rezbegin;
			const std::ptrdiff_t v1len = v1end - v1begin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			//assert(rezlen >= v1len && rezlen >= v2len);

			static_assert(sizeof(*rezbegin) == sizeof(*v1begin), "type have to be the same");
			static_assert(sizeof(*rezbegin) == sizeof(*v2begin), "type have to be the same");

			static_assert(sizeof(int) > sizeof(*rezbegin), "cant detect overflow");

			bool carry = false;
			for (std::ptrdiff_t i = 0; i < rezlen; i++)
			{
				const auto v1 = (i < v1len) ? v1begin[i] : 0;
				const auto v2 = (i < v2len) ? v2begin[i] : 0;
				carry = sbc(rezbegin[i], v1, v2, carry);
			}
			return carry;
		}
	}
	
	template<int N>
	class Num
	{
	public:
		typedef unsigned char Digit;

		enum
		{
			DigitSizeBits = sizeof(Digit) * 8,
			Size = (N + (DigitSizeBits) - 1) / (DigitSizeBits),
		};
	private:
		Digit data[Size];
	public:
		Num() noexcept
		{
			for (int i = 0; i < Size; i++)
				data[i] = 0;
		}

		template<typename T>
		Num(T num) noexcept
		{
			for (int i = 0; i < Size; i++)
			{
				data[i] = num & ((1 << DigitSizeBits) - 1);
				num >>= DigitSizeBits;
			}
		}

		template<int M>
		Num(const Num<M>& other) noexcept
		{
			*this = other;
		}

		Num& operator = (const Num& other) noexcept
		{
			if (this != &other)
			{
				for (int i = 0; i < Size; i++)
					data[i] = other[i];
			}
			return *this;
		}

		template<int M>
		Num& operator = (const Num<M>& other) noexcept
		{
			for (int i = 0; i < std::min<int>(Size, Num<M>::Size); i++)
				data[i] = other[i];

			for (int i = std::min<int>(Size, Num<M>::Size); i < Size; i++)
				data[i] = 0;
				
			assert(std::none_of(&other[std::min<int>(Size, Num<M>::Size)],&other[Num<M>::Size], [](unsigned char x){ return x != 0; }));

			return *this;
		}

		unsigned char& operator[](int i) noexcept { return data[i]; }
		const unsigned char& operator[](int i) const noexcept { return data[i]; }

		bool bit(int i) const noexcept
		{
			//assert(i >= 0 && i < Size);
			bool ret = false;
			if (i >= 0 && i < N)
			{
				const int byte = i / DigitSizeBits;
				const int bit = i % DigitSizeBits;
				ret = ((*this)[byte] & (1 << bit)) != 0;
			}
			return ret;
		}

		Digit* begin() noexcept { return &data[0]; }
		Digit* end() noexcept { return &data[Size]; }

		const Digit* begin() const noexcept { return &data[0]; }
		const Digit* end() const noexcept { return &data[Size]; }
	};

	template<int N, int M>
	const Num<N> operator + (const Num<N>& v1, const Num<M>& v2)
	{
		Num<N> rez;
		Core::add(rez.begin(), rez.end(), v1.begin(), v1.end(), v2.begin(), v2.end());
		return rez;
	}

	template<int N>
	const Num<N> operator + (Num<N> v1, unsigned char v2)
	{
		Core::add(v1.begin(), v1.end(), v1.begin(), v1.end(), &v2, (&v2)+sizeof(v2));
		return v1;
	}

	template<int N>
	const Num<N> operator - (const Num<N>& v)
	{
		return Num<N>(0) - v;
	}

	template<int N>
	const Num<N> operator - (const Num<N>& v1, const Num<N>& v2)
	{
		Num<N> rez;
		Core::sub(rez.begin(), rez.end(), v1.begin(), v1.end(), v2.begin(), v2.end());
		return rez;
	}

	template<int N>
	const Num<N> operator - (Num<N> v1, unsigned char v2)
	{
		Core::sub(v1.begin(), v1.end(), v1.begin(), v1.end(), &v2, (&v2) + sizeof(v2));
		return v1;
	}

	template<int N>
	const Num<2*N> mul2N (const Num<N>& v1, const Num<N>& v2)
	{
		unsigned int poly[Num<N>::Size * 2] = { 0 };

		for (int i = 0; i < Num<N>::Size; i++)
			for (int j = 0; j < Num<N>::Size; j++)
				poly[i + j] += v1[i] * v2[j];

		Num<2*N> rez;
		unsigned int carry = 0;
		for (int i = 0; i < Num<2*N>::Size; i++)
		{
			int sum = carry + poly[i];
			rez[i] = static_cast<unsigned char>(sum & 0xFF);
			carry = sum >> 8;
		}
		return rez;
	}

	template<int N>
	const Num<N> operator * (const Num<N>& v1, const Num<N>& v2)
	{
		return mul2N(v1, v2);
	}

	template<int N>
	const Num<N> operator * (const Num<N>& v1, unsigned char v2)
	{
		Num<N> rez;
		unsigned int carry = 0;
		for (int i = 0; i < Num<N>::Size; i++)
		{
			unsigned int sum = carry + v1[i] * v2;
			rez[i] = sum & 0xFF;
			carry = sum >> 8;
		}
		return rez;
	}

	template<int N, int M>
	const bool operator == (const Num<N>& v1, const Num<M>& v2)
	{
		return Core::cmp(&v1[0], &v1[v1.Size], &v2[0], &v2[v2.Size]) == 0;
	}

	template<int N>
	const bool operator == (const Num<N>& v1, unsigned char v2)
	{
		return Core::cmp(&v1[0], &v1[v1.Size], &v2, (&v2) + sizeof(v2)) == 0;
	}

	template<int N, int M>
	const bool operator != (const Num<N>& v1, const Num<M>& v2)
	{
		return !(v1 == v2);
	}
	
	template<int N>
	const bool operator != (const Num<N>& v1, unsigned char v2)
	{
		return !(v1 == v2);
	}

	template<int N, int M>
	const bool operator > (const Num<N>& v1, const Num<M>& v2)
	{
		return Core::cmp(&v1[0], &v1[v1.Size], &v2[0], &v2[v2.Size]) > 0;
	}

	template<int N, int M>
	const bool operator >= (const Num<N>& v1, const Num<M>& v2)
	{
		return Core::cmp(&v1[0], &v1[v1.Size], &v2[0], &v2[v2.Size]) >= 0;
	}

	template<int N, int M>
	const bool operator < (const Num<N>& v1, const Num<M>& v2)
	{
		return Core::cmp(&v1[0], &v1[v1.Size], &v2[0], &v2[v2.Size]) < 0;
	}

	template<int N, int M>
	const bool operator <= (const Num<N>& v1, const Num<M>& v2)
	{
		return Core::cmp(&v1[0], &v1[v1.Size], &v2[0], &v2[v2.Size]) <= 0;
	}

	template<int N>
	const Num<N> shift_left_bytes(const Num<N>& v, int count)
	{
		Num<N> ret;
		for (int i = Num<N>::Size - 1; i >= int(count); i--)
		{
			ret[i] = v[i - count];
		}
		return ret;
	}

	template<int N>
	const Num<N> shift_left(const Num<N>& v, int count)
	{
		Num<N> ret;
		if (count == 0)
		{
			ret = v;
		}
		else if (count > 0 && count < N)
		{
			const int bytes = count / Num<N>::DigitSizeBits;
			const int bits = count % Num<N>::DigitSizeBits;
			if (bits == 0)
			{
				ret = shift_left_bytes(v, bytes);
			}
			else
			{
				unsigned int carry = 0;
				for (int i = bytes; i < Num<N>::Size; i++)
				{
					unsigned int val = static_cast<unsigned int>(v[i - bytes]) << bits;
					ret[i] = static_cast<Num<N>::Digit>(val | carry);
					carry = val >> Num<N>::DigitSizeBits;
				}
			}
		}
		return ret;
	}

	template<int N>
	const Num<N> operator << (const Num<N>& v1, int bits)
	{
		return shift_left(v1, bits);
	}

	template<int N>
	const Num<N> shift_right_bytes(const Num<N>& v, int count)
	{
		Num<N> ret;
		for (int i = 0; i < Num<N>::Size - count; i++)
		{
			ret[i] = v[i + count];
		}
		return ret;
	}

	template<int N>
	const Num<N> shift_right(const Num<N>& v, int count)
	{
		Num<N> ret;
		if (count == 0)
		{
			ret = v;
		}
		else if (count > 0 && count < N)
		{
			const int bytes = count / Num<N>::DigitSizeBits;
			const int bits = count % Num<N>::DigitSizeBits;
			if (bits == 0)
			{
				ret = shift_right_bytes(v, bytes);
			}
			else
			{
				unsigned int carry = v[bytes];
				for (int i = 0; i < Num<N>::Size - bytes; i++)
				{
					unsigned int val = (static_cast<unsigned int>(((i + bytes + 1) < Num<N>::Size) ? v[i + bytes + 1] : 0) << Num<N>::DigitSizeBits);
					val = (val | carry);
					ret[i] = static_cast<Num<N>::Digit>(val >> bits);
					carry = val >> (Num<N>::DigitSizeBits);
				}
			}
		}
		return ret;
	}

	template<int N>
	const Num<N> operator >> (const Num<N>& v1, int bits)
	{
		return shift_right(v1, bits);
	}

	template<int N, int M>
	void div(const Num<N>& n, const Num<M>& d, Num<N>& q, Num<M>& r)
	{
		Num<M + Num<M>::DigitSizeBits> rest;

		const Num<M + Num<M>::DigitSizeBits> shiftedD[Num<N>::DigitSizeBits] = 
			{ 
				Num<M + 8>(d), 
				Num<M + 8>(d) << 1, 
				Num<M + 8>(d) << 2, 
				Num<M + 8>(d) << 3, 
				Num<M + 8>(d) << 4, 
				Num<M + 8>(d) << 5, 
				Num<M + 8>(d) << 6, 
				Num<M + 8>(d) << 7 
			};

		for (int i = Num<N>::Size - 1; i >= 0; i--)
		{
			rest = shift_left_bytes(rest, 1u);
			rest[0] = n[i];
			unsigned int val = 0;
			if (rest >= shiftedD[0])
			{
				for (int j = Num<M>::DigitSizeBits - 1; j >= 0; --j)
				{
					if (rest >= shiftedD[j])
					{
						rest = rest - shiftedD[j];
						val += 1 << j;
					}
				}
			}
			q[i] = static_cast<unsigned char>(val);
		}
		r = rest;
	}

	//template<int N>
	//void div(const Num<N>& n, unsigned char d, Num<N>& q, unsigned char& r)
	//{
	//	unsigned short rest = 0;
	//	for (int i = Num<N>::Size - 1; i >= 0; i--)
	//	{
	//		rest = rest<<8;
	//		rest+= n[i];
	//		q[i] = (unsigned char)(rest / d);
	//		rest = rest % d;
	//	}
	//	r = (unsigned char)rest;
	//}

	//template<int N>
	//void div(const Num<N>& n, unsigned short d, Num<N>& q, unsigned short& r)
	//{
	//	unsigned int rest = 0;
	//	for (int i = Num<N>::Size - 1; i >= 0; i--)
	//	{
	//		rest = rest << 8;
	//		rest += n[i];
	//		q[i] = (unsigned char)(rest / d);
	//		rest = rest % d;
	//	}
	//	r = (unsigned short)rest;
	//}

	template<int N>
	void div(const Num<N>& n, unsigned int d, Num<N>& q, unsigned int& r)
	{
		unsigned long long int rest = 0;
		for (int i = Num<N>::Size - 1; i >= 0; i--)
		{
			rest = rest << Num<N>::DigitSizeBits;
			rest += n[i];
			q[i] = (unsigned char)(rest / d);
			rest = rest % d;
		}
		r = (unsigned int)rest;
	}

	template<int N, int M>
	const Num<N> operator / (const Num<N>& v1, const Num<M>& v2)
	{
		if (v1 < v2)
			return Num<N>(0);

		Num<N> q;
		Num<M> r;
		div(v1, v2, q, r);
		return q;
	}

	//template<int N>
	//const Num<N> operator / (const Num<N>& v1, const unsigned char v2)
	//{
	//	Num<N> q;
	//	unsigned char r = 0;
	//	div(v1, v2, q, r);
	//	return q;
	//}

	//template<int N>
	//const Num<N> operator / (const Num<N>& v1, const unsigned short v2)
	//{
	//	Num<N> q;
	//	unsigned short r = 0;
	//	div(v1, v2, q, r);
	//	return q;
	//}

	template<int N>
	const Num<N> operator / (const Num<N>& v1, const unsigned int v2)
	{
		Num<N> q;
		unsigned int r = 0;
		div(v1, v2, q, r);
		return q;
	}

	template<int N, int M>
	const Num<M> operator % (const Num<N>& v1, const Num<M>& v2)
	{
		if (v1 < v2)
			return v1;

		Num<N> q;
		Num<M> r;
		div(v1, v2, q, r);
		return r;
	}

	//template<int N>
	//const unsigned char operator % (const Num<N>& v1, const unsigned char v2)
	//{
	//	Num<N> q;
	//	unsigned char r = 0;
	//	div(v1, v2, q, r);
	//	return r;
	//}

	//template<int N>
	//const unsigned short operator % (const Num<N>& v1, const unsigned short v2)
	//{
	//	Num<N> q;
	//	unsigned short r = 0;
	//	div(v1, v2, q, r);
	//	return r;
	//}

	template<int N>
	const unsigned int operator % (const Num<N>& v1, const unsigned int v2)
	{
		Num<N> q;
		unsigned int r = 0;
		div(v1, v2, q, r);
		return r;
	}

	template<int N>
	const Num<N> operator & (const Num<N>& v1, const Num<N>& v2)
	{
		Num<N> rez;

		for (int i = 0; i < Num<N>::Size; i++)
		{
			rez[i] = v1[i] & v2[i];
		}

		return rez;
	}

	template<int N>
	const Num<N> operator & (const Num<N>& v1, Num<N>&& v2)
	{
		Num<N> rez;

		for (int i = 0; i < Num<N>::Size; i++)
		{
			rez[i] = v1[i] & v2[i];
		}

		return rez;
	}

	template<int N>
	const Num<N> operator & (const Num<N>& v1, const unsigned char v2)
	{
		Num<N> rez(v1[0] & v2);
		return rez;
	}

	template<int N>
	const Num<N> operator | (const Num<N>& v1, const Num<N>& v2)
	{
		Num<N> rez;

		for (int i = 0; i < Num<N>::Size; i++)
		{
			rez[i] = v1[i] | v2[i];
		}

		return rez;
	}

	template<int N>
	const Num<N> operator | (const Num<N>& v1, const unsigned char v2)
	{
		Num<N> rez = v1;
		rez[0] = v1[0] | v2;
		return rez;
	}

	template<int N>
	const Num<N> operator ~ (const Num<N>& v)
	{
		Num<N> rez;

		for (int i = 0; i < Num<N>::Size; i++)
		{
			rez[i] = ~v[i];
		}

		return rez;
	}

	template<int N, int M>
	bool operator && (const Num<N>& v1, const Num<M>& v2)
	{
		return v1 != 0 && v2 != 0;
	}

	template<int N, int M>
	bool operator || (const Num<N>& v1, const Num<M>& v2)
	{
		return v1 != 0 || v2 != 0;
	}

	template<int N>
	bool operator !(const Num<N>& v)
	{
		return v == 0;
	}
	
	template<int N>
	const Num<N> gcd(const Num<N>& v1, const Num<N>& v2)
	{
		if (v1 == v2)
			return v1;
			
		Num<N> a = (v1 > v2) ? v1 : v2;
		Num<N> b = (v1 > v2) ? v2 : v1; 
		
		while(b != 0)
		{
			Num<N> t = b;
			b = a % b;
			a = t;
		}
		return a;
	}

	template<int N>
	const Num<2*N> lcm(const Num<N>& v1, const Num<N>& v2)
	{
		Num<2 * N> m = mul2N(v1, v2);
		Num<N> g = gcd(v1, v2);
		return m / g;
	}

	template<int N>
	const Num<N> modInv(const Num<N>& a, const Num<N>& n)
	{
	    Num<N> t(0), newt(1);    
	    Num<N> r = n, newr = a;    
	    while (newr != 0)
	    {
	        Num<N> quotient = r / newr;
	        {
	        	//(t, newt) := (newt, t - quotient * newt) 
	        	Num<N> tmp = newt;
	        	newt = t - quotient * newt;
	        	t = tmp;
	        }
	        {
	        	//(r, newr) := (newr, r - quotient * newr)
	        	Num<N> tmp = newr;
	        	newr = r - quotient * newr;
	        	r = tmp;
	        }
	    }
	    if (r > Num<N>(1)) return 0;
		if (t > n)
		{
			// signed < 0
			t = t + n;
		}
	    return t;
	}

	template<typename T>
	const T modInv(const T& a, const T& n)
	{
		T t(0), newt(1);
		T r = n, newr = a;
		while (newr != T(0))
		{
			T quotient = r / newr;
			{
				//(t, newt) := (newt, t - quotient * newt) 
				T tmp = newt;
				newt = t - quotient * newt;
				t = tmp;
			}
			{
				//(r, newr) := (newr, r - quotient * newr)
				T tmp = newr;
				newr = r - quotient * newr;
				r = tmp;
			}
		}
		if (r > T(1)) return 0;
		if (t < 0) t = t + n;
		return t;
	}

	template<int N>
	const Num<N> modExp(const Num<N>& base, const Num<N>& exp, const Num<N>& mod)
	{
		if (mod == 1)
			return Num<N>(0);

		Num<N> result = 1;
		Num<N> power = base;

		int realExpSize = Num<N>::Size;
		// skip leading zeros
		while (realExpSize >= 0 && exp[realExpSize - 1] == 0)
			realExpSize--;

		for (int i = realExpSize * 8; i-->0;)
		{
			{
				Num<2 * N> q;
				div(mul2N(result, result), mod, q, result);
			}
			if (exp.bit(i))
			{
				Num<2 * N> q;
				div(mul2N(result, power), mod, q, result);
			}
		}

		return Num<N>(result);
	}

	template<int N>
	struct MonMul
	{
		const Num<N> n;
		const Num<2 * N> r;
		Num<N> rinv;
		Num<N> ninv;

		MonMul(const Num<N>& _n) : n(_n), r(Num<2 * N>(1) << N)
		{
			assert(gcd(Num<2 * N>(n), r) == 1);

			rinv = Num<N>(modInv<2 * N>(r, n));

			assert((Num<2 * N>(rinv) << N) % Num<2 * N>(n) == 1);

			Num<N> tempr;
			Num<2 * N> q;
			div((Num<2 * N>(rinv)<<N) - 1, n, q, tempr);
			ninv = Num<N>(q);

			assert((Num<2 * N>(rinv) << N) - mul2N(ninv,n) == 1);
		}

		const Num<N> In(const Num<N>& x) const
		{
			Num<2 * N> ret(x);
			ret = ret << N;

			Num<N> tempr;
			Num<2 * N> q;
			div(ret, n, q, tempr);
			return tempr;
		}

		template<int M>
		const Num<N> In(const Num<M>& x) const
		{
			Num<N+M> ret(x);
			ret = ret << N;

			Num<N> tempr;
			Num<N+M> q;
			div(ret, n, q, tempr);
			return tempr;
		}

		const Num<N> Out(const Num<N>& x) const
		{
			return (*this)(x, 1);
		}

		const Num<N> operator()(const Num<N>& a, const Num<N>& b) const
		{
			Num<2*N> t = mul2N(a, b);
			Num<N> m = mul2N(Num<N>(t & (r - 1)), ninv) & (r - 1);
			Num<N> u = (t + mul2N(m, n)) >> N;
			while (u >= n)
				u = u - n;
			return Num<N>(u);
		}
	};

	template<int N>
	const Num<N> monModMul(const Num<N>& a, const Num<N>& b, const Num<N>& mod)
	{
		MonMul<N> monMod((mod));
		Num<N> xa = monMod.In(a);
		return monMod(xa, b);
	}

	template<int N, int M, int K>
	const Num<N> monModExp(const Num<K>& base, const Num<M>& exp, const Num<N>& mod)
	{
		MonMul<N> monMod((mod));
		Num<N> ret = monMod.In(Num<N>(1));
		Num<N> power = monMod.In(base);

		int realExpSize = Num<M>::Size;
		while (realExpSize >= 0 && exp[realExpSize - 1] == 0)
			realExpSize--;

		for (int i = realExpSize * Num<M>::DigitSizeBits; i -->0 ;)
		{
			ret = monMod(ret, ret);
			if (exp.bit(i))
				ret = monMod(ret, power);
		}
		return monMod.Out(ret);
	}

	template<int N>
	const Num<N> rand()
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
		BigNum::Num<N> a_to_power = BigNum::monModExp<N>(a, d, num);

		if (a_to_power == 1)
			return true;

		for (int i = 0; i < s - 1; i++)
		{
			if (a_to_power == num - 1)
				return true;

			BigNum::Num<N * 2> temp = BigNum::mul2N(a_to_power, a_to_power);

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
		constexpr int times = 3;
		for (int i = 0; i < times; i++)
		{
			BigNum::Num<N> a = (rand<N>() % (num - 2)) + 2;
			if (millerRabinPass(num, a) == false)
				return false;
		}
		return true;
	}

	constexpr unsigned short primes[] = { 3, 5, 7, 11, 13,  17, 19, 23, 29,
		31,     37,     41,     43,     47,     53,     59,     61,     67,     71,
		73,     79,     83,     89,     97,     101,    103,    107,    109,    113,
		127,    131,    137,    139,    149,    151,    157,    163,    167,    173,
		179,    181,    191,    193,    197,    199,    211,    223,	227,    229,
		233,    239,    241,    251,    257,    263,    269,    271,    277,    281,
		283,    293,    307,    311,    313,    317,    331,    337,    347,    349,
		353,    359,    367,    373,    379,    383,    389,    397,    401,    409,
		419,    421,    431,    433,    439,    443,    449,    457,    461,    463,
		467,    479,    487,    491,    499,    503,    509,    521,    523,    541,
		547,    557,    563,    569,    571,    577,    587,    593,    599,    601,
		607,    613,    617,    619,    631,    641,    643,    647,    653,    659,
		661,    673,    677,    683,    691,    701,    709,    719,    727,    733,
		739,    743,    751,    757,    761,    769,    773,    787,    797,    809,
		811,    821,    823,    827,    829,    839,    853,    857,    859,    863,
		877,    881,    883,    887,    907,    911,    919,    929,    937,    941,
		947,    953,    967,    971,    977,    983,    991,    997, };

	constexpr int simplechecks(int size)
	{
		return (size <= 256) ? 30 : ((size <= 512) ? 50 : ((size <= 1024) ? 80 : ((size <= 2048) ? 128 : (sizeof(primes) / sizeof(primes[0])))));
	}

	template<int N>
	const Num<N> findPrime(const BigNum::Num<N>& from)
	{
		BigNum::Num<N> prime = from | 1;

		unsigned short rests[simplechecks(N)] = { 0 };

		bool simpleTest = true;
		for (int i = 0; i < sizeof(rests) / sizeof(rests[0]); i++)
		{
			rests[i] = prime % primes[i];
			simpleTest &= (rests[i] != 0);
		}

		while (!(simpleTest && millerRabinTest(prime)))
		{
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
			prime = prime + 2;
		}

		return prime;
	}

	template<int N>
	const Num<N> randPrime()
	{
		BigNum::Num<N> num = rand<N>() | (BigNum::Num<N>(1) | (BigNum::Num<N>(1) << (N - 1)));
		return findPrime(num);
	}
};