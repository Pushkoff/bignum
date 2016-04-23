#pragma once

#include <algorithm>
#include <array>

namespace BigNum
{
	typedef unsigned short Digit;
	constexpr size_t DigitSizeBits = sizeof(Digit) * 8;
	constexpr unsigned int DigitMaskBits = ((1ull << DigitSizeBits) - 1ull);

	namespace Core
	{
		int cmp(const Digit* v1begin, const Digit* v1end, const Digit* v2begin, const Digit* v2end) noexcept
		{
			const std::ptrdiff_t v1len = v1end - v1begin;
			const std::ptrdiff_t v2len = v2end - v2begin;
			const std::ptrdiff_t vmax = (v1len > v2len)? v1len : v2len;

			for(std::ptrdiff_t i = vmax; i-->0;)
			{
				const Digit v1 = (i < v1len) ? v1begin[i] : 0;
				const Digit v2 = (i < v2len) ? v2begin[i] : 0;
				if (v1 != v2)
					return v1 > v2 ? 1 : -1;
			}
			return 0;
		}

		Digit adc(Digit& ret, const Digit v1, const Digit v2, Digit carry) noexcept
		{
			static_assert(sizeof(long long) > sizeof(Digit), "cant detect overflow");

			long long sum = v1 + v2 + carry;
			ret = sum & DigitMaskBits;
			return Digit(sum >> DigitSizeBits);
		}

		Digit add(Digit* rezbegin, Digit* rezend, const Digit* v1begin, const Digit* v1end, const Digit* v2begin, const Digit* v2end) noexcept
		{
			const std::ptrdiff_t rezlen = rezend - rezbegin;
			const std::ptrdiff_t v1len = v1end - v1begin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			Digit carry = 0;
			for (std::ptrdiff_t i = 0; i < rezlen; i++)
			{
				const Digit v1 = (i < v1len) ? v1begin[i] : 0;
				const Digit v2 = (i < v2len) ? v2begin[i] : 0;
				carry = adc(rezbegin[i], v1, v2, carry);
			}
			return carry;
		}

		Digit add(Digit* rezbegin, Digit* rezend, const Digit* v2begin, const Digit* v2end) noexcept
		{
			const std::ptrdiff_t rezlen = rezend - rezbegin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			Digit carry = 0;
			for (std::ptrdiff_t i = 0; i < rezlen && (i < v2len || carry != 0); i++)
			{
				const Digit v2 = (i < v2len) ? v2begin[i] : 0;
				carry = adc(rezbegin[i], rezbegin[i], v2, carry);
			}
			return carry;
		}

		Digit sbc(Digit& ret, const Digit v1, const Digit v2, Digit carry) noexcept
		{
			static_assert(sizeof(long long) > sizeof(Digit), "cant detect overflow");

			long long sum = v1 - (v2 + carry);
			ret = sum & DigitMaskBits;
			return Digit(-(sum >> DigitSizeBits));
		}

		Digit sub(Digit* rezbegin, Digit* rezend, const Digit* v1begin, const Digit* v1end, const Digit* v2begin, const Digit* v2end) noexcept
		{
			const std::ptrdiff_t rezlen = rezend - rezbegin;
			const std::ptrdiff_t v1len = v1end - v1begin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			Digit carry = 0;
			for (std::ptrdiff_t i = 0; i < rezlen; i++)
			{
				const Digit v1 = (i < v1len) ? v1begin[i] : 0;
				const Digit v2 = (i < v2len) ? v2begin[i] : 0;
				carry = sbc(rezbegin[i], v1, v2, carry);
			}
			return carry;
		}

		Digit sub(Digit* rezbegin, Digit* rezend, const Digit* v2begin, const Digit* v2end) noexcept
		{
			const std::ptrdiff_t rezlen = rezend - rezbegin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			Digit carry = 0;
			for (std::ptrdiff_t i = 0; i < rezlen && (i < v2len || carry != 0); i++)
			{
				const Digit v2 = (i < v2len) ? v2begin[i] : 0;
				carry = sbc(rezbegin[i], rezbegin[i], v2, carry);
			}
			return carry;
		}

		void mul(Digit* rezbegin, Digit* rezend, const Digit* v1begin, const Digit* v1end, const Digit* v2begin, const Digit* v2end) noexcept
		{
			//const std::ptrdiff_t rezlen = rezend - rezbegin;
			const std::ptrdiff_t v1len = v1end - v1begin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			assert((rezend - rezbegin) >= (v1end - v1begin) + (v2end - v2begin));
						
			for (std::ptrdiff_t v1it = 0; v1it < v1len; ++v1it)
			{
				unsigned long long int carry = 0;
				for (std::ptrdiff_t v2it = 0; v2it < v2len; ++v2it)
				{
					carry += (unsigned long long int)(rezbegin[v1it + v2it]) + (unsigned long long int)(v1begin[v1it]) * (unsigned long long int)(v2begin[v2it]);
					rezbegin[v1it + v2it] = Digit(carry & DigitMaskBits);
					carry >>= DigitSizeBits;
				}
				rezbegin[v1it + v2len] = Digit(carry & DigitMaskBits);
			}
		}
	}
	
	template<int N>
	class Num
	{
	public:
		enum
		{
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

		explicit Num(int num) noexcept
		{
			for (int i = 0; i < Size; i++)
			{
				data[i] = num & ((1ull << DigitSizeBits) - 1ull);
				num >>= DigitSizeBits;
			}
		}

		template<typename T>
		explicit Num(T num) noexcept
		{
			for (int i = 0; i < Size; i++)
			{
				data[i] = num & ((1ull << DigitSizeBits) - 1ull);
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
				
			//assert(std::none_of(&other[std::min<int>(Size, Num<M>::Size)],&other[Num<M>::Size], [](unsigned char x){ return x != 0; }));

			return *this;
		}

		Digit& operator[](int i) noexcept { return data[i]; }
		const Digit& operator[](int i) const noexcept { return data[i]; }

		bool bit(int i) const noexcept
		{
			assert(i >= 0 && i < N);
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
	const Num<N> operator + (Num<N> v1, const Num<M>& v2) noexcept
	{
		Core::add(v1.begin(), v1.end(), v2.begin(), v2.end());
		return v1;
	}

	template<int N>
	const Num<N> operator + (Num<N> v1, Digit v2) noexcept
	{
		Core::add(v1.begin(), v1.end(), &v2, (&v2) + 1);
		return v1;
	}

	template<int N>
	const Num<N> operator - (const Num<N>& v) noexcept
	{
		return Num<N>(0) - v;
	}

	template<int N, int M>
	const Num<N> operator - (Num<N> v1, const Num<M>& v2) noexcept
	{
		Core::sub(v1.begin(), v1.end(), v2.begin(), v2.end());
		return v1;
	}

	template<int N>
	const Num<N> operator - (Num<N> v1, Digit v2) noexcept
	{
		Core::sub(v1.begin(), v1.end(), &v2, (&v2) + 1);
		return v1;
	}

	template<int N, int M>
	const Num<N+M> operator * (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		Num<N + M> rez(0);
		Core::mul(rez.begin(), rez.end(), v1.begin(), v1.end(), v2.begin(), v2.end());
		return rez;
	}

	template<int N>
	const Num<N + DigitSizeBits> operator * (const Num<N>& v1, Digit v2) noexcept
	{
		Num<N + DigitSizeBits> rez(0);
		Core::mul(rez.begin(), rez.end(), v1.begin(), v1.end(), &v2, (&v2) + 1);
		return rez;
	}

	template<int N, int M>
	const bool operator == (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(&v1[0], &v1[v1.Size], &v2[0], &v2[v2.Size]) == 0;
	}

	template<int N>
	const bool operator == (const Num<N>& v1, Digit v2) noexcept
	{
		return Core::cmp(&v1[0], &v1[v1.Size], &v2, (&v2) + 1) == 0;
	}

	template<int N, int M>
	const bool operator != (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return !(v1 == v2);
	}
	
	template<int N>
	const bool operator != (const Num<N>& v1, Digit v2) noexcept
	{
		return !(v1 == v2);
	}

	template<int N, int M>
	const bool operator > (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(&v1[0], &v1[v1.Size], &v2[0], &v2[v2.Size]) > 0;
	}

	template<int N>
	const bool operator > (const Num<N>& v1, Digit v2) noexcept
	{
		return Core::cmp(&v1[0], &v1[v1.Size], &v2, (&v2) + 1) > 0;
	}

	template<int N, int M>
	const bool operator >= (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(&v1[0], &v1[v1.Size], &v2[0], &v2[v2.Size]) >= 0;
	}

	template<int N>
	const bool operator >= (const Num<N>& v1, Digit v2) noexcept
	{
		return Core::cmp(&v1[0], &v1[v1.Size], &v2, (&v2) + 1) >= 0;
	}

	template<int N, int M>
	const bool operator < (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(&v1[0], &v1[v1.Size], &v2[0], &v2[v2.Size]) < 0;
	}

	template<int N>
	const bool operator < (const Num<N>& v1, Digit v2) noexcept
	{
		return Core::cmp(&v1[0], &v1[v1.Size], &v2, (&v2) + 1) < 0;
	}

	template<int N, int M>
	const bool operator <= (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(&v1[0], &v1[v1.Size], &v2[0], &v2[v2.Size]) <= 0;
	}

	template<int N>
	const bool operator <= (const Num<N>& v1, Digit v2) noexcept
	{
		return Core::cmp(&v1[0], &v1[v1.Size], &v2, (&v2) + 1) <= 0;
	}

	template<int N>
	const Num<N> shift_left(const Num<N>& v, int count) noexcept
	{
		Num<N> ret;
		if (count == 0)
		{
			ret = v;
		}
		else if (count > 0 && count < N)
		{
			const int bytes = count / DigitSizeBits;
			const int bits = count % DigitSizeBits;
			if (bits == 0)
			{
				for (int i = Num<N>::Size - 1; i >= int(bytes); i--)
				{
					ret[i] = v[i - bytes];
				}
			}
			else
			{
				unsigned long long int carry = 0;
				for (int i = bytes; i < Num<N>::Size; i++)
				{
					unsigned long long int val = static_cast<unsigned long long int>(v[i - bytes]) << bits;
					ret[i] = static_cast<Digit>(val | carry);
					carry = val >> DigitSizeBits;
				}
			}
		}
		return ret;
	}

	template<int N>
	const Num<N> operator << (const Num<N>& v1, int bits) noexcept
	{
		return shift_left(v1, bits);
	}

	template<int N>
	const Num<N> shift_right(const Num<N>& v, int count) noexcept
	{
		Num<N> ret;
		if (count == 0)
		{
			ret = v;
		}
		else if (count > 0 && count < N)
		{
			const int bytes = count / DigitSizeBits;
			const int bits = count % DigitSizeBits;
			if (bits == 0)
			{
				for (int i = 0; i < Num<N>::Size - bytes; i++)
				{
					ret[i] = v[i + bytes];
				}
			}
			else
			{
				unsigned long long int carry = v[bytes];
				for (int i = 0; i < Num<N>::Size - bytes; i++)
				{
					unsigned long long int val = (static_cast<unsigned long long int>(((i + bytes + 1) < Num<N>::Size) ? v[i + bytes + 1] : 0) << DigitSizeBits);
					val = (val | carry);
					ret[i] = static_cast<Digit>(val >> bits);
					carry = val >> (DigitSizeBits);
				}
			}
		}
		return ret;
	}

	template<int N>
	const Num<N> operator >> (const Num<N>& v1, int bits) noexcept
	{
		return shift_right(v1, bits);
	}

	template<int N>
	std::array<Num<N + DigitSizeBits>, DigitSizeBits> calcShiftedD(const Num<N>& d) noexcept
	{
		std::array<Num<N + DigitSizeBits>, DigitSizeBits> ret;
		ret[0] = d;
		for (int i = 1; i < DigitSizeBits; i++)
			ret[i] = ret[0] << i;
		return ret;
	}

	template<int N, int M>
	Num<N> div(Num<N>& n, const Num<M>& d) noexcept
	{
		Num<N> q;

		const std::array<Num<M + DigitSizeBits>, DigitSizeBits> shiftedD = calcShiftedD(d);

		Digit* qbeg = n.end()-1;
		Digit* qend = n.end();

		for (int i = Num<N>::Size; i-->0;)
		{
			Digit val = 0;
			if (Core::cmp(qbeg, qend, d.begin(), d.end()) >= 0)
			{
				for (int j = DigitSizeBits; j-- > 0;)
				{
					if (Core::cmp(qbeg, qend, shiftedD[j].begin(), shiftedD[j].end()) >= 0)
					{
						Core::sub(qbeg, qend, shiftedD[j].begin(), shiftedD[j].end());
						val += 1 << j;
					}
				}
			}
			q[i] = val;
			qbeg--;
		}
		return q;
	}

	template<int N, int M>
	Num<M> mod(Num<N> n, const Num<M>& d) noexcept
	{
		const std::array<Num<M + DigitSizeBits>, DigitSizeBits> shiftedD = calcShiftedD(d);

		Digit* qbeg = n.end() - 1;
		Digit* qend = n.end();

		for (int i = Num<N>::Size; i-->0;)
		{
			if (Core::cmp(qbeg, qend, d.begin(), d.end()) >= 0)
			{
				for (int j = DigitSizeBits; j-- > 0;)
				{
					if (Core::cmp(qbeg, qend, shiftedD[j].begin(), shiftedD[j].end()) >= 0)
					{
						Core::sub(qbeg, qend, shiftedD[j].begin(), shiftedD[j].end());
					}
				}
			}
			qbeg--;
		}
		return n;
	}

	template<int N>
	Num<N> div(Num<N>& n, const Digit d) noexcept
	{
		Num<N> q(0);
		unsigned long long int rest = 0;
		for (int i = Num<N>::Size; i--> 0;)
		{
			rest = rest << DigitSizeBits;
			rest += n[i];
			n[i] = 0;
			q[i] = (Digit)(rest / d);
			rest = rest % d;
		}
		n[0] = (Digit)rest;
		return q;
	}

	template<int N>
	Digit mod(const Num<N>& n, const Digit d) noexcept
	{
		unsigned long long int rest = 0;
		for (int i = Num<N>::Size; i--> 0;)
		{
			rest = rest << DigitSizeBits;
			rest += n[i];
			rest = rest % d;
		}
		return (Digit)rest;
	}

	template<int N, int M>
	const Num<N> operator / (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		if (v1 < v2)
			return Num<N>(0);

		Num<N> q = v1;
		return div(q, v2);
	}

	template<int N>
	const Num<N> operator / (const Num<N>& v1, const Digit v2) noexcept
	{
		Num<N> q = v1;
		return div(q, v2);
	}

	template<int N, int M>
	const Num<M> operator % (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return mod(v1, v2);
	}

	template<int N>
	const Digit operator % (const Num<N>& v1, const Digit v2) noexcept
	{
		return mod(v1, v2);
	}

	template<int N, int M>
	const Num<N> operator & (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		Num<N> rez(0);

		for (int i = 0; i < Num<N>::Size; i++)
		{
			rez[i] = v1[i] & ((i < M) ? v2[i] : 0);
		}

		return rez;
	}

	template<int N, int M>
	const Num<N> operator & (const Num<N>& v1, Num<M>&& v2) noexcept
	{
		Num<N> rez(0);

		for (int i = 0; i < Num<N>::Size; i++)
		{
			rez[i] = v1[i] & ((i < M) ? v2[i] : 0);
		}

		return rez;
	}

	template<int N>
	const Num<N> operator & (const Num<N>& v1, const Digit v2) noexcept
	{
		Num<N> rez(v1[0] & v2);
		return rez;
	}

	template<int N>
	const Num<N> operator | (const Num<N>& v1, const Num<N>& v2) noexcept
	{
		Num<N> rez;

		for (int i = 0; i < Num<N>::Size; i++)
		{
			rez[i] = v1[i] | v2[i];
		}

		return rez;
	}

	template<int N>
	const Num<N> operator | (const Num<N>& v1, const Digit v2) noexcept
	{
		Num<N> rez = v1;
		rez[0] = v1[0] | v2;
		return rez;
	}

	template<int N>
	const Num<N> operator ~ (const Num<N>& v) noexcept
	{
		Num<N> rez;

		for (int i = 0; i < Num<N>::Size; i++)
		{
			rez[i] = ~v[i];
		}

		return rez;
	}

	template<int N, int M>
	bool operator && (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return v1 != 0 && v2 != 0;
	}

	template<int N, int M>
	bool operator || (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return v1 != 0 || v2 != 0;
	}

	template<int N>
	bool operator !(const Num<N>& v) noexcept
	{
		return v == 0;
	}
	
	template<int N>
	const Num<N> gcd(const Num<N>& v1, const Num<N>& v2) noexcept
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
	const Num<2*N> lcm(const Num<N>& v1, const Num<N>& v2) noexcept
	{
		Num<2 * N> m = v1 * v2;
		Num<N> g = gcd(v1, v2);
		return m / g;
	}

	template<int N>
	const Num<N> modInv(const Num<N>& a, const Num<N>& n) noexcept
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
	    if (r > Num<N>(1)) return Num<N>(0);
		if (t > n)
		{
			// signed < 0
			t = t + n;
		}
	    return t;
	}

	template<typename T>
	const T modInv(const T& a, const T& n) noexcept
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
		if (r > T(1)) return T(0);
		if (t < 0) t = t + n;
		return t;
	}

	template<int N>
	const Num<N> modExp(const Num<N>& base, const Num<N>& exp, const Num<N>& modulo) noexcept
	{
		if (modulo == 1)
			return Num<N>(0);

		Num<N> result(1);
		Num<N> power = base;

		int realExpSize = Num<N>::Size;
		// skip leading zeros
		while (realExpSize >= 0 && exp[realExpSize - 1] == 0)
			realExpSize--;

		for (int i = realExpSize * DigitSizeBits; i-->0;)
		{
			result = (result * result) % modulo;
			if (exp.bit(i))
				result = (result * power) % modulo;
		}

		return result;
	}

	template<int N>
	struct MonMul
	{
		const Num<N> n;
		Num<N> rinv;
		Num<N> ninv;
		//Num<N> rmask;

		MonMul(const Num<N>& _n) noexcept : n(_n)
		{
			const Num<2 * N> r(Num<2 * N>(1) << N);

			assert(gcd(Num<2 * N>(n), r) == 1);

			rinv = Num<N>(modInv<2 * N>(r, n));

			assert((Num<2 * N>(rinv) << N) % Num<2 * N>(n) == 1);

			ninv = ((Num<2 * N>(rinv) << N) - 1) / n;

			assert((Num<2 * N>(rinv) << N) - (ninv*n) == 1);

			//rmask = (r - 1);
		}

		const Num<N> In(const Num<N>& x) const noexcept
		{
			Num<2 * N> ret(x);
			ret = ret << N;

			return ret % n;
		}

		template<int M>
		const Num<N> In(const Num<M>& x) const noexcept
		{
			Num<N+M> ret(x);
			ret = ret << N;

			return ret % n;
		}

		const Num<N> Out(const Num<N>& x) const noexcept
		{
			return (*this)(x, Num<N>(1));
		}

		const Num<N> operator()(const Num<N>& a, const Num<N>& b) const noexcept
		{
			Num<2*N> t = a * b;
			//Num<N> m = (Num<N>(t & rmask) * ninv) & rmask;
			Num<N> m = Num<N>(t) * ninv;
			Num<N> u = (t + m* n) >> N;

			//Num<2 * N> t1 = t;
			//Core::mul(t1.begin(), t1.end(), m.begin(), m.end(), n.begin(), n.end());
			//Num<N> u1 = (t1) >> N;
			//assert(u == u1);

			while (u >= n)
				u = u - n;
			return Num<N>(u);
		}
	};

	template<int N>
	const Num<N> monModMul(const Num<N>& a, const Num<N>& b, const Num<N>& mod) noexcept
	{
		MonMul<N> monMod((mod));
		Num<N> xa = monMod.In(a);
		return monMod(xa, b);
	}

	template<int N, int M, int K>
	const Num<N> monModExp(const Num<K>& base, const Num<M>& exp, const Num<N>& mod) noexcept
	{
		MonMul<N> monMod((mod));
		Num<N> ret = monMod.In(Num<N>(1));
		Num<N> power = monMod.In(base);

		int realExpSize = Num<M>::Size;
		while (realExpSize > 0 && exp[realExpSize - 1] == 0)
			realExpSize--;

		for (int i = realExpSize * DigitSizeBits; i --> 0;)
		{
			ret = monMod(ret, ret);
			if (exp.bit(i))
				ret = monMod(ret, power);
		}
		return monMod.Out(ret);
	}

	template<int N, int M, int K>
	const Num<N> monModExp2ary(const Num<K>& base, const Num<M>& exp, const Num<N>& mod) noexcept
	{
		MonMul<N> monMod((mod));
		Num<N> ret = monMod.In(Num<N>(1));
		Num<N> baseIn = monMod.In(base);
		Num<N> power[3];
		power[0] = baseIn;
		power[1] = monMod(baseIn, baseIn);
		power[2] = monMod(power[1], baseIn);

		int realExpSize = Num<M>::Size;
		while (realExpSize > 0 && exp[realExpSize - 1] == 0)
			realExpSize--;

		for (int i = realExpSize * DigitSizeBits/2; i--> 0;)
		{
			ret = monMod(ret, ret);
			ret = monMod(ret, ret);
			int p = (exp.bit(i * 2 + 1) ? 1 : 0) * 2 + (exp.bit(i * 2) ? 1 : 0);
			if (p > 0)
				ret = monMod(ret, power[p - 1]);
		}
		return monMod.Out(ret);
	}

	template<int N>
	const Num<N> rand() noexcept
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
	bool millerRabinPass(const BigNum::Num<N>& num, const BigNum::Num<N>& a) noexcept
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

			a_to_power = (a_to_power*a_to_power) % num;
		}

		if (a_to_power == num - 1)
			return true;

		return false;
	}

	//numbers from HAC table 4.3
	int millerRabinProbes(int N)
	{
		if (N>=600) return 2; 
		if (N>=550) return 4;
		if (N>=500) return 5;
		if (N>=400) return 6;
		if (N>=350) return 7;
		if (N>=300) return 9;
		if (N>=250) return 12;
		if (N>=200) return 15;
		if (N>=150) return 18;
		if (N>=100) return 27;
		            return 40;
	}

	template<int N>
	bool millerRabinTest(const BigNum::Num<N>& num) noexcept
	{
		int times = millerRabinProbes(N);
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
		947,    953,    967,    971,    977,    983,    991,    997,   1009,   1013, 
	   1019,   1021,   1031,   1033,   1039,   1049,   1051,   1061,   1063,   1069, 
	   1087,   1091,   1093,   1097,   1103,   1109,   1117,   1123,   1129,   1151, 
	   1153,   1163,   1171,   1181,   1187,   1193,   1201,   1213,   1217,   1223, 
	   1229,   1231,   1237,   1249,   1259,   1277,   1279,   1283,   1289,   1291, 
	   1297,   1301,   1303,   1307,   1319,   1321,   1327,   1361,   1367,   1373, 
	   1381,   1399,   1409,   1423,   1427,   1429,   1433,   1439,   1447,   1451, 
	   1453,   1459,   1471,   1481,   1483,   1487,   1489,   1493,   1499,   1511, 
	   1523,   1531,   1543,   1549,   1553,   1559,   1567,   1571,   1579,   1583, 
	   1597,   1601,   1607,   1609,   1613,   1619,   1621,   1627,   1637,   1657, 
	   1663,   1667,   1669,   1693,   1697,   1699,   1709,   1721,   1723,   1733, 
	   1741,   1747,   1753,   1759,   1777,   1783,   1787,   1789,   1801,   1811, 
	   1823,   1831,   1847,   1861,   1867,   1871,   1873,   1877,   1879,   1889, 
	   1901,   1907,   1913,   1931,   1933,   1949,   1951,   1973,   1979,   1987, 
	   1993,   1997,   1999,   2003,   2011,   2017,   2027,   2029,   2039,   2053,
	   2063,   2069,   2081,   2083,   2087,   2089,   2099,   2111,   2113,   2129, 
	   2131,   2137,   2141,   2143,   2153,   2161,   2179,   2203,   2207,   2213, 
	   2221,   2237,   2239,   2243,   2251,   2267,   2269,   2273,   2281,   2287, 
	   2293,   2297,   2309,   2311,   2333,   2339,   2341,   2347,   2351,   2357, 
	   2371,   2377,   2381,   2383,   2389,   2393,   2399,   2411,   2417,   2423, 
	   2437,   2441,   2447,   2459,   2467,   2473,   2477,   2503,   2521,   2531, 
	   2539,   2543,   2549,   2551,   2557,   2579,   2591,   2593,   2609,   2617, 
	   2621,   2633,   2647,   2657,   2659,   2663,   2671,   2677,   2683,   2687, 
	   2689,   2693,   2699,   2707,   2711,   2713,   2719,   2729,   2731,   2741, 
	   2749,   2753,   2767,   2777,   2789,   2791,   2797,   2801,   2803,   2819, 
	   2833,   2837,   2843,   2851,   2857,   2861,   2879,   2887,   2897,   2903, 
	   2909,   2917,   2927,   2939,   2953,   2957,   2963,   2969,   2971,   2999,};

	constexpr int simplechecks(int size) noexcept
	{
		//return sizeof(primes) / sizeof(primes[0]);
		return (size <= 256) ? 30 : ((size <= 512) ? 50 : ((size <= 1024) ? 80 : ((size <= 2048) ? 128 : (sizeof(primes) / sizeof(primes[0])))));
	}

	template<int N>
	bool simpleTest(const BigNum::Num<N>& num)
	{
		for (int i = 0; i < sizeof(primes) / sizeof(primes[0]); i++)
		{
			if (num % primes[i] == 0)
				return false;
		}
		return true;
	}

	template<int N>
	const Num<N> nextPrimeOpt(const BigNum::Num<N>& from) noexcept
	{
		BigNum::Num<N> prime = from | 1;

		unsigned int rests[simplechecks(N)] = { 0 };

		bool simpleTestPassed = true;
		for (int i = 0; i < sizeof(rests) / sizeof(rests[0]); i++)
		{
			rests[i] = prime % primes[i];
			simpleTestPassed = simpleTestPassed && (rests[i] != 0);
		}

		while (!(simpleTestPassed && millerRabinTest(prime)))
		{
			simpleTestPassed = true;
			for (int i = 0; i < sizeof(rests) / sizeof(rests[0]); i++)
			{
				rests[i] = (rests[i] + 2) % primes[i];
				simpleTestPassed = simpleTestPassed && (rests[i] != 0);
			}
			prime = prime + 2;
		}

		return prime;
	}

	template<int N>
	const Num<N> nextPrime(const BigNum::Num<N>& from) noexcept
	{
		BigNum::Num<N> prime = from | 1;

		while (!(simpleTest(prime) && millerRabinTest(prime)))
		{
			prime = prime + 2;
		}

		return prime;
	}

	enum class RandType
	{
		Simple,
		FindNext,
		FindNextOpt,
	};

	template<int N>
	const Num<N> randPrime(RandType type = RandType::Simple) noexcept
	{
		static BigNum::Num<N> mask = (BigNum::Num<N>(1) | (BigNum::Num<N>(1) << (N - 1)));
		BigNum::Num<N> num = rand<N>() | mask;
		switch(type)
		{
		case RandType::FindNext:
			return nextPrime(num);
		case RandType::FindNextOpt:
			return nextPrimeOpt(num);
		case RandType::Simple:
		default:
			{
				while (!(simpleTest(num) && millerRabinTest(num)))
				{
					num = rand<N>() | mask;
				}
				return num;
			}
		}
		//return nextPrime(num);
	}

	namespace operators
	{
		template <char... args>
		BigNum::Num<2048> operator "" _bn2048()
		{
			BigNum::Num<2048> ret(0);
			char buf[sizeof...(args)] = { args... };
			for (char i : buf)
			{
				assert(isdigit(i));
				ret = ret * 10 + BigNum::Digit(i - '0');
			}
			return ret;
		}

		template <char... args>
		BigNum::Num<1024> operator "" _bn1024()
		{
			BigNum::Num<1024> ret(0);
			char buf[sizeof...(args)] = { args... };
			for (char i : buf)
			{
				assert(isdigit(i));
				ret = ret * 10 + BigNum::Digit(i - '0');
			}
			return ret;
		}
		
		template <char... args>
		BigNum::Num<512> operator "" _bn512()
		{
			BigNum::Num<512> ret(0);
			char buf[sizeof...(args)] = { args... };
			for (char i : buf)
			{
				assert(isdigit(i));
				ret = ret * 10 + BigNum::Digit(i - '0');
			}
			return ret;
		}
		
		template <char... args>
		BigNum::Num<256> operator "" _bn256()
		{
			BigNum::Num<256> ret(0);
			char buf[sizeof...(args)] = { args... };
			for (char i : buf)
			{
				assert(isdigit(i));
				ret = ret * 10 + BigNum::Digit(i - '0');
			}
			return ret;
		}
		
		template <char... args>
		BigNum::Num<128> operator "" _bn128()
		{
			BigNum::Num<128> ret(0);
			char buf[sizeof...(args)] = { args... };
			for (char i : buf)
			{
				assert(isdigit(i));
				ret = ret * 10 + BigNum::Digit(i - '0');
			}
			return ret;
		}
		
		template <char... args>
		BigNum::Num<64> operator "" _bn64()
		{
			BigNum::Num<64> ret(0);
			char buf[sizeof...(args)] = { args... };
			for (char i : buf)
			{
				assert(isdigit(i));
				ret = ret * 10 + BigNum::Digit(i - '0');
			}
			return ret;
		}
	}
};

using namespace BigNum::operators;