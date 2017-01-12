#pragma once

#include <algorithm>
#include <array>

namespace BigNum
{
	typedef unsigned int Word;
	typedef unsigned long long DWord;
	constexpr size_t WordSizeBits = sizeof(Word) * 8;
	constexpr DWord WordMaskBits = ((1ull << WordSizeBits) - 1ull);

	static_assert(sizeof(DWord) > sizeof(Word), "cant detect overflow");

	namespace Core
	{
		int cmp(const Word* v1begin, const Word* v1end, const Word* v2begin, const Word* v2end) noexcept
		{
			const std::ptrdiff_t v1len = v1end - v1begin;
			const std::ptrdiff_t v2len = v2end - v2begin;
			const std::ptrdiff_t vmax = (v1len > v2len)? v1len : v2len;

			for(std::ptrdiff_t i = vmax; i-->0;)
			{
				const Word v1 = (i < v1len) ? v1begin[i] : 0;
				const Word v2 = (i < v2len) ? v2begin[i] : 0;
				if (v1 != v2)
					return v1 > v2 ? 1 : -1;
			}
			return 0;
		}

        template<int N, int M>
      	int cmp(const Word (&v1)[N], const Word (&v2)[M]) noexcept
		{
			const int v1len = N;
			const int v2len = M;
			const int vmax = (v1len > v2len)? v1len : v2len;

			for(int i = vmax; i-->0;)
			{
				const Word v1v = (i < v1len) ? v1[i] : 0;
				const Word v2v = (i < v2len) ? v2[i] : 0;
				if (v1v != v2v)
					return v1v > v2v ? 1 : -1;
			}
			return 0;
		}

		Word adc(Word& ret, const Word v1, const Word v2, const Word carry) noexcept
		{
			const DWord sum = DWord(v1) + v2 + carry;
			ret = sum & WordMaskBits;
			return Word(sum >> WordSizeBits);
		}

		Word add(Word* rezbegin, Word* rezend, const Word* v1begin, const Word* v1end, const Word* v2begin, const Word* v2end) noexcept
		{
			const std::ptrdiff_t rezlen = rezend - rezbegin;
			const std::ptrdiff_t v1len = v1end - v1begin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			Word carry = 0;
			for (std::ptrdiff_t i = 0; i < rezlen; i++)
			{
				const Word v1 = (i < v1len) ? v1begin[i] : 0;
				const Word v2 = (i < v2len) ? v2begin[i] : 0;
				carry = adc(rezbegin[i], v1, v2, carry);
			}
			return carry;
		}
		
		template<int N, int M, int K>
		Word add(Word (&rez)[N], const Word (&v1)[M], const Word (&v2)[K]) noexcept
		{
			Word carry = 0;
			for (int i = 0; i < N; i++)
			{
				const Word v1v = (i < M) ? v1[i] : 0;
				const Word v2v = (i < K) ? v2[i] : 0;
				carry = adc(rez[i], v1v, v2v, carry);
			}
			return carry;
		}

		Word add(Word* rezbegin, Word* rezend, const Word* v2begin, const Word* v2end) noexcept
		{
			const std::ptrdiff_t rezlen = rezend - rezbegin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			Word carry = 0;
			for (std::ptrdiff_t i = 0; i < rezlen; i++)
			{
				const Word v2 = (i < v2len) ? v2begin[i] : 0;
				carry = adc(rezbegin[i], rezbegin[i], v2, carry);
			}
			return carry;
		}

		Word sbc(Word& ret, const Word v1, const Word v2, const Word carry) noexcept
		{
			const DWord sum = DWord(v1) - v2 - carry;
			ret = sum & WordMaskBits;
			return (sum >> WordSizeBits) ? 1 : 0;
		}

		Word sub(Word* rezbegin, Word* rezend, const Word* v1begin, const Word* v1end, const Word* v2begin, const Word* v2end) noexcept
		{
			const std::ptrdiff_t rezlen = rezend - rezbegin;
			const std::ptrdiff_t v1len = v1end - v1begin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			Word carry = 0;
			for (std::ptrdiff_t i = 0; i < rezlen; i++)
			{
				const Word v1 = (i < v1len) ? v1begin[i] : 0;
				const Word v2 = (i < v2len) ? v2begin[i] : 0;
				carry = sbc(rezbegin[i], v1, v2, carry);
			}
			return carry;
		}

		template<int N, int M, int K>
		Word sub(Word (&rez)[N], const Word (&v1)[M], const Word (&v2)[K]) noexcept
		{
			Word carry = 0;
			for (int i = 0; i < N; i++)
			{
				const Word v1v = (i < M) ? v1[i] : 0;
				const Word v2v = (i < K) ? v2[i] : 0;
				carry = sbc(rez[i], v1v, v2v, carry);
			}
			return carry;
		}

		Word sub(Word* rezbegin, Word* rezend, const Word* v2begin, const Word* v2end) noexcept
		{
			const std::ptrdiff_t rezlen = rezend - rezbegin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			Word carry = 0;
			for (std::ptrdiff_t i = 0; i < rezlen; i++)
			{
				const Word v2 = (i < v2len) ? v2begin[i] : 0;
				carry = sbc(rezbegin[i], rezbegin[i], v2, carry);
			}
			return carry;
		}

		void mul(Word* rezbegin, Word* rezend, const Word* v1begin, const Word* v1end, const Word* v2begin, const Word* v2end) noexcept
		{
			const std::ptrdiff_t v1len = v1end - v1begin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			assert((rezend - rezbegin) >= (v1end - v1begin) + (v2end - v2begin));
			assert(std::any_of(rezbegin, rezend, [&](Word w) { return w > 0; }) == false);

			for (std::ptrdiff_t v1it = 0; v1it < v1len; ++v1it)
			{
				DWord carry = 0;
				for (std::ptrdiff_t v2it = 0; v2it < v2len; ++v2it)
				{
					carry += DWord(rezbegin[v1it + v2it]) + DWord(v1begin[v1it]) * DWord(v2begin[v2it]);
					rezbegin[v1it + v2it] = Word(carry & WordMaskBits);
					carry >>= WordSizeBits;
				}
				rezbegin[v1it + v2len] = Word(carry & WordMaskBits);
			}
		}
		
		template<int N, int M, int K>
        void mul(Word (&rez)[N], const Word (&v1)[M], const Word (&v2)[K]) noexcept
		{
			static_assert((N) >= (M) + (K), "");
			assert(std::any_of(std::begin(rez), std::end(rez), [&](Word w) { return w > 0; }) == false);

			for (int v1it = 0; v1it < M; ++v1it)
			{
				DWord carry = 0;
				for (int v2it = 0; v2it < K; ++v2it)
				{
					carry += DWord(rez[v1it + v2it]) + DWord(v1[v1it]) * DWord(v2[v2it]);
					rez[v1it + v2it] = Word(carry & WordMaskBits);
					carry >>= WordSizeBits;
				}
				rez[v1it + K] = Word(carry & WordMaskBits);
			}
		}
		
		template<unsigned int N>
		void shl(Word (&rez)[N], const Word (&v)[N], unsigned int count) noexcept
		{
			assert(count < N*sizeof(Word)*8);

			const int bytes = count / WordSizeBits;
			const int bits = count % WordSizeBits;
			if (bits == 0)
			{
				for (int i = N; i --> bytes;)
				{
					rez[i] = v[i - bytes];
				}
			}
			else
			{
				int i = 0;
				for (; i < bytes; i++)
					rez[i] = 0;
				for (; i < N; i++)
					rez[i] = (v[i - bytes] << (bits)) | ( ((i - bytes - 1 >= 0) ? v[i - bytes - 1] : 0) >> (sizeof(Word) * 8 - bits));
			}
		}
		
		template<int N>
		void shr(Word (&rez)[N], const Word (&v)[N], unsigned int count) noexcept
		{
			assert(count < N*sizeof(Word)*8);
			
			const int bytes = count / WordSizeBits;
			const int bits = count % WordSizeBits;
			if (bits == 0)
			{
				for (int i = 0; i < N - bytes; i++)
					rez[i] = v[i + bytes];
			}
			else
			{
				int i = 0;
				for (; i < N - bytes; i++)
					rez[i] = (((i + bytes + 1 < N) ? v[i + bytes + 1] : 0) << (sizeof(Word) * 8 - bits)) | (v[i + bytes] >> bits);
				for (; i < N; i++)
					rez[i] = 0;
			}
		}
	}
	
	template<int N>
	class Num
	{
	public:
		enum
		{
			Size = (N + (WordSizeBits) - 1) / (WordSizeBits),
		};

		Word data[Size];

		Num() noexcept
		{
			for (int i = 0; i < Size; i++)
				data[i] = 0;
		}

		explicit Num(int num) noexcept
		{
			for (int i = 0; i < Size; i++)
			{
				data[i] = num & ((1ull << WordSizeBits) - 1ull);
				num >>= WordSizeBits/2;
				num >>= WordSizeBits/2;
			}
		}

		template<typename T>
		explicit Num(T num) noexcept
		{
			for (int i = 0; i < Size; i++)
			{
				data[i] = num & ((1ull << WordSizeBits) - 1ull);
				num >>= WordSizeBits/2;
				num >>= WordSizeBits/2;
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

		Word& operator[](int i) noexcept { return data[i]; }
		const Word& operator[](int i) const noexcept { return data[i]; }

		bool bit(int i) const noexcept
		{
			assert(i >= 0 && i < N);
			bool ret = false;
			if (i >= 0 && i < N)
			{
				const int byte = i / WordSizeBits;
				const int bit = i % WordSizeBits;
				ret = ((*this)[byte] & (1 << bit)) != 0;
			}
			return ret;
		}

		Word* begin() noexcept { return &data[0]; }
		Word* end() noexcept { return &data[Size]; }

		const Word* begin() const noexcept { return &data[0]; }
		const Word* end() const noexcept { return &data[Size]; }
	};

	template<int N, int M>
	const Num<N> operator + (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		Num<N> ret;
		Core::add(ret.data, v1.data, v2.data);
		return ret;
	}

	template<int N>
	const Num<N> operator + (const Num<N>& v1, Word v2) noexcept
	{
		Num<N> ret;
		const Word a2[] = {v2};
		Core::add(ret.data, v1.data, a2);
		return ret;
	}

	template<int N>
	const Num<N> operator - (const Num<N>& v) noexcept
	{
		return Num<N>(0) - v;
	}

	template<int N, int M>
	const Num<N> operator - (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		Num<N> ret;
		Core::sub(ret.data, v1.data, v2.data);
		return ret;
	}

	template<int N>
	const Num<N> operator - (const Num<N>& v1, Word v2) noexcept
	{
		Num<N> ret;
		const Word a2[] = { v2 };
		Core::sub(ret.data, v1.data, a2);
		return ret;
	}

	template<int N, int M>
	const Num<N+M> operator * (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		Num<N + M> rez(0);
		Core::mul(rez.data, v1.data, v2.data);
		return rez;
	}

	template<int N>
	const Num<N + WordSizeBits> operator * (const Num<N>& v1, Word v2) noexcept
	{
		Num<N + WordSizeBits> rez(0);
		const Word a2[] = { v2 };
		Core::mul(rez.data, v1.data, a2);
		return rez;
	}

	template<int N, int M>
	const bool operator == (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(v1.data, v2.data) == 0;
	}

	template<int N>
	const bool operator == (const Num<N>& v1, Word v2) noexcept
	{
		const Word a2[] = { v2 };
		return Core::cmp(v1.data, a2) == 0;
	}

	template<int N, int M>
	const bool operator != (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return !(v1 == v2);
	}
	
	template<int N>
	const bool operator != (const Num<N>& v1, Word v2) noexcept
	{
		return !(v1 == v2);
	}

	template<int N, int M>
	const bool operator > (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(v1.data, v2.data) > 0;
	}

	template<int N>
	const bool operator > (const Num<N>& v1, Word v2) noexcept
	{
		const Word a2[] = { v2 };
		return Core::cmp(v1.data, a2) > 0;
	}

	template<int N, int M>
	const bool operator >= (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(v1.data, v2.data) >= 0;
	}

	template<int N>
	const bool operator >= (const Num<N>& v1, Word v2) noexcept
	{
		const Word a2[] = { v2 };
		return Core::cmp(v1.data, a2) >= 0;
	}

	template<int N, int M>
	const bool operator < (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(v1.data, v2.data) < 0;
	}

	template<int N>
	const bool operator < (const Num<N>& v1, Word v2) noexcept
	{
		const Word a2[] = { v2 };
		return Core::cmp(v1.data, a2) < 0;
	}

	template<int N, int M>
	const bool operator <= (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(v1.data, v2.data) <= 0;
	}

	template<int N>
	const bool operator <= (const Num<N>& v1, Word v2) noexcept
	{
		const Word a2[] = { v2 };
		return Core::cmp(v1.data, a2) <= 0;
	}

	template<int N>
	const Num<N> operator << (const Num<N>& v1, unsigned int bits) noexcept
	{
		Num<N> ret;
		if (bits < N)
			Core::shl(ret.data, v1.data, bits);
		return ret;
	}

	template<int N>
	const Num<N> operator >> (const Num<N>& v1, unsigned int bits) noexcept
	{
		Num<N> ret;
		if (bits < N)
			Core::shr(ret.data, v1.data, bits);
		return ret;
	}

	template<int N>
	std::array<Num<N + WordSizeBits>, WordSizeBits> calcShiftedD(const Num<N>& d) noexcept
	{
		std::array<Num<N + WordSizeBits>, WordSizeBits> ret;
		ret[0] = d;
		for (unsigned int i = 1; i < WordSizeBits; i++)
			ret[i] = ret[0] << i;
		return ret;
	}

	template<int N, int M>
	Num<N> div(Num<N>& n, const Num<M>& d) noexcept
	{
		Num<N> q;

		const std::array<Num<M + WordSizeBits>, WordSizeBits> shiftedD = calcShiftedD(d);

		Word* qbeg = n.end();
		Word* qend = n.end();

		for (int i = Num<N>::Size; i-->0;)
		{
			qbeg--;
			Word val = 0;
			if (Core::cmp(qbeg, qend, d.begin(), d.end()) >= 0)
			{
				for (int j = WordSizeBits; j-- > 0;)
				{
					if (Core::cmp(qbeg, qend, shiftedD[j].begin(), shiftedD[j].end()) >= 0)
					{
						Core::sub(qbeg, qend, shiftedD[j].begin(), shiftedD[j].end());
						val += 1 << j;
					}
				}
			}
			q[i] = val;
		}
		return q;
	}

	template<int N, int M>
	Num<M> mod(Num<N> n, const Num<M>& d) noexcept
	{
		const std::array<Num<M + WordSizeBits>, WordSizeBits> shiftedD = calcShiftedD(d);

		Word* qbeg = n.end();
		Word* qend = n.end();

		for (int i = Num<N>::Size; i-->0;)
		{
			qbeg--;
			if (Core::cmp(qbeg, qend, d.begin(), d.end()) >= 0)
			{
				for (int j = WordSizeBits; j-- > 0;)
				{
					if (Core::cmp(qbeg, qend, shiftedD[j].begin(), shiftedD[j].end()) >= 0)
					{
						Core::sub(qbeg, qend, shiftedD[j].begin(), shiftedD[j].end());
					}
				}
			}
		}
		assert(qbeg == n.begin());
		return n;
	}

	template<int N>
	Num<N> div(Num<N>& n, const Word d) noexcept
	{
		Num<N> q(0);
		DWord rest = 0;
		for (int i = Num<N>::Size; i--> 0;)
		{
			rest = rest << WordSizeBits;
			rest += n[i];
			n[i] = 0;
			q[i] = (Word)(rest / d);
			rest = rest % d;
		}
		n[0] = (Word)rest;
		return q;
	}

	template<int N>
	Word mod(const Num<N>& n, const Word d) noexcept
	{
		DWord rest = 0;
		for (int i = Num<N>::Size; i--> 0;)
		{
			rest = rest << WordSizeBits;
			rest += n[i];
			rest = rest % d;
		}
		return (Word)rest;
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
	const Num<N> operator / (const Num<N>& v1, const Word v2) noexcept
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
	const Word operator % (const Num<N>& v1, const Word v2) noexcept
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
	const Num<N> operator & (const Num<N>& v1, const Word v2) noexcept
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
	const Num<N> operator | (const Num<N>& v1, const Word v2) noexcept
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

		for (int i = realExpSize * WordSizeBits; i-->0;)
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

		explicit MonMul(const Num<N>& _n) noexcept : n(_n)
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
			Num<N> m = Num<N>(t) * ninv;
			Num<N> u = (t + m * n) >> N;

			//Num<2 * N> t1 = t;
			//Core::muladd(t1.data, m.data, n.data);
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

		for (int i = realExpSize * WordSizeBits; i --> 0;)
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

		for (int i = realExpSize * WordSizeBits/2; i--> 0;)
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
	void CRT_init(const Num<N>& p, const Num<N>& q, const Num<2*N>& d, Num<N>& dp, Num<N>& dq, Num<N>& qinvp)
	{
		qinvp = modInv(q, p);
		dp = d % (p - 1);
		dq = d % (q - 1);
	}

	template<int N>
	const Num<2*N> CRT(const Num<2*N>& m, const Num<N>& p, const Num<N>& q, const Num<N>& dp, const Num<N>& dq, const Num<N>& qinvp) noexcept
	{
		BigNum::Num<N> mp = BigNum::monModExp2ary(m, dp, p);
		BigNum::Num<N> mq = BigNum::monModExp2ary(m, dq, q);

		return BigNum::Num<2*N>(mq) + (((mp - mq) * qinvp) % p) * q;
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
	int millerRabinProbes(int N) noexcept
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
	   1019,   1021,   /*1031,   1033,   1039,   1049,   1051,   1061,   1063,   1069, 
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
	   2909,   2917,   2927,   2939,   2953,   2957,   2963,   2969,   2971,   2999,
       3001,   3011,   3019,   3023,   3037,   3041,   3049,   3061,   3067,   3079, 
	   3083,   3089,   3109,   3119,   3121,   3137,   3163,   3167,   3169,   3181, 
	   3187,   3191,   3203,   3209,   3217,   3221,   3229,   3251,   3253,   3257, 
	   3259,   3271,   3299,   3301,   3307,   3313,   3319,   3323,   3329,   3331, 
	   3343,   3347,   3359,   3361,   3371,   3373,   3389,   3391,   3407,   3413, 
	   3433,   3449,   3457,   3461,   3463,   3467,   3469,   3491,   3499,   3511, 
	   3517,   3527,   3529,   3533,   3539,   3541,   3547,   3557,   3559,   3571, 
	   3581,   3583,   3593,   3607,   3613,   3617,   3623,   3631,   3637,   3643, 
	   3659,   3671,   3673,   3677,   3691,   3697,   3701,   3709,   3719,   3727, 
	   3733,   3739,   3761,   3767,   3769,   3779,   3793,   3797,   3803,   3821, 
	   3823,   3833,   3847,   3851,   3853,   3863,   3877,   3881,   3889,   3907, 
	   3911,   3917,   3919,   3923,   3929,   3931,   3943,   3947,   3967,   3989,
	   4001,   4003,   4007,   4013,   4019,   4021,   4027,   4049,   4051,   4057, 
	   4073,   4079,   4091,   4093,   4099,   4111,   4127,   4129,   4133,   4139, 
	   4153,   4157,   4159,   4177,   4201,   4211,   4217,   4219,   4229,   4231, 
	   4241,   4243,   4253,   4259,   4261,   4271,   4273,   4283,   4289,   4297, 
	   4327,   4337,   4339,   4349,   4357,   4363,   4373,   4391,   4397,   4409, 
	   4421,   4423,   4441,   4447,   4451,   4457,   4463,   4481,   4483,   4493, 
	   4507,   4513,   4517,   4519,   4523,   4547,   4549,   4561,   4567,   4583, 
	   4591,   4597,   4603,   4621,   4637,   4639,   4643,   4649,   4651,   4657, 
	   4663,   4673,   4679,   4691,   4703,   4721,   4723,   4729,   4733,   4751, 
	   4759,   4783,   4787,   4789,   4793,   4799,   4801,   4813,   4817,   4831, 
	   4861,   4871,   4877,   4889,   4903,   4909,   4919,   4931,   4933,   4937, 
	   4943,   4951,   4957,   4967,   4969,   4973,   4987,   4993,   4999,   5003,*/
	};

	constexpr unsigned int primesCount = sizeof(primes) / sizeof(primes[0]);

	constexpr unsigned int simplechecks(unsigned int size) noexcept
	{
		return ((size / 16) < primesCount) ? (size / 16) : primesCount;
	}

	template<int N>
	bool simpleTest(const BigNum::Num<N>& num, unsigned skip = 0) noexcept
	{
		for (unsigned int i = skip; i < simplechecks(N); i++)
		{
			if (num % primes[i] == 0)
				return false;
		}
		return true;
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

	template<int N>
	const Num<N> nextPrimeOpt3(const BigNum::Num<N>& from) noexcept
	{
		BigNum::Num<N> prime = from | 1;
		
		constexpr unsigned char steps3[3] = {2, 4, 2};

		unsigned char rest3 = prime % 3;
		if (rest3 == 0)
		{
			rest3 = 2;
			prime = prime + 2;
		}

		while (!(simpleTest(prime, 1) && millerRabinTest(prime)))
		{
			assert(prime % 3 != 0);
			unsigned char step = steps3[rest3];
			rest3 = (rest3 + step) % 3;
			prime = prime + step;
		}

		return prime;
	}

	template<int N>
	const Num<N> nextPrimeOpt35(const BigNum::Num<N>& from) noexcept
	{
		BigNum::Num<N> prime = from | 1;
		
		constexpr unsigned char steps[3][5] = {{2, 2, 2, 4, 2}, {4, 6, 4, 4, 4}, {2, 2, 2, 6, 2}};

		unsigned char rest3 = prime % 3, rest5 = prime % 5;
		if (rest3 == 0 || rest5 == 0)
		{
			unsigned char step = steps[rest3][rest5];
			rest3 = (rest3 + step) % 3;
			rest5 = (rest5 + step) % 5;
			prime = prime + step;
		}

		while (!(simpleTest(prime, 2) && millerRabinTest(prime)))
		{
			assert(prime % 3 != 0);
			assert(prime % 5 != 0);
			unsigned char step = steps[rest3][rest5];
			rest3 = (rest3 + step) % 3;
			rest5 = (rest5 + step) % 5;
			prime = prime + step;
		}

		return prime;
	}
	
	template<int N>
	const Num<N> nextPrimeOpt357(const BigNum::Num<N>& from) noexcept
	{
		BigNum::Num<N> prime = from | 1;
		
		constexpr unsigned char steps[3][5][7] = {
			{{2, 2, 2, 2, 2, 4, 2}, {2, 2, 2, 2, 2, 8, 2}, {2, 2, 2, 2, 2, 4, 2}, {4, 4, 4, 8, 4, 4, 4}, {2, 2, 2, 2, 2, 4, 2}},
			{{4, 4, 4, 6, 4, 4, 4}, {6, 10, 6, 6, 6, 6, 6}, {4, 4, 4, 6, 4, 4, 4}, {4, 4, 4, 6, 4, 4, 4}, {4, 4, 4, 10, 4, 4, 4}},
			{{2, 2, 2, 2, 2, 6, 2}, {2, 2, 2, 2, 2, 6, 2}, {2, 2, 2, 2, 2, 6, 2}, {6, 8, 6, 6, 6, 6, 6}, {2, 2, 2, 2, 2, 8, 2}}
		};

		unsigned rest3 = prime % 3, rest5 = prime % 5, rest7 = prime % 7;
		
		if (rest3 == 0 || rest5 == 0 || rest7 == 0)
		{
			const unsigned step = steps[rest3][rest5][rest7];
			rest3 = (rest3 + step) % 3;
			rest5 = (rest5 + step) % 5;
			rest7 = (rest7 + step) % 7;
			prime = prime + step;
		}

		while (!(simpleTest(prime, 3) && millerRabinTest(prime)))
		{
			assert(prime % 3 != 0);
			assert(prime % 5 != 0);
			assert(prime % 7 != 0);
			const unsigned step = steps[rest3][rest5][rest7];
			rest3 = (rest3 + step) % 3;
			rest5 = (rest5 + step) % 5;
			rest7 = (rest7 + step) % 7;
			prime = prime + step;
		}

		return prime;
	}
	
	template<int N>
	const Num<N> nextPrimeOpt35711(const BigNum::Num<N>& from) noexcept
	{
		BigNum::Num<N> prime = from | 1;
		
		constexpr unsigned char steps[3][5][7][11] = {
			{{{2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},  {4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2}},  {  {2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},  {8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2}},  {  {2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},  {4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2}},  {  {4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},  {8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},  {4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4}},  {  {2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},  {4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2}}},
			{{{4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},  {6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},  {4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4}},  {  {6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},  {10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},  {6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},  {6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},  {6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},  {6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},  {6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6}},  {  {4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},  {6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},  {4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4}},  {  {4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},  {6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},  {4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4}},  {  {4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},  {10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},  {4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},  {4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4}}},
			{{{2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},  {6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2}},  {  {2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},  {6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2}},  {  {2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},  {6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2}},  {  {6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},  {8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8},  {6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},  {6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},  {6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},  {6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},  {6, 6, 6, 6, 6, 14, 6, 6, 6, 6, 6}},  {  {2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},  {8, 8, 8, 12, 8, 8, 8, 8, 8, 8, 8},  {2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2}}}
		};

		unsigned rest3 = prime % 3, rest5 = prime % 5, rest7 = prime % 7, rest11 = prime % 11;
		
		if (rest3 == 0 || rest5 == 0 || rest7 == 0 || rest11 == 0)
		{
			const unsigned step = steps[rest3][rest5][rest7][rest11];
			rest3 = (rest3 + step) % 3;
			rest5 = (rest5 + step) % 5;
			rest7 = (rest7 + step) % 7;
			rest11 = (rest11 + step) % 11;
			prime = prime + step;
		}

		while (!(simpleTest(prime, 4) && millerRabinTest(prime)))
		{
			assert(prime % 3 != 0);
			assert(prime % 5 != 0);
			assert(prime % 7 != 0);
			assert(prime % 11 != 0);
			const unsigned step = steps[rest3][rest5][rest7][rest11];
			rest3 = (rest3 + step) % 3;
			rest5 = (rest5 + step) % 5;
			rest7 = (rest7 + step) % 7;
			rest11 = (rest11 + step) % 11;
			prime = prime + step;
		}

		return prime;
	}
	
	template<int N>
	const Num<N> nextPrimeOpt3571113(const BigNum::Num<N>& from) noexcept
	{
		BigNum::Num<N> prime = from | 1;
		
		constexpr unsigned char steps[3][5][7][11][13] = {
		{{{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 14, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{8, 8, 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 14, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{8, 8, 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 14, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{4, 4, 4, 4, 4, 4, 4, 4, 4, 14, 4, 4, 4},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2}}
		},{{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 16, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{8, 8, 8, 8, 8, 16, 8, 8, 8, 8, 8, 8, 8},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2}},
		{{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 20, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{10, 10, 10, 20, 10, 10, 10, 10, 10, 10, 10, 10, 10},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 16, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{10, 10, 10, 16, 10, 10, 10, 10, 10, 10, 10, 10, 10},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2}}
		},{{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 14, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{10, 10, 10, 14, 10, 10, 10, 10, 10, 10, 10, 10, 10},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 14, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{4, 4, 4, 4, 4, 4, 4, 4, 4, 14, 4, 4, 4},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 14, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{10, 10, 10, 14, 10, 10, 10, 10, 10, 10, 10, 10, 10},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2}}
		},{{{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4}},
		{{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{10, 10, 10, 14, 10, 10, 10, 10, 10, 10, 10, 10, 10},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 14, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{8, 8, 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 14, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{10, 10, 10, 14, 10, 10, 10, 10, 10, 10, 10, 10, 10},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4}}
		},{{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{8, 8, 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2}}
		}},
		{{{{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{12, 16, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 16, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 16, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{6, 6, 6, 6, 6, 6, 6, 16, 6, 6, 6, 6, 6},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4}},
		{{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{12, 16, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 16, 6, 6, 6, 6, 6}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4}}
		},{{{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6}},
		{{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{12, 16, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{10, 10, 10, 16, 10, 10, 10, 10, 10, 10, 10, 10, 10}},
		{{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 16, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{10, 10, 10, 16, 10, 10, 10, 10, 10, 10, 10, 10, 10},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6}},
		{{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6}},
		{{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{12, 16, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 16, 6, 6, 6, 6, 6}},
		{{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6}},
		{{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6}}
		},{{{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4}},
		{{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4}}
		},{{{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 16, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{10, 10, 10, 16, 10, 10, 10, 10, 10, 10, 10, 10, 10},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4}},
		{{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 16, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{10, 10, 10, 16, 10, 10, 10, 10, 10, 10, 10, 10, 10},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 16, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{6, 6, 6, 6, 6, 6, 6, 16, 6, 6, 6, 6, 6},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4}}
		},{{{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 18, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{10, 10, 10, 18, 10, 10, 10, 10, 10, 10, 10, 10, 10},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4}},
		{{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{12, 22, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{10, 10, 10, 22, 10, 10, 10, 10, 10, 10, 10, 10, 10}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{12, 18, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 18, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4}},
		{{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{10, 10, 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4},{4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4}}
		}},
		{{{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{8, 8, 8, 8, 8, 12, 8, 8, 8, 8, 8, 8, 8},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2}},
		{{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{8, 8, 8, 8, 8, 12, 8, 8, 8, 8, 8, 8, 8},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2}}
		},{{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{8, 8, 8, 8, 8, 12, 8, 8, 8, 8, 8, 8, 8},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2}},
		{{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{8, 8, 8, 8, 8, 12, 8, 8, 8, 8, 8, 8, 8},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2}}
		},{{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{12, 14, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 14, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 14, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{6, 6, 6, 6, 6, 6, 6, 14, 6, 6, 6, 6, 6},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2}},
		{{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{12, 14, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 14, 6, 6, 6, 6, 6}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2},{6, 6, 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2}}
		},{{{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 18, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{8, 8, 8, 8, 8, 18, 8, 8, 8, 8, 8, 8, 8},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6}},
		{{8, 8, 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8},{14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 18},{8, 8, 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 18, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8}},
		{{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 14, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{8, 8, 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6}},
		{{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 14, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{8, 8, 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6}},
		{{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 14, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{8, 8, 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6}},
		{{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 14, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{8, 8, 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6}},
		{{6, 6, 6, 6, 6, 6, 6, 14, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 14, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 14, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 14, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 14, 6, 6, 6, 6, 6},{14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 18},{6, 6, 6, 6, 6, 6, 6, 14, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 14, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 18, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 14, 6, 6, 6, 6, 6},{6, 6, 6, 6, 6, 6, 6, 14, 6, 6, 6, 6, 6}}
		},{{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{8, 8, 8, 8, 8, 12, 8, 8, 8, 8, 8, 8, 8},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{8, 8, 8, 8, 8, 12, 8, 8, 8, 8, 8, 8, 8},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 14, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{8, 8, 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{8, 8, 8, 8, 8, 12, 8, 8, 8, 8, 8, 8, 8},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2},{8, 8, 8, 8, 8, 12, 8, 8, 8, 8, 8, 8, 8},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2}},
		{{8, 8, 8, 8, 8, 12, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 12, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 12, 8, 8, 8, 8, 8, 8, 8},{12, 14, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12},{8, 8, 8, 8, 8, 12, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 12, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 12, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 12, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 12, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 12, 8, 8, 8, 8, 8, 8, 8},{8, 8, 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8}},
		{{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2},{12, 14, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12},{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 14, 2}}
		}}};

		unsigned rest3 = prime % 3, rest5 = prime % 5, rest7 = prime % 7, rest11 = prime % 11, rest13 = prime % 13;
		
		if (rest3 == 0 || rest5 == 0 || rest7 == 0 || rest11 == 0 || rest13 == 0)
		{
			const unsigned step = steps[rest3][rest5][rest7][rest11][rest13];
			rest3 = (rest3 + step) % 3;
			rest5 = (rest5 + step) % 5;
			rest7 = (rest7 + step) % 7;
			rest11 = (rest11 + step) % 11;
			rest13 = (rest13 + step) % 13;
			prime = prime + step;
		}

		while (!(simpleTest(prime, 5) && millerRabinTest(prime)))
		{
			assert(prime % 3 != 0);
			assert(prime % 5 != 0);
			assert(prime % 7 != 0);
			assert(prime % 11 != 0);
			assert(prime % 13 != 0);
			const unsigned step = steps[rest3][rest5][rest7][rest11][rest13];
			rest3 = (rest3 + step) % 3;
			rest5 = (rest5 + step) % 5;
			rest7 = (rest7 + step) % 7;
			rest11 = (rest11 + step) % 11;
			rest13 = (rest13 + step) % 13;
			prime = prime + step;
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
	const Num<N> randPrime(RandType type = RandType::FindNext) noexcept
	{
		static BigNum::Num<N> mask = (BigNum::Num<N>(1) | (BigNum::Num<N>(1) << (N - 1)));
		BigNum::Num<N> num = rand<N>() | mask;
		switch(type)
		{
		case RandType::FindNext:
			return nextPrime(num);
		case RandType::FindNextOpt:
			return nextPrimeOpt35711(num);
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
		template <int N>
		void m10add(BigNum::Num<N>& ret, const char arg) noexcept
		{
			assert(arg >= '0' && arg <= '9' && "Must be number");
			ret = ret * 10 + (arg - '0');
		}

		template <int N>
		void num_traverse(BigNum::Num<N>& ret, const char arg) noexcept
		{
			m10add(ret, arg);
		}

		template <int N, class ... Args>
		void num_traverse(BigNum::Num<N>& ret, const char arg, Args... args) noexcept
		{
			m10add(ret, arg);
			num_traverse(ret, args...);
		}

#define DECLARE_BIGNUM_OPERATOR(size) 	template <char... args> BigNum::Num<size> operator "" _bn ## size() { BigNum::Num<size> ret(0); num_traverse(ret, args...); return ret; }
		DECLARE_BIGNUM_OPERATOR(8192)
		DECLARE_BIGNUM_OPERATOR(4096)
		DECLARE_BIGNUM_OPERATOR(2048)
		DECLARE_BIGNUM_OPERATOR(1024)
		DECLARE_BIGNUM_OPERATOR(512)
		DECLARE_BIGNUM_OPERATOR(256)
		DECLARE_BIGNUM_OPERATOR(128)
		DECLARE_BIGNUM_OPERATOR(64)
#undef DECLARE_BIGNUM_OPERATOR
	}
};

using namespace BigNum::operators;