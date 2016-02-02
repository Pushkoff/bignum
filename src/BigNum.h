#pragma once

#include <algorithm>

namespace BigNum
{
	template<int N>
	class Num
	{
	public:
		enum
		{
			Size = (N + (sizeof(unsigned char) * 8) - 1) / (sizeof(unsigned char) * 8),
		};
	private:
		unsigned char data[Size];
	public:
		Num(unsigned int num = 0) noexcept
		{
			for (int i = 0; i < Size; i++)
			{
				data[i] = num & 0xFF;
				num >>= 8;
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

			return *this;
		}

		unsigned char& operator[](int i) noexcept { return data[i]; }
		const unsigned char& operator[](int i) const noexcept { return data[i]; }

		bool bit(int i) const noexcept
		{
			//assert(i >= 0 && i < Size);
			bool ret = false;
			if (i >= 0 && i < Size)
			{
				const int byte = i / 8;
				const int bit = i % 8;
				ret = ((*this)[byte] & (1 << bit)) != 0;
			}
			return ret;
		}
	};

	template<int N>
	const Num<N> operator + (const Num<N>& v1, const Num<N>& v2)
	{
		Num<N> rez;

		unsigned int carry = 0;
		for (int i = 0; i < Num<N>::Size; i++)
		{
			unsigned int sum = carry + v1[i] + v2[i];
			rez[i] = sum & 0xFF;
			carry = sum >> 8;
		}

		return rez;
	}

	template<int N>
	const Num<N> operator + (Num<N> v1, unsigned char v2)
	{
		unsigned int carry = v2;
		for (int i = 0; i < Num<N>::Size && carry != 0; i++)
		{
			unsigned int sum = carry + v1[i];
			v1[i] = sum & 0xFF;
			carry = sum >> 8;
		}

		return v1;
	}

	//template<int N>
	//const Num<N> operator + (Num<N> v1, unsigned int v2)
	//{
	//	unsigned long long int carry = v2;
	//	for (int i = 0; i < Num<N>::Size && carry != 0; i++)
	//	{
	//		unsigned long long int sum = carry + v1[i];
	//		v1[i] = sum & 0xFF;
	//		carry = sum >> 8;
	//	}

	//	return v1;
	//}

	template<int N>
	const Num<N> operator - (const Num<N>& v)
	{
		return Num<N>(0) - v;
	}

	template<int N>
	const Num<N> operator - (const Num<N>& v1, const Num<N>& v2)
	{
		Num<N> rez;

		int carry = 0;
		for (int i = 0; i < Num<N>::Size; i++)
		{
			int sum = v1[i] + carry - v2[i];
			rez[i] = static_cast<unsigned char>(sum & 0xFF);
			carry = sum >> 8;
		}

		return rez;
	}

	template<int N>
	const Num<N> operator - (Num<N> v1, unsigned char v2)
	{
		int carry = -v2;
		for (int i = 0; i < Num<N>::Size && carry != 0; i++)
		{
			int sum = v1[i] + carry;
			v1[i] = static_cast<unsigned char>(sum & 0xFF);
			carry = sum >> 8;
		}

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
		unsigned int poly[Num<N>::Size * 2] = { 0 };

		for (int i = 0; i < Num<N>::Size; i++)
			for (int j = 0; j < Num<N>::Size; j++)
				poly[i + j] += v1[i] * v2[j];

		Num<N> rez;
		unsigned int carry = 0;
		for (int i = 0; i < Num<N>::Size; i++)
		{
			int sum = carry + poly[i];
			rez[i] = static_cast<unsigned char>(sum & 0xFF);
			carry = sum >> 8;
		}
		return rez;
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

	template<int N>
	const bool operator == (const Num<N>& v1, const Num<N>& v2)
	{
		bool ret = true;
		for (int i = 0; (i < Num<N>::Size) && (ret); i++)
		{
			ret = (v1[i] == v2[i]);
		}
		return ret;
	}

	template<int N>
	const bool operator == (const Num<N>& v1, unsigned char v2)
	{
		bool ret = v1[0] == v2;
		for (int i = 1; (i < Num<N>::Size) && (ret); i++)
		{
			ret = (v1[i] == 0);
		}
		return ret;
	}

	template<int N>
	const bool operator != (const Num<N>& v1, const Num<N>& v2)
	{
		return !(v1 == v2);
	}
	
	template<int N>
	const bool operator != (const Num<N>& v1, unsigned char v2)
	{
		return !(v1 == v2);
	}

	template<int N>
	const bool operator > (const Num<N>& v1, const Num<N>& v2)
	{
		for (int i = Num<N>::Size - 1; (i >=0 ); i--)
		{
			if (v1[i] != v2[i])
			{
				return v1[i] > v2[i];
			}
		}

		return false;
	}

	template<int N>
	const bool operator >= (const Num<N>& v1, const Num<N>& v2)
	{
		for (int i = Num<N>::Size - 1; (i >= 0); i--)
		{
			if (v1[i] != v2[i])
			{
				return v1[i] > v2[i];
			}
		}
		return true;
	}

	template<int N>
	const bool operator < (const Num<N>& v1, const Num<N>& v2)
	{
		for (int i = Num<N>::Size - 1; (i >= 0); i--)
		{
			if (v1[i] != v2[i])
			{
				return v1[i] < v2[i];
			}
		}
		return false;
	}

	template<int N>
	const bool operator <= (const Num<N>& v1, const Num<N>& v2)
	{
		for (int i = Num<N>::Size - 1; (i >= 0); i--)
		{
			if (v1[i] != v2[i])
			{
				return v1[i] <= v2[i];
			}
		}
		return true;
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
			const int bytes = count / 8;
			const int bits = count % 8;
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
					ret[i] = static_cast<unsigned char>(val | carry);
					carry = val >> 8;
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
		for (int i = count; i < Num<N>::Size - 1; i++)
		{
			ret[i - count] = v[i];
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
			const int bytes = count / 8;
			const int bits = count % 8;
			if (bits == 0)
			{
				ret = shift_right_bytes(v, bytes);
			}
			else
			{
				unsigned int carry = v[bytes];
				for (int i = 0; i < Num<N>::Size - bytes; i++)
				{
					unsigned int val = (static_cast<unsigned int>(((i + bytes + 1) < Num<N>::Size) ? v[i + bytes + 1] : 0) << 8 );
					val = (val | carry);
					ret[i] = static_cast<unsigned char>(val >> bits);
					carry = val >> (8);
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

	template<int N>
	void div(const Num<N>& n, const Num<N>& d, Num<N>& q, Num<N>& r)
	{
		Num<N + 8> rest;

		const Num<N+8> shiftedD[8] = { Num<N + 8>(d), Num<N + 8>(d) << 1, Num<N + 8>(d) << 2, Num<N + 8>(d) << 3, Num<N + 8>(d) << 4, Num<N + 8>(d) << 5, Num<N + 8>(d) << 6, Num<N + 8>(d) << 7 };

		for (int i = Num<N>::Size - 1; i >= 0; i--)
		{
			rest = shift_left_bytes(rest, 1u);
			rest[0] = n[i];
			unsigned int val = 0;
			if (rest >= shiftedD[0])
			{
				Num<N+8> current;
				for (int j = 7; j >= 0; --j)
				{
					Num<N+8> test = current + shiftedD[j];
					if (test <= rest)
					{
						val += 1 << j;
						current = test;
					}
				}
				rest = rest - current;
			}
			q[i] = static_cast<unsigned char>(val);
		}
		r = rest;
	}

	template<int N>
	void div(const Num<2*N>& n, const Num<N>& d, Num<2*N>& q, Num<N>& r)
	{
		const Num<N + 8> shiftedD[8] = { Num<N + 8>(d), Num<N + 8>(d) << 1, Num<N + 8>(d) << 2, Num<N + 8>(d) << 3, Num<N + 8>(d) << 4, Num<N + 8>(d) << 5, Num<N + 8>(d) << 6, Num<N + 8>(d) << 7 };
		Num<N + 8> rest(0);
		for (int i = Num<2*N>::Size - 1; i >= 0; i--)
		{
			rest = shift_left_bytes(rest, 1u);
			rest[0] = n[i];
			unsigned int val = 0;
			if (rest >= shiftedD[0])
			{
				Num<N+8> current(0);
				for (int j = 7; j >= 0; --j)
				{
					Num<N+8> test = current + shiftedD[j];
					if (test <= rest)
					{
						val += 1 << j;
						current = test;
					}
				}
				rest = rest - current;
			}
			q[i] = static_cast<unsigned char>(val);
		}
		r = rest;
	}

	template<int N>
	void div(const Num<N>& n, unsigned char d, Num<N>& q, unsigned char& r)
	{
		unsigned short rest = 0;
		for (int i = Num<N>::Size - 1; i >= 0; i--)
		{
			rest = rest<<8;
			rest+= n[i];
			q[i] = (unsigned char)(rest / d);
			rest = rest % d;
		}
		r = (unsigned char)rest;
	}

	template<int N>
	void div(const Num<N>& n, unsigned short d, Num<N>& q, unsigned short& r)
	{
		unsigned int rest = 0;
		for (int i = Num<N>::Size - 1; i >= 0; i--)
		{
			rest = rest << 8;
			rest += n[i];
			q[i] = (unsigned char)(rest / d);
			rest = rest % d;
		}
		r = (unsigned short)rest;
	}

	template<int N>
	void div(const Num<N>& n, unsigned int d, Num<N>& q, unsigned int& r)
	{
		unsigned long long int rest = 0;
		for (int i = Num<N>::Size - 1; i >= 0; i--)
		{
			rest = rest << 8;
			rest += n[i];
			q[i] = (unsigned char)(rest / d);
			rest = rest % d;
		}
		r = (unsigned int)rest;
	}

	template<int N>
	const Num<N> operator / (const Num<N>& v1, const Num<N>& v2)
	{
		if (v1 < v2)
			return Num<N>(0);

		Num<N> q, r;
		div(v1, v2, q, r);
		return q;
	}

	template<int N>
	const Num<N> operator / (const Num<N>& v1, const unsigned char v2)
	{
		Num<N> q;
		unsigned char r = 0;
		div(v1, v2, q, r);
		return q;
	}

	template<int N>
	const Num<N> operator / (const Num<N>& v1, const unsigned short v2)
	{
		Num<N> q;
		unsigned short r = 0;
		div(v1, v2, q, r);
		return q;
	}

	template<int N>
	const Num<N> operator / (const Num<N>& v1, const unsigned int v2)
	{
		Num<N> q;
		unsigned int r = 0;
		div(v1, v2, q, r);
		return q;
	}

	template<int N>
	const Num<N> operator % (const Num<N>& v1, const Num<N>& v2)
	{
		if (v1 < v2)
			return v1;

		Num<N> q, r;
		div(v1, v2, q, r);
		return r;
	}

	template<int N>
	const unsigned char operator % (const Num<N>& v1, const unsigned char v2)
	{
		Num<N> q;
		unsigned char r = 0;
		div(v1, v2, q, r);
		return r;
	}

	template<int N>
	const unsigned short operator % (const Num<N>& v1, const unsigned short v2)
	{
		Num<N> q;
		unsigned short r = 0;
		div(v1, v2, q, r);
		return r;
	}

	template<int N>
	const unsigned short operator % (const Num<N>& v1, const unsigned int v2)
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

	template<int N>
	bool operator && (const Num<N>& v1, const Num<N>& v2)
	{
		return v1 != 0 && v2 != 0;
	}

	template<int N>
	bool operator || (const Num<N>& v1, const Num<N>& v2)
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
	struct BarretReducer
	{
		const Num<N>& n;
		Num<N> x;
		
		BarretReducer(const Num<N>& _n) : n(_n)
		{
			Num<2 * N> q;
			Num<N> r;
			div(Num<2 * N>(1) << N, n, q, r);
			x = Num<N>(q);
		}
		
		const Num<N> operator()(Num<N> a)
		{
			Num<N> q = Num<N>(mul2N(a, x) >> N);
			Num<N> m = a - q * n;
			while (m >= n)
				m = m - n;
			return m;
		}
	};

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

		//BarretReducer<2 * N> modRed(mod);

		for (int i = 0; i < realExpSize; ++i)
		{
			for (int bit = 0; bit < 8; bit++)
			{
				if ((exp[i] & (1u << bit)) != 0)
				{
					Num<2 * N> q;
					div(mul2N(result, power), mod, q, result);
					//result = modRed(mul2N(result, power));
				}
				{
					Num<2 * N> q;
					div(mul2N(power, power), mod, q, power);
					//result = modRed(mul2N(power, power));
				}
			}
		}

		return Num<N>(result);
	}
};