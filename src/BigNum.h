#pragma once

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
				data[i] = other.value(i);

			for (int i = std::min<int>(Size, Num<M>::Size); i < Size; i++)
				data[i] = 0;

			return *this;
		}

		unsigned char& operator[](int i) noexcept { return data[i]; }
		const unsigned char& operator[](int i) const noexcept { return data[i]; }
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
	const Num<N> operator + (const Num<N>& v1, unsigned char v2)
	{
		Num<N> rez;

		unsigned int carry = v2;
		for (int i = 0; i < Num<N>::Size && carry != 0; i++)
		{
			unsigned int sum = carry + v1[i];
			rez[i] = sum & 0xFF;
			carry = sum >> 8;
		}

		return rez;
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
	const Num<N> operator - (const Num<N>& v1, unsigned char v2)
	{
		Num<N> rez;

		int carry = -v2;
		for (int i = 0; i < Num<N>::Size && carry != 0; i++)
		{
			int sum = v1[i] + carry;
			rez[i] = static_cast<unsigned char>(sum & 0xFF);
			carry = sum >> 8;
		}

		return rez;
	}

	template<int N>
	const Num<N> operator * (const Num<N>& v1, const Num<N>& v2)
	{
		int poly[Num<N>::Size * 2] = { 0 };

		for (int i = 0; i < Num<N>::Size; i++)
			for (int j = 0; j < Num<N>::Size; j++)
				poly[i + j] += v1[i] * v2[j];

		Num<N> rez;
		int carry = 0;
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
				for (int i = bytes+1; i < Num<N>::Size; i++)
				{
					unsigned int val = static_cast<unsigned int>(v[i]) << 8;
					ret[i - bytes - 1] = static_cast<unsigned char>((val | carry) >> bits);
					carry = val >> 8;
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
		r = Num<N>(0);
		for (int i = Num<N>::Size - 1; i >= 0; i--)
		{
			r = shift_left_bytes(r, 1u);
			r[0] = n[i];
			unsigned int val = 0;
			if (r >= d)
			{
				unsigned int lo = 1u, hi = 255u;
				val = (hi + lo) / 2u;
				for (;;)
				{
					Num<N> test = d * static_cast<unsigned char>(val);
					if (test > r)
					{
						hi = val;	
					}
					else
					{
						Num<N> rest = r - test;
						if (rest < d)
						{
							r = rest;
							break;
						}
						else
						{
							lo = val+1u;
						}
					}
					val = (hi + lo) / 2u;
				}	
			}
			q[i] = static_cast<unsigned char>(val);
		}
	}

	template<int N>
	const Num<N> operator / (const Num<N>& v1, const Num<N>& v2)
	{
		Num<N> q, r;
		div(v1, v2, q, r);
		return q;
	}

	template<int N>
	const Num<N> operator % (const Num<N>& v1, const Num<N>& v2)
	{
		Num<N> q, r;
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
};