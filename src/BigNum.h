#pragma once
#include <assert.h>
#include <algorithm>
#include <array>
#include <type_traits>
#include <cstddef>
#include <climits>
#include <stdint.h>

#define COMPILER_UNKNOWN 0
#define COMPILER_MSVC 1
#define COMPILER_GCC 2
#define COMPILER_ICC 4
#define COMPILER_CLANG 8

#define COMPILER COMPILER_UNKNOWN

#ifndef COMPILER
	#if defined(__clang__)
	#define COMPILER COMPILER_CLANG
	#elif defined(__ICC) || defined(__INTEL_COMPILER)
	#define COMPILER COMPILER_ICC
	#elif defined(__GNUC__) || defined(__GNUG__)
	#define COMPILER COMPILER_GCC
	#elif defined(_MSC_VER)
	#define COMPILER COMPILER_MSVC
	#else
	#define COMPILER COMPILER_UNKNOWN
	#endif
#endif

#if COMPILER & COMPILER_MSVC
#include <intrin.h>
#pragma intrinsic(_umul128)
#elif COMPILER & COMPILER_GCC
#include <x86intrin.h>
#include <immintrin.h>
#endif


namespace BigNum
{
	namespace Core
	{
		template<typename T>
		struct bitsize
		{
			static constexpr size_t value = sizeof(T) * CHAR_BIT;
		};

		template<typename T, size_t N>
		void zero(T(&to)[N]) noexcept
		{
			for (size_t i = 0; i < N; i++)
				to[i] = 0;
		}

		template<typename T, size_t N>
		void one(T(&to)[N]) noexcept
		{
			to[0] = 1;
			for (size_t i = 1; i < N; i++)
				to[i] = 0;
		}

		template<typename T>
		void one(T(&to)[1]) noexcept
		{
			to[0] = 1;
		}

		template<typename T, size_t N>
		T get(const T(&from)[N], size_t pos) noexcept
		{
			return (pos < N) ? from[pos] : T(0);
		}

		template<typename T, size_t N, size_t M>
		void copy(T(&to)[N], const T(&from)[M]) noexcept
		{
			for (size_t i = 0; i < N; i++)
				to[i] = get(from, i);
		}

		template<typename T, size_t N>
		void swap(T(&v1)[N], T(&v2)[N]) noexcept
		{
			for (size_t i = 0; i < N; i++)
				std::swap(v1[i], v2[i]);
		}

		template<typename T>
		int cmp(const T* v1begin, const T* v1end, const T* v2begin, const T* v2end) noexcept
		{
			const std::ptrdiff_t v1len = v1end - v1begin;
			const std::ptrdiff_t v2len = v2end - v2begin;
			const std::ptrdiff_t vmax = (v1len > v2len) ? v1len : v2len;

			for (std::ptrdiff_t i = vmax; i-- > 0;)
			{
				const T v1 = (i < v1len) ? v1begin[i] : 0;
				const T v2 = (i < v2len) ? v2begin[i] : 0;
				if (v1 != v2)
					return v1 > v2 ? 1 : -1;
			}
			return 0;
		}

		template<typename T, size_t N, size_t M>
		int cmp(const T(&v1)[N], const T(&v2)[M]) noexcept
		{
			constexpr size_t v1len = N;
			constexpr size_t v2len = M;
			constexpr size_t vmax = (v1len > v2len) ? v1len : v2len;

			for (size_t i = vmax; i-- > 0;)
			{
				const T v1v = get(v1, i);
				const T v2v = get(v2, i);
				if (v1v != v2v)
					return v1v > v2v ? 1 : -1;
			}
			return 0;
		}

		template<typename T, size_t N>
		int cmp(const T(&v1)[N], const T v2) noexcept
		{
			const T v2arr[] = { v2 };
			return cmp(v1, v2arr);
		}

		template<typename T>
		T adc(T& ret, const T a, const T b) noexcept
		{
			ret = a + b;
			return (ret < a) ? 1 : 0;
		}

		inline uint64_t adc(uint64_t& ret, const uint64_t v1, const uint64_t v2, const uint64_t carry) noexcept
		{
			assert(carry == 0 || carry == 1);
#if COMPILER & (COMPILER_MSVC | COMPILER_GCC)
			return uint64_t(_addcarry_u64((unsigned char)carry, (unsigned long long)v1, (unsigned long long)v2, (unsigned long long*)&ret));
#else
			uint64_t retCarry = adc(ret, v1, v2);
			retCarry = adc(ret, ret, carry) + retCarry;
			return retCarry;
#endif
		}

		template <typename T>
		T adc(T& ret, const T v1, const T v2, const T carry, typename std::enable_if<(sizeof(uint64_t) > sizeof(T))>::type* = 0) noexcept
		{
			constexpr size_t kWordSizeBits = bitsize<T>::value;
			constexpr uint64_t kWordMaskBits = (uint64_t(1) << (kWordSizeBits)) - 1;
			static_assert(sizeof(uint64_t) > sizeof(T), "cant detect overflow");

			const uint64_t sum = uint64_t(v1) + v2 + carry;
			ret = sum & kWordMaskBits;
			return T(sum >> kWordSizeBits);
		}

		template <typename T>
		T add(T* rezbegin, T* rezend, const T* v1begin, const T* v1end, const T* v2begin, const T* v2end) noexcept
		{
			const std::ptrdiff_t rezlen = rezend - rezbegin;
			const std::ptrdiff_t v1len = v1end - v1begin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			T carry = 0;
			for (std::ptrdiff_t i = 0; i < rezlen; i++)
			{
				const T v1 = (i < v1len) ? v1begin[i] : 0;
				const T v2 = (i < v2len) ? v2begin[i] : 0;
				carry = adc(rezbegin[i], v1, v2, carry);
			}
			return carry;
		}

		template<typename T, size_t N, size_t M, size_t K>
		T add(T(&rez)[N], const T(&v1)[M], const T(&v2)[K]) noexcept
		{
			T carry = 0;
			for (size_t i = 0; i < N; i++)
			{
				const T v1v = (i < M) ? v1[i] : 0;
				const T v2v = (i < K) ? v2[i] : 0;
				carry = adc(rez[i], v1v, v2v, carry);
			}
			return carry;
		}

		template<typename T, size_t N, size_t K>
		T add(T(&rez)[N], const T(&v2)[K]) noexcept
		{
			return add(rez, rez, v2);
		}

		template<typename T, size_t N>
		T add(T(&rez)[N], T v2) noexcept
		{
			const T v2arr[1] = { v2 };
			return add(rez, rez, v2arr);
		}

		template<typename T>
		T add(T* rezbegin, T* rezend, const T* v2begin, const T* v2end) noexcept
		{
			const std::ptrdiff_t rezlen = rezend - rezbegin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			T carry = 0;
			for (std::ptrdiff_t i = 0; i < rezlen; i++)
			{
				const T v2 = (i < v2len) ? v2begin[i] : 0;
				carry = adc(rezbegin[i], rezbegin[i], v2, carry);
			}
			return carry;
		}

		template<typename T>
		T sbc(T& ret, const T a, const T b) noexcept
		{
			ret = a - b;
			return (ret > a) ? 1 : 0;
		}

		inline uint64_t sbc(uint64_t& ret, const uint64_t v1, const uint64_t v2, const uint64_t carry) noexcept
		{
			assert(carry == 0 || carry == 1);
#if COMPILER == COMPILER_MSVC
			return uint64_t(_subborrow_u64((unsigned char)carry, (unsigned long long)v1, (unsigned long long)v2, (unsigned long long*)&ret));
#else
			uint64_t retCarry = sbc(ret, v1, v2);
			retCarry = sbc(ret, ret, carry) + retCarry;
			return retCarry;
#endif
		}

		template<typename T>
		T sbc(T& ret, const T v1, const T v2, const T carry, typename std::enable_if<(sizeof(uint64_t) > sizeof(T))>::type* = 0) noexcept
		{
			constexpr size_t kWordSizeBits = bitsize<T>::value;
			constexpr uint64_t kWordMaskBits = (1ull << (kWordSizeBits)) - 1ull;
			static_assert(sizeof(uint64_t) > sizeof(T), "cant detect overflow");

			const uint64_t sum = uint64_t(v1) - v2 - carry;
			ret = sum & kWordMaskBits;
			return (sum >> kWordSizeBits) ? 1u : 0u;
		}

		template<typename T>
		T sub(T* rezbegin, T* rezend, const T* v1begin, const T* v1end, const T* v2begin, const T* v2end) noexcept
		{
			const std::ptrdiff_t rezlen = rezend - rezbegin;
			const std::ptrdiff_t v1len = v1end - v1begin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			T carry = 0;
			for (std::ptrdiff_t i = 0; i < rezlen; i++)
			{
				const T v1 = (i < v1len) ? v1begin[i] : 0;
				const T v2 = (i < v2len) ? v2begin[i] : 0;
				carry = sbc(rezbegin[i], v1, v2, carry);
			}
			return carry;
		}

		template<typename T, size_t N, size_t M, size_t K>
		T sub(T(&rez)[N], const T(&v1)[M], const T(&v2)[K]) noexcept
		{
			T carry = 0;
			for (size_t i = 0; i < N; i++)
			{
				const T v1v = (i < M) ? v1[i] : 0;
				const T v2v = (i < K) ? v2[i] : 0;
				carry = sbc(rez[i], v1v, v2v, carry);
			}
			return carry;
		}

		template<typename T, size_t N, size_t K>
		T sub(T(&rez)[N], const T(&v2)[K]) noexcept
		{
			return sub(rez, rez, v2);
		}

		template<typename T, size_t N>
		T sub(T(&rez)[N], const T v2) noexcept
		{
			const T v2arr[] = { v2 };
			return sub(rez, rez, v2arr);
		}

		template<typename T>
		T sub(T* rezbegin, T* rezend, const T* v2begin, const T* v2end) noexcept
		{
			const std::ptrdiff_t rezlen = rezend - rezbegin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			T carry = 0;
			for (std::ptrdiff_t i = 0; i < rezlen; i++)
			{
				const T v2 = (i < v2len) ? v2begin[i] : 0;
				carry = sbc(rezbegin[i], rezbegin[i], v2, carry);
			}
			return carry;
		}

		inline uint64_t mul128(uint64_t& retLo, const uint64_t v1, const uint64_t v2) noexcept
		{
#if COMPILER == COMPILER_MSVC
			uint64_t mh = 0;
			retLo = _umul128(v1, v2, &mh);
			return mh;
#elif COMPILER == COMPILER_GCC
			const uint128_t r = static_cast<uint128_t>(a) * b;
			retLo = uint64_t(r);
			return uint64_t(r >> 64);
#else
			constexpr size_t kWordSizeBits = bitsize<uint64_t>::value;
			constexpr size_t kHalfWordSizeBits = kWordSizeBits / 2;
			constexpr uint64_t kHalfWordMaskBits = (uint64_t(1) << (kHalfWordSizeBits)) - 1;

			const uint64_t v1l = v1 & kHalfWordMaskBits;
			const uint64_t v1h = v1 >> kHalfWordSizeBits;
			const uint64_t v2l = v2 & kHalfWordMaskBits;
			const uint64_t v2h = v2 >> kHalfWordSizeBits;

			const uint64_t z0 = v1l * v2l;
			const uint64_t z1 = v1l * v2h;
			const uint64_t z2 = v1h * v2l;
			const uint64_t z3 = v1h * v2h;

			const uint64_t carry1 = adc(retLo, z0, z1 << kHalfWordSizeBits);
			const uint64_t carry2 = adc(retLo, retLo, z2 << kHalfWordSizeBits);
			uint64_t hi = 0;
			adc(hi, z1 >> kHalfWordSizeBits, z2 >> kHalfWordSizeBits, carry1);
			adc(hi, hi, z3, carry2);
			return hi;
#endif
		}

		inline uint64_t madd(uint64_t& ret, const uint64_t v1, const uint64_t v2, const uint64_t carry) noexcept
		{
			uint64_t ml = 0;
			uint64_t mh = mul128(ml, v1, v2);
			uint64_t mc = adc(ret, ret, ml);
			mc += adc(ret, ret, carry);
			return mh + mc;
		}

		template<typename T>
		T madd(T& ret, const T v1, const T v2, const T carry, typename std::enable_if<(sizeof(unsigned long long int) > sizeof(T))>::type* = 0) noexcept
		{
			uint64_t m = uint64_t(ret) + uint64_t(v1) * uint64_t(v2) + carry;
			ret = T(m);
			return T(m >> bitsize<T>::value);
		}

		template<typename T>
		void mul(T* rezbegin, T* rezend, const T* v1begin, const T* v1end, const T* v2begin, const T* v2end) noexcept
		{
			const std::ptrdiff_t v1len = v1end - v1begin;
			const std::ptrdiff_t v2len = v2end - v2begin;
			const std::ptrdiff_t rezlen = rezend - rezbegin;

			assert((rezlen) >= (v1len)+(v2len));
			assert(std::any_of(rezbegin, rezend, [&](T w) { return w > 0u; }) == false);

			const ptrdiff_t v1last = std::min(rezlen, v1len);
			for (std::ptrdiff_t v1it = 0; v1it < v1last; ++v1it)
			{
				const ptrdiff_t v2last = std::min(v2len, rezlen - v1it);
				T carry = 0;
				for (std::ptrdiff_t v2it = 0; v2it < v2last; ++v2it)
					carry = madd(rezbegin[v1it + v2it], v1begin[v1it], v2begin[v2it], carry);
				rezbegin[v1it + v2last] = carry;
			}
		}

		template<typename T>
		void mul(T(&rez)[1], const T(&v1)[1], const T(&v2)[1]) noexcept
		{
			rez[0] = T(v1[0] * v2[0]);
		}

		template<typename T, size_t N, size_t M, size_t K>
		void mul(T(&rez)[N], const T(&v1)[M], const T(&v2)[K]) noexcept
		{
			assert(std::any_of(std::begin(rez), std::end(rez), [&](T w) { return w != 0; }) == false);

			if constexpr (N >= M + K)
			{
				for (size_t v1it = 0; v1it < M; ++v1it)
				{
					T carry = 0;
					for (size_t v2it = 0; v2it < K; ++v2it)
						carry = madd(rez[v1it + v2it], v1[v1it], v2[v2it], carry);
					rez[v1it + K] = carry;
				}
			}
			else
			{
				T ret[M + K] = { 0 };
				mul(ret, v1, v2);
				copy(rez, ret);
				/* const int v1last = std::min(N, M);
				for (int v1it = 0; v1it < v1last; ++v1it)
				{
					const int v2last = std::min(K, N - v1it);
					T carry = 0;
					for (int v2it = 0; v2it < v2last; ++v2it)
						carry = madd(rez[v1it + v2it], v1[v1it], v2[v2it], carry);
					rez[v1it + v2last] = carry;
				}*/
			}
		}

		template<typename T, size_t N, size_t M>
		void shl(T(&rez)[N], const T(&v)[M], size_t count) noexcept
		{
			// assert(count < N*kWordSizeBits);
			const size_t words = count / bitsize<T>::value;
			const size_t bits = count % bitsize<T>::value;
			if (bits == 0)
			{
				for (size_t i = 0; i < N; i++)
					rez[i] = get(v, i - words); // unsigned - unsigned -> unsigned 
			}
			else
			{
				for (size_t i = 0; i < N; i++)
					rez[i] = (get(v, i - words) << (bits)) | (get(v, i - words - 1) >> (bitsize<T>::value - bits));
			}
		}

		template<typename T>
		void shl(T* rezbegin, T* rezend, const T* vbegin, const T* vend, size_t count) noexcept
		{
			const size_t vcount = vend - vbegin;
			const size_t words = count / bitsize<T>::value;
			const size_t bits = count % bitsize<T>::value;
			if (bits == 0)
			{
				for (T* rezIt = rezbegin; rezIt != rezend; ++rezIt)
				{
					const size_t offset = (rezIt - rezbegin) - words;
					*rezIt = (offset >= 0 && offset < vcount) ? *(vbegin + offset) : 0;
				}
			}
			else
			{
				for (T* rezIt = rezbegin; rezIt != rezend; ++rezIt)
				{
					const size_t offset1 = (rezIt - rezbegin) - words;
					const size_t offset2 = offset1 - 1;
					const T w1 = (offset1 >= 0 && offset1 < vcount) ? *(vbegin + offset1) : 0;
					const T w2 = (offset2 >= 0 && offset2 < vcount) ? *(vbegin + offset2) : 0;
					*rezIt = (w1 << bits) | (w2 >> (bitsize<T>::value - bits));
				}
			}
		}

		template<typename T, size_t N>
		void shl(T(&rez)[N], size_t count) noexcept
		{
			if (count == 0) return;

			const size_t words = count / bitsize<T>::value;
			const size_t bits = count % bitsize<T>::value;
			if (bits == 0)
			{
				for (size_t i = N; i --> 0;)
					rez[i] = get(rez, i - words);
			}
			else
			{
				T carry = get(rez, N - words - 1) << (bits);
				for (size_t i = N; i --> 0;)
				{
					T v = get(rez, i - words - 1);
					rez[i] = (v >> (bitsize<T>::value - bits)) | carry;
					carry = v << (bits);
				}
			}
		}

		template<typename T, size_t N, size_t M>
		void shr(T(&rez)[N], const T(&v)[M], size_t count) noexcept
		{
			// assert(count < N*kWordSizeBits);
			const size_t words = count / bitsize<T>::value;
			const size_t bits = count % bitsize<T>::value;
			if (bits == 0)
			{
				for (size_t i = 0; i < N; i++)
					rez[i] = get(v, i + words);
			}
			else
			{
				for (size_t i = 0; i < N; i++)
					rez[i] = (get(v, i + words + 1) << (bitsize<T>::value - bits)) | (get(v, i + words) >> bits);
			}
		}

		template<typename T>
		void shr(T* rezbegin, T* rezend, const T* vbegin, const T* vend, size_t count) noexcept
		{
			const size_t vcount = vend - vbegin;
			const size_t words = count / bitsize<T>::value;
			const size_t bits = count % bitsize<T>::value;
			if (bits == 0)
			{
				for (T* rezIt = rezbegin; rezIt != rezend; ++rezIt)
				{
					const size_t offset = (rezIt - rezbegin) + words;
					*rezIt = (offset >= 0 && offset < vcount) ? *(vbegin + offset) : 0;
				}
			}
			else
			{
				for (T* rezIt = rezbegin; rezIt != rezend; ++rezIt)
				{
					const size_t offset1 = (rezIt - rezbegin) + words;
					const size_t offset2 = offset1 + 1;
					const T w1 = (offset1 >= 0 && offset1 < vcount) ? *(vbegin + offset1) : 0;
					const T w2 = (offset2 >= 0 && offset2 < vcount) ? *(vbegin + offset2) : 0;
					*rezIt = (w1 >> bits) | (w2 << (bitsize<T>::value - bits));
				}
			}
		}

		template<typename T, size_t N>
		void shr(T(&rez)[N], size_t count) noexcept
		{
			if (count == 0) return;

			const size_t words = count / bitsize<T>::value;
			const size_t bits = count % bitsize<T>::value;
			if (bits == 0)
			{
				for (size_t i = 0; i< N; ++i)
					rez[i] = get(rez, i + words);
			}
			else
			{
				T carry = get(rez, words) >> (bits);
				for (size_t i = 0; i < N; ++i)
				{
					T v = get(rez, i + words + 1);
					rez[i] = (v << (bitsize<T>::value - bits)) | carry;
					carry = v >> bits ;
				}
			}
		}

		template<typename T, size_t N, size_t M>
		void div(T(&q)[N], T(&n)[N], const T(&d)[M]) noexcept
		{
			if (cmp(n, d) < 0)
			{
				zero(q);
				return;
			}

			T* rbeg = std::end(n);
			T* rend = std::end(n);

			for (size_t i = N; i-- > 0;)
			{
				rbeg--;
				T val = 0;
				if (Core::cmp(rbeg, rend, std::begin(d), std::end(d)) >= 0)
				{
					for (size_t j = bitsize<T>::value; j-- > 0;)
					{
						T shiftedD[M + 1] = { 0 };
						shl(shiftedD, d, j);
						if (Core::cmp(rbeg, rend, std::begin(shiftedD), std::end(shiftedD)) >= 0)
						{
							Core::sub(rbeg, rend, std::begin(shiftedD), std::end(shiftedD));
							val += T(1) << j;
						}
					}
				}
				q[i] = val;
			}
		}

		template<typename T, size_t N>
		void div(T(&q)[N], T(&n)[N], const T d, typename std::enable_if<(sizeof(uint64_t) > sizeof(T))>::type* = 0) noexcept
		{
			static_assert(sizeof(uint64_t) > sizeof(T), "T have to be smaller size than unsigned long long");

			uint64_t rest = 0;
			for (size_t i = N; i-- > 0;)
			{
				rest = rest << bitsize<T>::value;
				rest += n[i];
				n[i] = 0;
				q[i] = T(rest / d);
				rest = rest % d;
			}
			n[0] = (T)rest;
		}

		template<typename T, size_t N>
		void div(T(&q)[N], T(&n)[N], const T d, typename std::enable_if<(sizeof(uint64_t) == sizeof(T))>::type* = 0) noexcept
		{
			const T converted_d[] = { d };
			div(q, n, converted_d);
		}

		template<typename T, size_t N, size_t M>
		void mod(T(&n)[N], const T(&d)[M]) noexcept
		{
			T* rbeg = std::end(n);
			T* rend = std::end(n);

			for (size_t i = N; i-- > 0;)
			{
				rbeg--;
				if (Core::cmp(rbeg, rend, std::begin(d), std::end(d)) >= 0)
				{
					for (size_t j = bitsize<T>::value; j-- > 0;)
					{
						T shiftedD[M + 1] = { 0 };
						shl(shiftedD, d, j);
						if (Core::cmp(rbeg, rend, std::cbegin(shiftedD), std::cend(shiftedD)) >= 0)
							Core::sub(rbeg, rend, std::cbegin(shiftedD), std::cend(shiftedD));
					}
				}
			}
		}

		template<typename T, size_t N>
		T mod(const T(&n)[N], const T d, typename std::enable_if<(sizeof(uint64_t) > sizeof(T))>::type* = 0) noexcept
		{
			static_assert(sizeof(uint64_t) > sizeof(T), "T have to be smaller size than unsigned long long");
			uint64_t rest = 0;
			for (size_t i = N; i-- > 0;)
				rest = (rest << bitsize<T>::value + n[i]) % d;
			return (T)rest;
		}

		template<typename T, typename R, size_t N>
		R mod(const T(&n)[N], const R d, typename std::enable_if<(sizeof(T) > sizeof(R))>::type* = 0) noexcept
		{
			static_assert(bitsize<T>::value % bitsize<R>::value == 0);

			constexpr size_t kRestSizeBits = bitsize<T>::value / 2;
			constexpr T kRestMaskBits = (T(1) << (kRestSizeBits)) - 1;

			T rest = 0;
			for (size_t i = N; i-- > 0;)
			{
				rest = ((rest << kRestSizeBits) + (n[i] >> kRestSizeBits)) % d;
				rest = ((rest << kRestSizeBits) + (n[i] & kRestMaskBits)) % d;
			}
			return (R)rest;
		}

		template<typename T, size_t N>
		T mod(const T(&n)[N], const T d, typename std::enable_if<(sizeof(uint64_t) == sizeof(T))>::type* = 0) noexcept
		{
			const T tmp_d[] = { d };
			T tmp_n[N];
			copy(tmp_n, n);
			mod(tmp_n, tmp_d);
			return tmp_n[0];
		}

		template<typename T, size_t N>
		size_t zeroBitsCount(const T(&v)[N]) noexcept
		{
			size_t ret = 0;
			size_t i = 0;
			for (; i < N && v[i] == 0; ++i)
				ret += bitsize<T>::value;
			
			if (i < N)
			{
				T tmp = v[i];
				while ((tmp & T(1)) == 0)
				{
					ret++;
					tmp >>= 1;
				}
			}
			return ret;
		}

		template<typename T, size_t N, size_t M>
		void gcd(T(&rez)[N], const T(&v1)[N], const T(&v2)[M], typename std::enable_if<(M <= N)>::type* = 0) noexcept
		{
			size_t v1shift = zeroBitsCount(v1);
			size_t v2shift = zeroBitsCount(v2);
			const size_t rezShift = std::min(v1shift, v2shift);

			T v1t[N], v2t[N];
			shr(v1t, v1, v1shift);
			shr(v2t, v2, v2shift);

			while (int cmpRez = cmp(v1t, v2t))
			{
				if (cmpRez < 0)
					swap(v1t, v2t);
				
				sub(v1t, v2t);
				shr(v1t, zeroBitsCount(v1t));
			}
			shl(rez, v1t, rezShift);
		}

		template<typename T, size_t N>
		void lcm(T(&rez)[2 * N], const T(&v1)[N], const T(&v2)[N]) noexcept
		{
			T m[2 * N] = { 0 };
			//zero(m);
			mul(m, v1, v2);
			T g[N];
			gcd(g, v1, v2);
			div(rez, m, g);
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
			if (t < 0 || t > n) t = t + n;
			return t;
		}

		template<typename T, size_t M, size_t N>
		void modInv(T(&rez)[M], const T(&a)[N], const T(&n)[N]) noexcept
		{
			T t[N] = { 0 }, newt[N];
			one(newt);

			T r[N], newr[N];
			copy(r, n);
			copy(newr, a);

			// newr ≠ 0
			while (cmp(newr, T(0)) != 0)
			{
				// quotient := r div newr
				T quotient[N];
				T tmpr[N];
				copy(tmpr, r);
				div(quotient, tmpr, newr);
				{
					//(t, newt) := (newt, t - quotient * newt) 
					T tmp[N];
					copy(tmp, newt);
					T qn[2 * N] = { 0 };
					mul(qn, quotient, newt);
					sub(newt, t, qn);
					copy(t, tmp);
				}
				{
					//(r, newr) := (newr, r - quotient * newr)
					T tmp[N];
					copy(tmp, newr);
					T qn[2 * N] = { 0 };
					mul(qn, quotient, newr);
					sub(newr, r, qn);
					copy(r, tmp);
				}
			}
			if (cmp(r, T(1)) > 0)
			{
				zero(rez);
				return;
			}
			if (cmp(t, n) > 0)
			{
				add(t, t, n);
			}
			copy(rez, t);
			assert(cmp(t, rez) == 0);
		}

		template<typename T, size_t N, size_t M, size_t K>
		void modExp(T(&result)[N], const T(&base)[N], const T(&exp)[K], const T(&modulo)[M]) noexcept
		{
			one(result);
			if (cmp(modulo, T(1)) == 0)
				return;

			size_t realExpSize = K;
			// skip leading zeros
			while (realExpSize > 0 && exp[realExpSize - 1] == 0)
				realExpSize--;

			for (size_t i = realExpSize; i-- > 0;)
			{
				for (size_t j = (bitsize<T>::value); j-- > 0;)
				{
					T tmp[2 * N] = { 0 }, q[2 * N] = { 0 };
					mul(tmp, result, result);
					div(q, tmp, modulo);
					copy(result, tmp);

					if ((exp[i] >> j) & T(1))
					{
						zero(tmp); zero(q);
						mul(tmp, result, base);
						div(q, tmp, modulo);
						copy(result, tmp);
					}
				}
			}
		}

		template<typename T, size_t N>
		struct MonMul
		{
			static constexpr size_t Len = N;
			T n[Len];
			T ninv[Len];

			explicit MonMul(const T(&_n)[N]) noexcept
			{
				copy(n, _n);
				assert(mod(_n, T(2)) != 0);

				T n2N[2 * Len];
				copy(n2N, n);

				const T val_one[] = { 1 };
				// const T B[2] = { 0, 1 };
				T R[2 * Len] = { 0 };
				shl(R, val_one, N * bitsize<T>::value);

				// const Num<2 * N> r(Num<2 * N>(1) << N);
				// assert(gcd(Num<2 * N>(n), r) == 1);

				T Rinv[2 * Len] = { 0 };
				modInv(Rinv, R, n2N);
				// Num<2*N> rinv = modInv<2 * N>(r, n);
				// assert(((rinv) << N) % Num<2 * N>(n) == 1);

				T Ninv2N[2 * Len] = { 0 }, tmp2N[2 * Len] = { 0 };
				shl(tmp2N, Rinv, N*bitsize<T>::value);

				sub(tmp2N, T(1));
				div(Ninv2N, tmp2N, n);

				copy(ninv, Ninv2N);
				assert(cmp(ninv, Ninv2N) == 0);
				// ninv = (((rinv) << N) - 1) / n;
				// assert(((rinv) << N) - (ninv*n) == 1);
			}

			void In(T(&ret)[N], const T(&x)[N]) const noexcept
			{
				T tmp[2 * N] = { 0 };
				shl(tmp, x, N*bitsize<T>::value);
				mod(tmp, n);
				copy(ret, tmp);
			}

			template<int M>
			void In(T(&ret)[N], const T(&x)[M]) const noexcept
			{
				T tmp[N + M] = { 0 };
				shl(tmp, x, N*bitsize<T>::value);
				mod(tmp, n);
				copy(ret, tmp);
			}

			void Out(T(&ret)[N], const T(&x)[N]) const noexcept
			{
				T val[2 * N];
				copy(val, x);
				REDC(ret, val);
			}

			void REDC(T(&ret)[N], const T(&t)[2 * N]) const noexcept
			{
				T m[Len] = { 0 };
				mul(m, t, ninv);
				//Num<N> m = Num<N>(t) * ninv;

				T mn[2 * Len] = { 0 };
				mul(mn, m, n);
				add(mn, t);
				shr(ret, mn, N*bitsize<T>::value);
				//Num<N> u = (t + m * n) >> N;

				while (cmp(ret, n) >= 0)
					sub(ret, n);
				//if (u >= n)
				//	u = u % n;
				//return Num<N>(u);
			}

			void Mul(T(&ret)[N], const T(&a)[N], const T(&b)[N])
			{
				T tmp[N + N] = { 0 };
				mul(tmp, a, b);
				REDC(ret, tmp);
			}
		};

		template<typename T, size_t N>
		void monModMul(T(&ret)[N], const T(&a)[N], const T(&b)[N], const T(&mod)[N]) noexcept
		{
			MonMul<T, N> monMod(mod);
			T xa[N];
			monMod.In(xa, a);
			monMod.Mul(ret, xa, b);
		}

		template<typename T, size_t N, size_t M, size_t K>
		void monModExp(T(&ret)[N], const T(&base)[K], const T(&exp)[M], const T(&mod)[N]) noexcept
		{
			MonMul<T, N> monMod(mod);
			T tmp[N] = { 0 }, power[N] = { 0 }, val1[N] = { 0 };
			one(val1);

			monMod.In(tmp, val1);
			monMod.In(power, base);

			size_t realExpSize = M;
			while (realExpSize > 0 && exp[realExpSize - 1] == 0)
				realExpSize--;

			for (size_t i = realExpSize; i-- > 0;)
			{
				for (size_t j = bitsize<T>::value; j-- > 0;)
				{
					monMod.Mul(tmp, tmp, tmp);
					if (T(exp[i] >> j) & T(1))
						monMod.Mul(tmp, tmp, power);
				}
			}
			monMod.Out(ret, tmp);
		}

		template<typename T, size_t N, size_t M, size_t K>
		void monModExp2ary(T(&ret)[N], const T(&base)[K], const T(&exp)[M], const T(&mod)[N]) noexcept
		{
			MonMul<T, N> monMod(mod);
			T tmp[N] = { 0 }, power[3][N] = { { 0 },{ 0 },{ 0 } }, val1[N] = { 0 };
			one(val1);

			monMod.In(tmp, val1);
			monMod.In(power[0], base);
			monMod.Mul(power[1], power[0], power[0]);
			monMod.Mul(power[2], power[1], power[0]);

			size_t realExpSize = M;
			while (realExpSize > 0 && exp[realExpSize - 1] == 0)
				realExpSize--;

			for (size_t i = realExpSize; i-- > 0;)
			{
				for (size_t j = bitsize<T>::value; j > 0; j -= 2)
				{
					monMod.Mul(tmp, tmp, tmp);
					monMod.Mul(tmp, tmp, tmp);
					int p = ((exp[i] >> (j - 1)) & T(1)) * 2 + ((exp[i] >> (j - 2)) & T(1));
					if (p > 0)
						monMod.Mul(tmp, tmp, power[p - 1]);
				}
			}
			monMod.Out(ret, tmp);
		}
	}

	using Core::bitsize;
	typedef uint64_t Word;
	constexpr size_t kWordSizeBits = bitsize<Word>::value;

	template<size_t N>
	class Num
	{
	public:
		enum: size_t {
			Size = (N + (kWordSizeBits)-1) / (kWordSizeBits)
		};

		Word data[Size] = { 0 };

		constexpr Num(unsigned long long int num = 0) noexcept
		{
			if constexpr (bitsize<unsigned long long int>::value > bitsize<Word>::value)
			{
				for (size_t i = 0; i < Size; i++)
				{
					data[i] = Word(num);
					num >>= bitsize<Word>::value;
				}
			}
			else
			{
				data[0] = Word(num);
			}
		}

		constexpr Num(std::initializer_list<Word> arg) noexcept
		{
			auto it = std::cbegin(arg), it_end = std::cend(arg);
			for (size_t i = 0; it != it_end && i < Size; ++i)
			{
				data[i] = *it;
				++it;
			}
		}

		template<size_t M>
		Num(const Num<M>& other) noexcept
		{
			*this = other;
		}

		template<size_t M>
		Num(Num<M>&& other) noexcept
		{
			*this = other;
		}

		Num<N>& operator = (const Num<N>& other) noexcept
		{
			if (this != &other)
				Core::copy(this->data, other.data);
			return *this;
		}

		template<size_t M>
		Num<N>& operator = (const Num<M>& other) noexcept
		{
			Core::copy(this->data, other.data);
			return *this;
		}

		Word& operator[](size_t i) noexcept { return data[i]; }
		const Word& operator[](size_t i) const noexcept { return data[i]; }

		bool bit(size_t i) const noexcept
		{
			assert(i >= 0 && i < N);
			bool ret = false;
			if (i >= 0 && i < N)
			{
				const size_t byte = i / kWordSizeBits;
				const size_t bit = i % kWordSizeBits;
				ret = ((*this)[byte] & (Word(1) << bit)) != 0;
			}
			return ret;
		}

		Word* begin() noexcept { return &data[0]; }
		Word* end() noexcept { return &data[Size]; }

		const Word* begin() const noexcept { return &data[0]; }
		const Word* end() const noexcept { return &data[Size]; }
	};

	template<size_t N, size_t M>
	const Num<N> operator + (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		Num<N> ret;
		Core::add(ret.data, v1.data, v2.data);
		return ret;
	}

	template<size_t N, size_t M>
	const Num<N> operator += (Num<N>& ret, const Num<M>& v2) noexcept
	{
		Core::add(ret.data, v2.data);
		return ret;
	}

	template<size_t N>
	const Num<N> operator + (const Num<N>& v1, Word v2) noexcept
	{
		Num<N> ret;
		const Word a2[] = { v2 };
		Core::add(ret.data, v1.data, a2);
		return ret;
	}

	template<size_t N>
	const Num<N> operator += (Num<N>& ret, Word v2) noexcept
	{
		const Word a2[] = { v2 };
		Core::add(ret.data, a2);
		return ret;
	}

	template<size_t N>
	const Num<N> operator - (const Num<N>& v) noexcept
	{
		return Num<N>(0) - v;
	}

	template<size_t N, size_t M>
	const Num<N> operator - (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		Num<N> ret;
		Core::sub(ret.data, v1.data, v2.data);
		return ret;
	}

	template<size_t N, size_t M>
	const Num<N> operator -= (Num<N>& ret, const Num<M>& v2) noexcept
	{
		Core::sub(ret.data, v2.data);
		return ret;
	}

	template<size_t N>
	const Num<N> operator - (const Num<N>& v1, Word v2) noexcept
	{
		Num<N> ret;
		const Word a2[] = { v2 };
		Core::sub(ret.data, v1.data, a2);
		return ret;
	}

	template<size_t N>
	const Num<N> operator -= (Num<N>& ret, Word v2) noexcept
	{
		const Word a2[] = { v2 };
		Core::sub(ret.data, a2);
		return ret;
	}

	template<size_t N, size_t M>
	const Num<N + M> operator * (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		Num<N + M> rez(0);
		Core::mul(rez.data, v1.data, v2.data);
		return rez;
	}

	template<size_t N>
	const Num<N + kWordSizeBits> operator * (const Num<N>& v1, Word v2) noexcept
	{
		Num<N + kWordSizeBits> rez(0);
		const Word a2[] = { v2 };
		Core::mul(rez.data, v1.data, a2);
		return rez;
	}

	template<size_t N, size_t M>
	bool operator == (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(v1.data, v2.data) == 0;
	}

	template<size_t N>
	bool operator == (const Num<N>& v1, Word v2) noexcept
	{
		const Word a2[] = { v2 };
		return Core::cmp(v1.data, a2) == 0;
	}

	template<size_t N, size_t M>
	bool operator != (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return !(v1 == v2);
	}

	template<size_t N>
	bool operator != (const Num<N>& v1, Word v2) noexcept
	{
		return !(v1 == v2);
	}

	template<size_t N, size_t M>
	bool operator > (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(v1.data, v2.data) > 0;
	}

	template<size_t N>
	bool operator > (const Num<N>& v1, Word v2) noexcept
	{
		const Word a2[] = { v2 };
		return Core::cmp(v1.data, a2) > 0;
	}

	template<size_t N, size_t M>
	bool operator >= (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(v1.data, v2.data) >= 0;
	}

	template<size_t N>
	bool operator >= (const Num<N>& v1, Word v2) noexcept
	{
		const Word a2[] = { v2 };
		return Core::cmp(v1.data, a2) >= 0;
	}

	template<size_t N, size_t M>
	bool operator < (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(v1.data, v2.data) < 0;
	}

	template<size_t N>
	bool operator < (const Num<N>& v1, Word v2) noexcept
	{
		const Word a2[] = { v2 };
		return Core::cmp(v1.data, a2) < 0;
	}

	template<size_t N, size_t M>
	bool operator <= (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(v1.data, v2.data) <= 0;
	}

	template<size_t N>
	bool operator <= (const Num<N>& v1, Word v2) noexcept
	{
		const Word a2[] = { v2 };
		return Core::cmp(v1.data, a2) <= 0;
	}

	template<size_t N>
	const Num<N> operator << (const Num<N>& v1, size_t bits) noexcept
	{
		Num<N> ret(0);
		if (bits < N)
			Core::shl(ret.data, v1.data, bits);
		return ret;
	}

	template<size_t N>
	Num<N>& operator <<= (Num<N>& v1, size_t bits) noexcept
	{
		if (bits < N)
			Core::shl(v1.data, bits);
		return v1;
	}

	template<size_t N>
	const Num<N> operator >> (const Num<N>& v1, size_t bits) noexcept
	{
		Num<N> ret(0);
		if (bits < N)
			Core::shr(ret.data, v1.data, bits);
		return ret;
	}

	template<size_t N>
	Num<N>& operator >>= (Num<N>& v1, size_t bits) noexcept
	{
		if (bits < N)
			Core::shr(v1.data, bits);
		return v1;
	}

	template<size_t N, size_t M>
	const Num<N> operator / (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		if (v1 < v2)
			return Num<N>(0);

		Num<N> q, r = v1;
		Core::div(q.data, r.data, v2.data);
		return q;
	}

	template<size_t N>
	const Num<N> operator / (const Num<N>& v1, const Word v2) noexcept
	{
		Num<N> q, r = v1;
		Core::div(q.data, r.data, v2);
		return q;
	}

	template<size_t N, size_t M>
	const Num<M> operator % (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		if (v1 < v2)
			return Num<M>(v1);

		Num<N> r = v1;
		Core::mod(r.data, v2.data);
		return r;
	}

	template<size_t N>
	Word operator % (const Num<N>& v1, const Word v2) noexcept
	{
		return Core::mod(v1.data, v2);
	}

	template<size_t N, size_t M>
	const Num<N> operator & (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		Num<N> rez(0);

		for (size_t i = 0; i < Num<N>::Size; i++)
		{
			rez[i] = v1[i] & ((i < M) ? v2[i] : 0);
		}

		return rez;
	}

	template<size_t N, size_t M>
	const Num<N> operator & (const Num<N>& v1, Num<M>&& v2) noexcept
	{
		Num<N> rez(0);

		for (size_t i = 0; i < Num<N>::Size; i++)
		{
			rez[i] = v1[i] & ((i < M) ? v2[i] : 0);
		}

		return rez;
	}

	template<size_t N>
	const Num<N> operator & (const Num<N>& v1, const Word v2) noexcept
	{
		Num<N> rez(v1[0] & v2);
		return rez;
	}

	template<size_t N>
	const Num<N> operator | (const Num<N>& v1, const Num<N>& v2) noexcept
	{
		Num<N> rez;

		for (size_t i = 0; i < Num<N>::Size; i++)
		{
			rez[i] = v1[i] | v2[i];
		}

		return rez;
	}

	template<size_t N>
	const Num<N> operator | (const Num<N>& v1, const Word v2) noexcept
	{
		Num<N> rez = v1;
		rez[0] = v1[0] | v2;
		return rez;
	}

	template<size_t N>
	const Num<N> operator ~ (const Num<N>& v) noexcept
	{
		Num<N> rez;

		for (size_t i = 0; i < Num<N>::Size; i++)
		{
			rez[i] = ~v[i];
		}

		return rez;
	}

	template<size_t N, size_t M>
	bool operator && (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return v1 != 0 && v2 != 0;
	}

	template<size_t N, size_t M>
	bool operator || (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return v1 != 0 || v2 != 0;
	}

	template<size_t N>
	bool operator !(const Num<N>& v) noexcept
	{
		return v == 0;
	}

	template<size_t N>
	const Num<N> gcd(const Num<N>& v1, const Num<N>& v2) noexcept
	{
		Num<N> ret;
		Core::gcd(ret.data, v1.data, v2.data);
		return ret;
	}

	template<size_t N>
	const Num<N> gcd(const Num<N>& v1, Word v2) noexcept
	{
		Num<N> ret;
		Word v2a[] = { v2 };
		Core::gcd(ret.data, v1.data, v2a);
		return ret;
	}

	template<size_t N>
	const Num<2 * N> lcm(const Num<N>& v1, const Num<N>& v2) noexcept
	{
		Num<2 * N> ret;
		Core::lcm(ret.data, v1.data, v2.data);
		return ret;
	}

	template<size_t N>
	const Num<N> modInv(const Num<N>& a, const Num<N>& n) noexcept
	{
		Num<N> rez;
		Core::modInv(rez.data, a.data, n.data);
		return rez;
	}

	template<typename T>
	const T modInv(const T& a, const T& n) noexcept
	{
		return Core::modInv(a, n);
	}

	template<size_t N>
	const Num<N> modExp(const Num<N>& base, const Num<N>& exp, const Num<N>& modulo) noexcept
	{
		Num<N> result;
		Core::modExp(result.data, base.data, exp.data, modulo.data);
		return result;
	}

	template<size_t N>
	struct MonMul
	{
		const Num<N> n;
		Num<N> ninv;

		explicit MonMul(const Num<N>& _n) noexcept : n(_n)
		{
			const Num<2 * N> r(Num<2 * N>(1) << N);

			assert(gcd(Num<2 * N>(n), r) == 1);

			Num<2 * N> rinv = modInv<2 * N>(r, n);

			assert(((rinv) << N) % Num<2 * N>(n) == 1);

			ninv = (((rinv) << N) - 1) / n;

			assert(((rinv) << N) - (ninv*n) == 1);
		}

		const Num<N> In(const Num<N>& x) const noexcept
		{
			Num<2 * N> ret(x);
			ret = ret << N;

			return ret % n;
		}

		template<size_t M>
		const Num<N> In(const Num<M>& x) const noexcept
		{
			Num<N + M> ret(x);
			ret = ret << N;

			return ret % n;
		}

		const Num<N> Out(const Num<N>& x) const noexcept
		{
			return (*this)(x, Num<N>(1));
		}

		const Num<N> operator()(const Num<2 * N>& t) const noexcept
		{
			Num<N> m = Num<N>(t) * ninv;
			Num<N> u = (t + m * n) >> N;

			if (u >= n)
				u = u % n;
			return Num<N>(u);
		}

		const Num<N> operator()(const Num<N>& a, const Num<N>& b) const noexcept
		{
			Num<2 * N> t = a * b;
			return (*this)(t);
		}
	};

	template<size_t N>
	const Num<N> monModMul(const Num<N>& a, const Num<N>& b, const Num<N>& mod) noexcept
	{
		Num<N> ret;
		Core::monModMul(ret.data, a.data, b.data, mod.data);
		return ret;
	}

	template<size_t N, size_t M, size_t K>
	const Num<N> monModExp(const Num<K>& base, const Num<M>& exp, const Num<N>& mod) noexcept
	{
		Num<N> ret;
		Core::monModExp(ret.data, base.data, exp.data, mod.data);
		return ret;
	}

	template<size_t N, size_t M, size_t K>
	const Num<N> monModExp2ary(const Num<K>& base, const Num<M>& exp, const Num<N>& mod) noexcept
	{
		Num<N> ret;
		Core::monModExp2ary(ret.data, base.data, exp.data, mod.data);
		return ret;
	}

	template<size_t N>
	void CRT_init(const Num<N>& p, const Num<N>& q, const Num<2 * N>& d, Num<N>& dp, Num<N>& dq, Num<N>& qinvp)
	{
		qinvp = modInv(q, p);
		dp = d % (p - 1);
		dq = d % (q - 1);
	}

	template<size_t N>
	const Num<2 * N> CRT(const Num<2 * N>& m, const Num<N>& p, const Num<N>& q, const Num<N>& dp, const Num<N>& dq, const Num<N>& qinvp) noexcept
	{
		Num<N> mp = monModExp2ary(m, dp, p);
		Num<N> mq = monModExp2ary(m, dq, q);

		return Num<2 * N>(mq) + (((mp - mq) * qinvp) % p) * q;
	}

	template<size_t N>
	const Num<N> rand() noexcept
	{
		Num<N> ret;
		for (size_t i = 0; i < Num<N>::Size; ++i)
		{
			ret[i] = (unsigned)(std::rand() % 256u);
		}
		return ret;
	}

	// http://en.literateprograms.org/index.php?title=Special:DownloadCode/Miller-Rabin_primality_test_(C)&oldid=18973

	template<size_t N>
	bool millerRabinPass(const Num<N>& num, const Num<N>& a) noexcept
	{
		Num<N> d = num - 1;

		size_t s = 0u;
		while (d.bit(s) == false)
			s++;

		d = d >> s;
		Num<N> a_to_power = monModExp2ary<N>(a, d, num);

		if (a_to_power == 1)
			return true;

		for (size_t i = 0; i < s - 1; i++)
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
	int millerRabinProbes(size_t N) noexcept
	{
		if (N >= 600) return 2;
		if (N >= 550) return 4;
		if (N >= 500) return 5;
		if (N >= 400) return 6;
		if (N >= 350) return 7;
		if (N >= 300) return 9;
		if (N >= 250) return 12;
		if (N >= 200) return 15;
		if (N >= 150) return 18;
		if (N >= 100) return 27;
		return 40;
	}

	template<size_t N>
	bool millerRabinTest(const Num<N>& num) noexcept
	{
		const int times = millerRabinProbes(N);
		for (int i = 0; i < times; i++)
		{
			Num<N> a = (rand<N>() % (num - 2)) + 2;
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
		4943,   4951,   4957,   4967,   4969,   4973,   4987,   4993,   4999,   5003,
	};

	constexpr size_t primesCount = sizeof(primes) / sizeof(primes[0]);

	constexpr size_t simplechecks(size_t size) noexcept
	{
		return std::min(std::max(size/64, size_t(4)), primesCount);
	}

	template<size_t N>
	bool simpleTest(const Num<N>& num, size_t skip = 0) noexcept
	{
		for (size_t i = skip; i < simplechecks(N); i++)
		{
			assert((num % Word(primes[i]) == 0) == (Core::mod(num.data, primes[i]) == 0));
			assert((num % Word(primes[i]) == 0) == (gcd(num, Word(primes[i])) != 1));
			//if (num % Word(primes[i]) == 0)
			if (Core::mod(num.data, primes[i]) == 0)
			//if (gcd(num, Word(primes[i])) != 1)
				return false;
		}
		return true;
	}

	template<size_t N>
	constexpr Num<N> primesProduct() noexcept
	{
		if constexpr (N >= 2048)
			return Num<N>{ 10408628748038725759ull, 1915073011592289563ull, 8763774421824524631ull, 2011063807001048133ull, 1031707741487195015ull, 10556206784421030708ull, 5144671606306936534ull,
			6265034155426823908ull, 6436778927562491100ull, 2992101046236765559ull, 6886016094704050588ull, 4836091997350331295ull, 11192339607704847528ull, 18185458679113843091ull,
			2237090061021923320ull, 6966528701537520681ull, 13676727934563798704ull, 35440123139362773ull, 1030616290511923130ull, 2598059789510833357ull, 11358731858419258879ull,
			17141514317551676213ull, 15674928412557522658ull, 5788948221837802822ull, 12584695127497503830ull, 6903558455145820541ull, 3240358321303087788ull, 12411608971331568287ull,
			11938157545863749836ull, 9434899357668762554ull, 11386550743579174540ull, 2622686790435282460ull };
		else if constexpr (N >= 1024)
			return Num<N>{ 439843546635370873ull, 5874615937175795343ull, 15114307389686269710ull, 3871889402991114301ull, 7173936838642267339ull, 9722714508587535958ull, 1930269961649009581ull,
			12607511240871844615ull, 10714462625239494289ull, 17992178854007237495ull, 12035886451572866178ull, 2022252832547602063ull, 9058264411918715906ull, 13274342168057743064ull,
			1108643872522681714ull, 200515704069442536ull };
		else if constexpr (N >= 512)
			return Num<N>{ 1003437303243261071ull, 11492846679185067928ull, 7993966689410686436ull, 2539380480667780101ull, 13958398667001434503ull, 540084209558748057ull, 12718512968486482079ull,
			1182944749624425070ull };
		else if constexpr (N >= 256)
			return Num<N>{ 13860834134742107767ull, 11463165765559708970ull, 17664302261596423816ull, 15848267622464664512ull };
		else if constexpr (N >= 128)
			return Num<N>{15507271189205528447ull, 6311747033189848393ull};
		else if constexpr (N >= 64)
			return Num<N>{16294579238595022365ull};
		else
			return Num<N>{ 2 * 3 * 5 * 7 };
	}

	template<size_t N>
	bool gcdTest(const Num<N>& num) noexcept
	{
		//return simpleTest(num);
		return gcd(num, primesProduct<N>()) == 1;
	}

	template<size_t N>
	const Num<N> nextPrime(const Num<N>& from) noexcept
	{
		Num<N> prime = from | 1u;

		while (!(gcdTest(prime) && millerRabinTest(prime)))
			prime = prime + 2u;

		return prime;
	}

	template<size_t N>
	const Num<N> nextPrimeOpt3(const Num<N>& from) noexcept
	{
		Num<N> prime = from | 1u;

		constexpr unsigned char steps[3] = { 2, 4, 2 };

		unsigned char rest3 = (unsigned char)(prime % 3u);
		if (rest3 == 0)
		{
			const unsigned char step = steps[rest3];
			rest3 = (rest3 + step) % 3u;
			prime = prime + step;
		}

		while (!(gcdTest(prime) && millerRabinTest(prime)))
		{
			assert(prime % 3u != 0);
			const unsigned char step = steps[rest3];
			rest3 = (rest3 + step) % 3u;
			prime = prime + step;
		}

		return prime;
	}

	template<size_t N>
	const Num<N> nextPrimeOpt35(const Num<N>& from) noexcept
	{
		Num<N> prime = from | 1u;

		constexpr unsigned char steps[3][5] = { { 2, 2, 2, 4, 2 },{ 4, 6, 4, 4, 4 },{ 2, 2, 2, 6, 2 } };

		unsigned char rest3 = (unsigned char)(prime % 3u), rest5 = (unsigned char)(prime % 5u);
		if (rest3 == 0 || rest5 == 0)
		{
			const unsigned char step = steps[rest3][rest5];
			rest3 = (rest3 + step) % 3u;
			rest5 = (rest5 + step) % 5u;
			prime = prime + step;
		}

		while (!(gcdTest(prime) && millerRabinTest(prime)))
		{
			assert(prime % 3u != 0);
			assert(prime % 5u != 0);
			const unsigned char step = steps[rest3][rest5];
			rest3 = (rest3 + step) % 3u;
			rest5 = (rest5 + step) % 5u;
			prime = prime + step;
		}

		return prime;
	}

	template<size_t N>
	const Num<N> nextPrimeOpt357(const Num<N>& from) noexcept
	{
		Num<N> prime = from | 1;

		constexpr unsigned char steps[3][5][7] = {
			{ { 2, 2, 2, 2, 2, 4, 2 },{ 2, 2, 2, 2, 2, 8, 2 },{ 2, 2, 2, 2, 2, 4, 2 },{ 4, 4, 4, 8, 4, 4, 4 },{ 2, 2, 2, 2, 2, 4, 2 } },
			{ { 4, 4, 4, 6, 4, 4, 4 },{ 6, 10, 6, 6, 6, 6, 6 },{ 4, 4, 4, 6, 4, 4, 4 },{ 4, 4, 4, 6, 4, 4, 4 },{ 4, 4, 4, 10, 4, 4, 4 } },
			{ { 2, 2, 2, 2, 2, 6, 2 },{ 2, 2, 2, 2, 2, 6, 2 },{ 2, 2, 2, 2, 2, 6, 2 },{ 6, 8, 6, 6, 6, 6, 6 },{ 2, 2, 2, 2, 2, 8, 2 } }
		};

		unsigned char rest3 = (unsigned char)(prime % 3u), rest5 = (unsigned char)(prime % 5u), rest7 = (unsigned char)(prime % 7u);

		if (rest3 == 0 || rest5 == 0 || rest7 == 0)
		{
			const unsigned char step = steps[rest3][rest5][rest7];
			rest3 = (rest3 + step) % 3u;
			rest5 = (rest5 + step) % 5u;
			rest7 = (rest7 + step) % 7u;
			prime = prime + step;
		}

		while (!(gcdTest(prime) && millerRabinTest(prime)))
		{
			assert(prime % 3u != 0);
			assert(prime % 5u != 0);
			assert(prime % 7u != 0);
			const unsigned char step = steps[rest3][rest5][rest7];
			rest3 = (rest3 + step) % 3u;
			rest5 = (rest5 + step) % 5u;
			rest7 = (rest7 + step) % 7u;
			prime = prime + step;
		}

		return prime;
	}

	template<size_t N>
	const Num<N> nextPrimeOpt35711(const Num<N>& from) noexcept
	{
		Num<N> prime = from | 1;

		constexpr unsigned char steps[3][5][7][11] = {
			{ { { 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2 },{ 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2 } },{ { 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2 },{ 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2 } },{ { 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 10, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2 },{ 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2 } },{ { 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4 },{ 8, 8, 8, 10, 8, 8, 8, 8, 8, 8, 8 },{ 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4 } },{ { 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2 },{ 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2 } } },
			{ { { 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4 },{ 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6 },{ 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4 } },{ { 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6 },{ 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10 },{ 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6 },{ 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6 },{ 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6 },{ 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6 },{ 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6 } },{ { 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4 },{ 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6 },{ 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4 } },{ { 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4 },{ 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6 },{ 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4 } },{ { 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4 },{ 10, 12, 10, 10, 10, 10, 10, 10, 10, 10, 10 },{ 4, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4 },{ 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4 } } },
			{ { { 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2 },{ 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2 } },{ { 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2 },{ 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2 } },{ { 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2 },{ 6, 6, 6, 6, 6, 12, 6, 6, 6, 6, 6 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 2 } },{ { 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6 },{ 8, 8, 8, 14, 8, 8, 8, 8, 8, 8, 8 },{ 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6 },{ 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6 },{ 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6 },{ 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6 },{ 6, 6, 6, 6, 6, 14, 6, 6, 6, 6, 6 } },{ { 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 2 },{ 8, 8, 8, 12, 8, 8, 8, 8, 8, 8, 8 },{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 2 } } }
		};

		unsigned char rest3 = (unsigned char)(prime % 3u), rest5 = (unsigned char)(prime % 5u), rest7 = (unsigned char)(prime % 7u), rest11 = (unsigned char)(prime % 11u);

		if (rest3 == 0 || rest5 == 0 || rest7 == 0 || rest11 == 0)
		{
			const unsigned char step = steps[rest3][rest5][rest7][rest11];
			rest3 = (rest3 + step) % 3u;
			rest5 = (rest5 + step) % 5u;
			rest7 = (rest7 + step) % 7u;
			rest11 = (rest11 + step) % 11u;
			prime = prime + step;
		}

		while (!(gcdTest(prime) && millerRabinTest(prime)))
		{
			assert(prime % 3u != 0);
			assert(prime % 5u != 0);
			assert(prime % 7u != 0);
			assert(prime % 11u != 0);
			const unsigned char step = steps[rest3][rest5][rest7][rest11];
			rest3 = (rest3 + step) % 3u;
			rest5 = (rest5 + step) % 5u;
			rest7 = (rest7 + step) % 7u;
			rest11 = (rest11 + step) % 11u;
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

	template<size_t N>
	const Num<N> randPrime(RandType type = RandType::FindNextOpt) noexcept
	{
		static Num<N> mask = (Num<N>(1) | (Num<N>(1) << (N - 1)));
		Num<N> num = rand<N>() | mask;
		switch (type)
		{
		case RandType::FindNext:
			return nextPrime(num);
		case RandType::FindNextOpt:
			return nextPrimeOpt35711(num);
		case RandType::Simple:
		default:
			{
				while (!(gcdTest(num) && millerRabinTest(num)))
				{
					num = rand<N>() | mask;
				}
				return num;
			}
		}
	}

	template<size_t N, typename Iter>
	const Num<N> d2i(Iter from_begin, Iter from_end) noexcept
	{
		Num<N> ret(0);
		auto it = from_begin;
		for (size_t i = 0; i < Num<N>::Size; i++)
		{
			Word val = 0;
			for (size_t j = 0; j < sizeof(Word); j++)
			{
				if (it != from_end)
				{
					val |= static_cast<Word>(*it) << (j * 8);
					it++;
				}
			}
			ret[i] = val;
		}
		return ret;
	}

	template<size_t N, typename Iter>
	void i2d(const Num<N>& from, Iter to_begin, Iter to_end) noexcept
	{
		auto it = to_begin;
		for (size_t i = 0; i < Num<N>::Size; i++)
		{
			Word val = from[i];
			for (size_t j = 0; j < sizeof(val); j++)
			{
				if (it != to_end)
				{
					(*it) = (val >> (j * 8)) & 0xFF;
					it++;
				}
			}
		}
	}

	namespace operators
	{
		template <size_t N>
		void m10add(Num<N>& ret, const char arg) noexcept
		{
			assert(arg >= '0' && arg <= '9' && "Must be number");
			ret = ret * 10u + (Word(arg) - '0');
		}

		template <size_t N>
		void num_traverse(Num<N>& ret, const char arg) noexcept
		{
			m10add(ret, arg);
		}

		template <size_t N, class ... Args>
		void num_traverse(Num<N>& ret, const char arg, Args... args) noexcept
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
}

using namespace BigNum::operators;