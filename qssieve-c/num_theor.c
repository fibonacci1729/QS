#include "init.h"

// Factorization of f nvolving maximal power of p, e is maximal exponent
void fact_max_pow(uint32_t* f, uint32_t* e, int p) {
	
	while ((*f % p) == 0) {
		*f >>= 1;
		*e += 1;		
	}
}

// Tonell's method for Square Roots Modulo a Prime p & (n = N mod p)
int sqrts_mod_p(uint32_t p, uint32_t n) {
	uint32_t Q;
	uint32_t S;
	uint32_t R;

	if (p % 4 == 3) { // p congruent to 3 modulo 4
		R = f_mod_exp(n,  (p+1)>>2, p);		
		return R;	
	}

	// factor out highest powers of 2
	Q = p - 1;
	S = 0;
	
	fact_max_pow(&Q, &S, 2);	

	//printf("(p-1): %d --> (%d,%d)\n", p - 1, Q, S);

	// quadratic non-residue modulo p
	uint32_t nr = 2; 
	while(legendre_symbol(nr, p) != -1) {
		nr++;
	}

	uint32_t i = 2;
	uint32_t c = ((uint32_t) (n * pow(nr, 2))) % p;
	uint32_t k = 1;
	uint32_t exp;
	uint32_t r;

	for (; k <= S - 1; k++) {
		exp = pow(2, S-k-1)*Q;
		r = f_mod_exp(c, exp, p);
		if ((r == -1) || ((r-p) == -1)) {
			uint32_t temp = pow(2,k);
			i += temp;
			c *= pow(nr, temp);
			c %= p;
		}
	}	


	uint32_t x1 = (i*Q) >> 1;
	uint32_t x2 = (Q+1) >> 1;
	
	R = f_mod_exp(nr, x1, p) * f_mod_exp(n, x2, p);

	return R;
}


// Fast modular exponentiation
int f_mod_exp(uint32_t b, uint32_t exp, uint32_t p) {
	
	int result = 1;

	while (exp != 0) {
		if (exp % 2 == 1) {
			result*= b;
			result %= p;
		}

		exp >>= 1;
		b *= b;
		b %= p;

	}
	
	return result;
}

// 1 if N quadratic residue modulo p, 0 if p | N, -1 non-quadratic residue modulo p
int legendre_symbol(uint32_t N, uint32_t prime) {
	if ((N % prime) == 0) return 0;
	
	int l_sym;
	int exp = prime;
		exp -= 1;
		exp >>= 1;

	// compute N^(p-1)/2 mod p -- Euler's Criterion
	//  if 1, N is quadratic residue modulo p
	l_sym = f_mod_exp(N, exp, prime);

	if (l_sym == (prime-1)) { return -1; }

	return 1;
}
