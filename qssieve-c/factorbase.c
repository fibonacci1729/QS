
#include "init.h"

/*----Sieve of Eratosthenes----*/
void sieve_of_eratosthenes(mpz_t* pbv, int BOUND) {
	int prime = 2; // start sieving with 2 to eliminate multiples of 2
	int multiple = 0;

	// the locations of non-primes i are now marked in pbv with 0 entry
	while((prime*prime) <= BOUND) {
		
		multiple = prime + prime;
		for (; multiple <= BOUND; multiple += prime) {
			mpz_setbit(*pbv, multiple);
		}
		
		prime += 1;

		while (mpz_tstbit(*pbv, prime)) {
			prime += 1;
		}
	}
}

/*----Construct Factor Base----*/
void set_factor_base(qs_controller_t* qs_obj) {
	
	f_base_i* factor_base;
	int fb_size; // this will depend on the number of quadratic residues mod p found
	int max = MAX_FB_SIZE;     // maximum size of factor base, resize if necessary
	int prime_index = 0;	   // current index in factor base
	int BOUND;
	mpz_t NUM;
	
	mpz_init_set(NUM, qs_obj->N);

	factor_base = (f_base_i*) calloc(max, sizeof(f_base_i)); // allocation of factor base	

	factor_base[0].prime = -1; // keep track of signs of smooth residue candidates of qudratic polynomial
	factor_base[0].index = prime_index++;
	factor_base[1].prime = 2;
	factor_base[1].index = prime_index++;
	factor_base[1].sol_1 = 1;
	factor_base[1].sol_2 = 1;
	factor_base[1].log_prime = log(2);
	

	// Sieve for primes up to smoothness bound
	mpz_t prime_bit_vec; // bitmap to keep track of primality of integers up to smoothness bound
		       	     // 0 = prime, 1 = composite
	mpz_init2(prime_bit_vec, qs_obj->smooth_bound);
	BOUND = qs_obj->smooth_bound;

	sieve_of_eratosthenes(&prime_bit_vec, BOUND);

	// We only want primes in prime_bit_vec that divide x^2 - N, thus we want
	//  to find primes that are moduli to which N is quadratic residue guaranteeing
	//  existence of x.

	int fb_length = 2; // length of factor base
	uint32_t residue;  // N mod p
	int soln;	   // solution to x^2 = N mod p for each prime in factor base
	uint32_t N = mpz_get_ui(qs_obj->N);
	int p_prime;

	for (p_prime = 3; p_prime <= BOUND; p_prime++) {
		if (!mpz_tstbit(prime_bit_vec, p_prime)) { // prime
	
			residue = mpz_fdiv_ui(NUM, p_prime); // N mod p
			if ((legendre_symbol(residue, p_prime)) == 1) { // N quadratic residue modulo p
				if (fb_length >= max) {
					max <<= 1;
					factor_base = realloc(factor_base, max*sizeof(f_base_i));
				}				

				// add index to factor base then solve x^2 = N mod p for solutions
				factor_base[prime_index].index = prime_index;
				factor_base[prime_index].prime = p_prime;	
				factor_base[prime_index].residue = residue;
				
				// solutions to x^2 = N mod p -- sol_1, and p - sol_1
				int soln1 = sqrts_mod_p(p_prime, residue);
				int soln2 = p_prime - (soln1 % p_prime);

				factor_base[prime_index].sol_1 = soln1 % p_prime;
				factor_base[prime_index].sol_2 = soln2;

				factor_base[prime_index].log_prime = log(p_prime);
				prime_index++;
				fb_length++;			
			}
		}
	}	
	
	factor_base = realloc(factor_base, fb_length*sizeof(f_base_i));
	qs_obj->f_base = &factor_base[0];
	qs_obj->fb_len = fb_length;
	qs_obj->fb_max = &factor_base[fb_length - 1];
}

