#define _GNU_SOURCE
#include<gmp.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<inttypes.h>
#include<time.h>
#include<pthread.h>
#include<sched.h>
#include<string.h>
#include<stdbool.h>
#include<sys/resource.h>
#include<sys/types.h>
#include<sys/syscall.h>
#include<unistd.h>
#include<errno.h>

#define MAX_THREADS 2 // using an intel core i5 dual-core processor with Hyperthreading (4 logical processors -- 2 per core)
#define POLICY SCHED_OTHER // Standard RR scheduling algorithm of Linux. Tried RT process, too slow.
#define NICENESS -20 // Global value. The lowest nicest translate to highest priority for static threads
// DIFFERENT POWERS OF 2: 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288
#define M_SIZE 16384 // Most optimal 16384 -- Minimizes run-time
#define EXTRA_RELATIONS 20
#define MAX_FB_SIZE 256 // Initial FB size -- Dynamically resizes as needed
#define BUFF_SZ 1024    // Initial size of kernel -- Dynamically resized if space is not used
#define T 1.85 // Threshold multiplier

// Stores smooth relations over Sieving interval
typedef struct smooth_relation {
	uint64_t Q; // Q(x) = x^2 - N
	uint64_t x;
	mpz_t QL;
	mpz_t xL;
	int* factorization;  // contains factorization of each smooth relation
} smooth_r;

// this structure defines the object that
//  holds a factor base prime
typedef struct factor_base {
	int prime;	// prime in factor base
	int index; 	// index in factor base
	int sol_1; 	// First solution to x^2 = N mod p
	int sol_2; 	// Second solution to x^2 = N mod p
	int residue;    // N mod p
	double log_prime;  // logarithm of prime
} f_base_i;

// exponent and identity matrix
typedef struct GF_matrix {
	mpz_t* exp_matrix; // exponent matrix modulo 2 -- |FB|x(smooth_count)
	mpz_t* id_matrix;  // record the index i 
	mpz_t* kernel;     // kernel of the matrix
	int null;	   // size of null space
	int num_r;	   // number of rows
	int num_c; 	   // number of columns
} GF_matrix_t;

// this structure defines the object that controls 
//  various details of the quadratic sieve
typedef struct qs_controller {
	mpz_t N;          	// number to be factored
	f_base_i* f_base; 	// pointer to factor base 
	f_base_i* fb_max; 	// pointer to max prime in factor base 
	int fb_len;	  	// length of factor base
	uint64_t smooth_bound;  // Smoothness bound (B)
	uint32_t M; 		// Sieving interval [-M,M]
	mpz_t sqrtN;		// Square root of N (Used for Sieving)
	smooth_r* s_rel;  	// Smooth relations discovered via Sieveing
	GF_matrix_t matrix; 
	double SIEVE_THRESH; 	// Threshold for sieving
	int smooth_count;
	mpz_t f1, f2;  		// non-trivial divisors of N
} qs_controller_t;

// Implement AKS primality test algorithm to test for primality before sieving
//clmake a Runs in polynomial time (Not yet)
int AKS_algo(mpz_t N);

// initializer functions of controller object / and factor base
void qs_init(qs_controller_t* qs_obj, mpz_t N);
void set_factor_base(qs_controller_t* qs_obj);
void sieve_of_eratosthenes(mpz_t* pbv, int BOUND);

// implemented in num_theor.c
// returns 1 if N is a quadratic residue modulp prime, 0 if prime | N, -1 otherwise
int legendre_symbol(uint32_t N, uint32_t prime);
int f_mod_exp(uint32_t b, uint32_t exp, uint32_t p); // fast modular exponentiation
int sqrts_mod_p(uint32_t p, uint32_t a); // square roots modulo a prime p (Tonelli's Method)
void fact_max_pow(uint32_t* f, uint32_t* e, int p);

// implemented in Sieve.c 
void Sieve(qs_controller_t* qs_obj);  // Sieving portion of Quadratic Sieve
void* Sieve_Interval(void* st_array); // Each thread's work on the sieving interval (4 Threads -- Flexible)
//int trial_divide(uint64_t x); // trial divide each x that meets the Sieve threshold to see if it is smooth

// In gauss.c, row reduce matrix to find squares
void matrix_init(qs_controller_t* qs_obj);
void GAUSS(GF_matrix_t* matrix); 
void Linear_Dep(qs_controller_t* qs_obj);

