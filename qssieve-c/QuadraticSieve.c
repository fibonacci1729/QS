
#include "init.h"	

/*----FUNCTION DECLARATIONS----*/
int quadratic_sieve(mpz_t N);
void Print_FB(qs_controller_t* qs_obj);
/*-----------------------------*/

void Print_FB(qs_controller_t* qs_obj) {
	printf("Factor base length: %d\n", qs_obj->fb_len);
	// DEBUG
	int i = 2;
	for (; i < qs_obj->fb_len; i++) {
		printf("Index: %d\n\tPrime: %d\n\t(N mod p): %d\n\tSoln1: %d\n\tSoln2: %d\n",
							 qs_obj->f_base[i].index, 
							 qs_obj->f_base[i].prime,
							 qs_obj->f_base[i].residue,
							 qs_obj->f_base[i].sol_1,
							 qs_obj->f_base[i].sol_2);
	}
}

/*----Initialize Controller Object----*/
void qs_init(qs_controller_t* qs_obj, mpz_t N) {
	mpz_init_set(qs_obj->N, N);
	double B = mpz_get_d(N);

	qs_obj->smooth_bound = (uint64_t) (exp(0.5*sqrt(log(B)*log(log(B)))));
	qs_obj->M = M_SIZE;

	mpz_init(qs_obj->sqrtN);
	mpz_sqrt(qs_obj->sqrtN, N);

	mpz_init(qs_obj->f1);
	mpz_init(qs_obj->f2);
}

int quadratic_sieve(mpz_t N) {
	// initalize controller
	qs_controller_t qs_obj;
	qs_init(&qs_obj, N);

	// Construct factor base
	//	1. Sieve for primes up to smooth bound, B
	//	2. Choose primes such that N is quadratic residue modulo p
	//	3. Solve for each prime in factor base (and record solutions t1, t2) x^2 = N mod p
	printf("Allocating factor base...\n");
	set_factor_base(&qs_obj);
	Print_FB(&qs_obj);

	printf("Sieving...\n");
	Sieve(&qs_obj);	// Sieve for smooth relations
	
	printf("Constructing Matrix...\n");
	matrix_init(&qs_obj);  // construct matrix
	printf("Gaussian Elimination...\n");
	GAUSS(&qs_obj.matrix);	 // row reduce to find squares

	printf("\nComputing gcd...\n");
	Linear_Dep(&qs_obj);

	gmp_printf("Two non-trivial divisors of %Zd: (%Zd, %Zd)\n", qs_obj.N, qs_obj.f1, qs_obj.f2);

	return 0;

}

int main(int argc, char* argv[]) {

	clock_t t = clock();
	mpz_t N;
	mpz_init_set_str(N, argv[1], 10);

	quadratic_sieve(N);	
	
	t = clock() - t;

	// Not actual time. clock() give CPU clicks however, process may be blocking
	printf("\nFinish Time: %f seconds\n", ((float) t)/CLOCKS_PER_SEC);
	return 0;
}

