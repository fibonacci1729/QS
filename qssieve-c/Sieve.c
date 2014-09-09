
#include"init.h"

typedef struct st_data {
	uint64_t low_b;
	uint64_t up_b;
	int c_id;
	int BLOCK_NUM;
	double LOG_P;
} s_data_t;

// GLOBALS: used for shared address space of threads
qs_controller_t* QS;
smooth_r* smooths;
uint64_t M;
int smooth_count = 0;
uint64_t SQRT;
uint64_t block_size;
uint64_t N_long;
pthread_mutex_t mutex_1 =  PTHREAD_MUTEX_INITIALIZER;

// Each thread sieves x = SQRT - M + i, where i runs from low_b to up_b
// 	1. Compute Q(x) from low_b to up_b
//	2. Check to see if x = t mod p where t^2 = N mod p
//	3. If 2, then add (log p) to sieve_log_p[i] where i is which x is congruent to t modulo p
//	4. Else 1.
// 	5. Do this until we have more relations than primes in factor base
// ADDED: (This approach has been sidelined for now until SIQS is implemented...)
//	Each thread sieves over 32768 values. (works nicely, smaller values work better -- investigate)
//	Within each 32768, the thread sieves in blocks of M/32 (if M = 32768, then block size = 2048)
//	I'm not sure if this is more efficient than just sieving the entire 32768 serially.
//	This is an attempt to circumvent the problem of each thread sieving for too long while
//	desired number of smooth relations may have already been found
//	TRADEOFFS:
//		- Uses a large sieve interval size for each thread while still operating on smaller chunks
//		- Sentinel variables smooth_c is updated after every block; results in less accesses to mutexed variable
//		  smooth_count. Otherwise, smooth_c is only updated for every smooth relation found, thus a block may
//		  sieve until a relation is found devoid of the fact that the desired number of relations has already
//		  been reached.
void* Sieve_Interval(void* st_array) {

	s_data_t* data = (s_data_t*) st_array;
	
	// Sieving
	int prime_index;
	uint64_t i;
	uint64_t max;
	uint64_t x;
	int soln_1;
	int soln_2;
	int prime;
	int smooth_c = 0;
	int complete = QS->fb_len + EXTRA_RELATIONS;
	
	while(smooth_c < (QS->fb_len + EXTRA_RELATIONS)) {

		data->low_b += (2*data->BLOCK_NUM*M);		
		data->up_b =  data->low_b + block_size; 
		data->BLOCK_NUM++;	
	
		i = data->low_b;
		max = data->up_b;

		for (; i <= max; i++) {
			x = abs(SQRT - M + i);
			prime_index = 2;
			for (; prime_index < QS->fb_len; prime_index++) {
				prime = QS->f_base[prime_index].prime;
				soln_1 = QS->f_base[prime_index].sol_1;
				soln_2 = QS->f_base[prime_index].sol_2;

				if ((x % prime) == soln_1) {
					// x congruent to soln1 modulo p => p | Q(x)
					// add log p to  sieve array position i
					data->LOG_P += QS->f_base[prime_index].log_prime;
			
				} else if ((x % prime) == soln_2) {
					// x congruent to soln2 modulo p => p | Q(x)
					// add log p to sieve array position at i
					data->LOG_P += QS->f_base[prime_index].log_prime;
				}
			}

			if (data->LOG_P >= QS->SIEVE_THRESH) { // trial divide to check for smoothness
				// if remaining unfactored portion is less than sqrt of largest prime in FB then 
				//  the candidate is smooth (Large Primes -- not implemented)

				// trial divide to find smooth relations that factor over factor base
				int prime_i = 1;
				int large_p = QS->fb_max[QS->fb_len-1].prime;

				uint64_t Qx = pow(x,2)- N_long;
				uint64_t Qx_copy = Qx;	
	
				int prime;
				int factorization[QS->fb_len];
				mpz_t temp;
				mpz_init(temp);

				int i_p = 0;

				for (; i_p < QS->fb_len; i_p++) {
					factorization[i_p] = 0;
				}

				mpz_set_si(temp, Qx);
			
				for (; prime_i < QS->fb_len; prime_i++) {
	
					prime = QS->f_base[prime_i].prime;
					while ((Qx % prime) == 0) {
						Qx =  Qx/prime;
					 	factorization[prime_i]++;
					}		
				}


				if (Qx == 1) {

					// smooth_count is shared
					pthread_mutex_lock( &mutex_1 );
					smooth_c = smooth_count++;		
					pthread_mutex_unlock( &mutex_1 );

					if (smooth_c < QS->fb_len+EXTRA_RELATIONS) {
						smooths[smooth_c].Q = Qx_copy;
						smooths[smooth_c].x = x;
						mpz_init_set_ui(smooths[smooth_c].xL, smooths[smooth_c].x);
						mpz_init_set_ui(smooths[smooth_c].QL, smooths[smooth_c].Q);
						smooths[smooth_c].factorization = malloc(QS->fb_len*sizeof(int));
			smooths[smooth_c].factorization = memcpy(smooths[smooth_c].factorization, &factorization[0], sizeof(factorization));
			printf("\tSmooth (%d/%d): Q(x): %ld\t[%ld,%ld]\n", smooth_c, QS->fb_len+EXTRA_RELATIONS, smooths[smooth_c].Q, data->low_b, max);
						printf("\t\tX: %ld\n\t\ti: %ld\n\t\tTID: %d\n",	smooths[smooth_c].x, i, data->c_id);
						printf("\t\tBlock#: %d\n\t\tLOG_P: %f\n", data->BLOCK_NUM, data->LOG_P);
					} else { break; }				
				} 
			} 	

			data->LOG_P = 0.0;
		}
	}

	pthread_exit((void*) st_array);
}

void Sieve(qs_controller_t* qs_obj) {


	pthread_t sieve_thread[MAX_THREADS];
	pthread_attr_t attr[MAX_THREADS];
	struct sched_param p_sched[MAX_THREADS+1];
	cpu_set_t cpuset[MAX_THREADS];
	
	int i;
	
	QS = qs_obj;

	N_long = mpz_get_ui(qs_obj->N);
	
	M = qs_obj->M;

	SQRT = mpz_get_ui(qs_obj->sqrtN);
	block_size = (uint64_t) ((2*M) / MAX_THREADS);	

	s_data_t st_array[MAX_THREADS];
	smooths =  calloc((qs_obj->fb_len + EXTRA_RELATIONS), sizeof(smooth_r));
	
	QS->SIEVE_THRESH = 0.5*log(mpz_get_d(QS->N)) + log(M);
	QS->SIEVE_THRESH -= T*log(QS->fb_max->prime);

	printf("Block Size: %ld M: %ld B: %ld |FB|: %d THRESH: %f\n", block_size, M, QS->smooth_bound, QS->fb_len, QS->SIEVE_THRESH);

	int nt = 0;
	int low_b = 0;
	int up_b = block_size;
	void* status; // for pthread_join

	// set bounds for each sieving thread (lower, upper)
	for (; nt < MAX_THREADS - 1; nt++) {
		st_array[nt].c_id = nt;
		st_array[nt].low_b = low_b;
		st_array[nt].up_b = up_b;
		st_array[nt].LOG_P = 0.0;
		st_array[nt].BLOCK_NUM = 0;
		low_b += block_size + 1;
		up_b += block_size + 1;
	}

	st_array[nt].low_b = low_b;
	st_array[nt].up_b = up_b;
	st_array[nt].c_id = nt;
	st_array[nt].LOG_P = 0.0;
	st_array[nt].BLOCK_NUM = 0;

	//printf("--------------------------------\n");


	// Set maximum niceness
	pid_t pid = getpid();

	if (setpriority(PRIO_PROCESS, pid, NICENESS) < 0 ) { 
		fprintf(stderr, "Error setting niceness for process %d: (%s)\n", pid, strerror(errno)); 
		exit(-1); 
	}
	
	int st;
	int t = 0;
	
	for (; t < MAX_THREADS; t++) {
			
		// clear the cpu set
		CPU_ZERO(&cpuset[t]);
		CPU_SET(2*t, &cpuset[t]);  // each thread assigned to unique cpu set (2 threads)

		pthread_attr_init(&attr[t]);
		p_sched[t].sched_priority = sched_get_priority_max(POLICY);

		pthread_attr_setschedparam(&attr[t], &p_sched[t]);
	
		if (pthread_attr_setaffinity_np(&attr[t], sizeof(cpu_set_t), &cpuset[t]) < 0) { 
			printf("Error setting CPU affinity for thread: %d\n", t);	
			exit(-1); 
		} 
		
		st = pthread_create( &sieve_thread[t], &attr[t], Sieve_Interval, (void*) &st_array[t]);
		if (st) { printf("Error: return code from pthread_create() is %d\n", st); exit(-1); }
	}

	// wait for each thread to finish sieving
	t = 0;
	
	int ret;
	for (; t < MAX_THREADS; t++) {
		if (ret = pthread_attr_destroy(&attr[t])) { printf("Error: return code from pthread_attr_destroy() is: %d\n", ret); exit(-1); }
		st = pthread_join(sieve_thread[t], &status);
		if (st) { printf("Error: return code from pthread_join() is: %d\n", st); exit(-1); }
	}

	printf("FB#: %d SMOOTHS FOUND: %d\n", QS->fb_len, smooth_count);

	i = 0;
	uint64_t product;
	int j;

	for (; i < QS->fb_len + EXTRA_RELATIONS; i++) {
		j = 0;
		product = 1;

		gmp_printf("Relation #%d: Q(x): %Zd\tx: %ld\tProduct: ", i,  smooths[i].QL, smooths[i].xL);
		for (; j < QS->fb_len; j++) {
			product *= pow(QS->f_base[j].prime, smooths[i].factorization[j]);
		}
		printf("%ld\n", product);

		j = 0;
		for (; j < QS->fb_len; j++) {
			printf(" (%d,%d) ", QS->f_base[j].prime, smooths[i].factorization[j]);
		}

		printf("\n");
	}

	qs_obj->smooth_count = smooth_count;
	qs_obj->s_rel = &smooths[0];
}
