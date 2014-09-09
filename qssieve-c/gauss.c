
#include "init.h"

void print_matrix(GF_matrix_t* matrix) {

	int row, col;
	
	/*row = 0;

	printf("COL: ");
	for (; row < matrix->num_r; row++) {
		printf("%d ", row);
	}

	printf("\n");*/
	
	printf("Exponent Matrix: \n");
	row = 0;
	for (; row < matrix->num_r; row++) {
	
		col = 0;
		for (; col < matrix->num_c; col++) {
			printf(" %d ", mpz_tstbit(matrix->exp_matrix[row], col));
		}

		printf("\n");
	}
	
	printf("\nIdentity Matrix: \n");

	row = 0;
	for (; row < matrix->num_r; row++) {
	
		col = 0;
		for (; col < matrix->num_r; col++) {
			printf(" %d ", mpz_tstbit(matrix->id_matrix[row], col));
		}

		printf("\n");
	}
	
	printf("\n");

}

void matrix_init(qs_controller_t* qs_obj) {
	
	int row, col, num_row, num_col;

	num_row = qs_obj->fb_len + EXTRA_RELATIONS;
	num_col = qs_obj->fb_len;
	
	qs_obj->matrix.num_r = num_row;
	qs_obj->matrix.num_c = num_col;
	// initialize number of rows of exponent bit-vector 
	qs_obj->matrix.exp_matrix = (mpz_t*) calloc(num_row, sizeof(mpz_t));
	// initialize number of row of identity bit-vector
	qs_obj->matrix.id_matrix = (mpz_t*) calloc(num_row, sizeof(mpz_t));

	row = 0;
	
	for (; row < num_row; row++) {
		mpz_init2(qs_obj->matrix.exp_matrix[row], num_row);
		mpz_init2(qs_obj->matrix.id_matrix[row], num_row);
		
		mpz_setbit(qs_obj->matrix.id_matrix[row], row);

		col = 0;
		for (; col < num_col; col++) {
			if (qs_obj->s_rel[row].factorization[col] % 2 == 1) {
				mpz_setbit(qs_obj->matrix.exp_matrix[row], col);
			}
		}
	}

	print_matrix(&qs_obj->matrix);
}


// Gaussian-Jordan Elimination
void GAUSS(GF_matrix_t* matrix) {
	int i, j;

	mpz_t sw1;
	mpz_t sw2;
	
	mpz_init(sw1);
	mpz_init(sw2);

	bool pivot_f;

	i = 0;

	for (; i < matrix->num_r; i++) {
		
		pivot_f = false;		
		
		j = i;
		for (; j < matrix->num_r; j++) {
			if(mpz_tstbit(matrix->exp_matrix[j], i)) {
				
				mpz_set(sw1, matrix->exp_matrix[j]);
				mpz_set(sw2, matrix->id_matrix[j]);

				mpz_set(matrix->exp_matrix[j], matrix->exp_matrix[i]);
				mpz_set(matrix->id_matrix[j], matrix->id_matrix[i]);

				mpz_set(matrix->exp_matrix[i], sw1);
				mpz_set(matrix->id_matrix[i], sw2);

				pivot_f = true;
				break;
			}
		}

		if (pivot_f) {
			
			j = i+1;
			for (; j < matrix->num_r; j++) {
				if (mpz_tstbit(matrix->exp_matrix[j],i)) {
					mpz_xor(matrix->exp_matrix[j], matrix->exp_matrix[i], matrix->exp_matrix[j]);
					mpz_xor(matrix->id_matrix[j], matrix->id_matrix[i], matrix->id_matrix[j]);
				}
			}
		}
		
		//printf("(%02.2f %% complete)\r", ((double) (i+1) / (double) matrix->num_r) * 100);			
	}


	matrix->kernel = (mpz_t*) calloc(BUFF_SZ, sizeof(mpz_t));
	matrix->null = 0;

	i = 0;
	for (; i < matrix->num_r; i++) {
		if (mpz_cmp_ui(matrix->exp_matrix[i], 0) == 0) {
			mpz_init_set(matrix->kernel[matrix->null], matrix->id_matrix[i]);
			matrix->null++;
		}
	}

	matrix->kernel = (mpz_t*) realloc(matrix->kernel, matrix->null*sizeof(mpz_t));
	print_matrix(matrix);
}

void Linear_Dep(qs_controller_t* qs_obj) {

	mpz_t f1, f2, X, Y, temp2;

	mpz_init(f1);
	mpz_init(f2);
	mpz_init(X);
	mpz_init(Y);
	mpz_init(temp2);
	
	int* prime_powers = (int*) calloc(qs_obj->fb_len, sizeof(int));
	
	int i, j ,k;

	for (; i < qs_obj->matrix.null; i++) {
		mpz_set_si(X, 1);
		mpz_set_si(Y, 1);

		j = 0;
		for (; j < qs_obj->matrix.num_r; j++) {
			if (mpz_tstbit(qs_obj->matrix.kernel[i], j)) {
				mpz_mul(X, X, qs_obj->s_rel[j].xL);
				mpz_mod(X, X, qs_obj->N);

				k = 0;
				for (; k < qs_obj->fb_len; k++) {
					prime_powers[k] += qs_obj->s_rel[j].factorization[k];
				}
			}		
		}
		
		if ((prime_powers[0] >> 1) % 2) {
			mpz_set_si(Y, -1);
		}

		j = 1;
		for (; j < qs_obj->fb_len; j++) {
			if (prime_powers[j] > 0) {
				prime_powers[j] >>= 1;
				mpz_ui_pow_ui(temp2, qs_obj->f_base[j].prime, prime_powers[j]);
				mpz_mul(Y, Y, temp2);
				mpz_mod(Y, Y, qs_obj->N);

				prime_powers[j] = 0;				
			}
		}

		mpz_sub(f1, X, Y);
		mpz_gcd(f1, f1, qs_obj->N);

		if (mpz_cmp(f1, qs_obj->N) != 0 && mpz_cmp_ui(f1, 1) != 0) { // non-trivial divisors found (f1, f2)			
			mpz_set(qs_obj->f1, f1);
			mpz_set(qs_obj->f2, f1);
			mpz_div(qs_obj->f2, qs_obj->N, qs_obj->f1);
			break;
		}
	}
}












