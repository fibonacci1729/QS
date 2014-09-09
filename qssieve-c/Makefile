CC = gcc

CFLAGS = -std=gnu99 $(DEBUG) $(WARNFLAGS) $(LIBS)
LIBS = -lgmp -lm -lpthread 

open:
	gedit QuadraticSieve.c init.h factorbase.c Sieve.c num_theor.c gauss.c &

all:
	$(CC) factorbase.c num_theor.c Sieve.c gauss.c QuadraticSieve.c -o qs $(CFLAGS) 
	#$(CC) factorbase.c num_theor.c Sieve_ut.c QuadraticSieve.c -o qs $(CFLAGS) 
clean:
	rm qs

# gcc *.c -o qs -std=gnu99 -lgmp -lm -lpthread
