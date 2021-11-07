#CC=gcc
CFLAGS=-O2
CFLAGS=-gdwarf-2 -g

fasta-genome-io.o : fasta-genome-io.h fasta-genome-io.c
	echo "Making fasta-genome-io.o..."
	$(CC) $(CFLAGS) fasta-genome-io.c -c -lz -o fasta-genome-io.o

test-fasta-genome : test-fasta-genome.c fasta-genome-io.o
	echo "Making test-fasta-genome..."
	$(CC) $(CFLAGS) fasta-genome-io.o test-fasta-genome.c -lz -o test-fasta-genome

sam-parse.o : sam-parse.h sam-parse.c
	echo "Making sam-parse.o..."
	$(CC) $(CFLAGS) -c sam-parse.c -o sam-parse.o -lm

