CC=gcc
CFLAGS=-Og
OBJS=fasta-genome-io.o sam-parse.o
LDFLAGS=-lz

all: main

main: $(OBJS)
	@echo Building main...
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) pss-bam.c -o pss-bam
	@echo Done.

fasta-genome-io.o: fasta-genome-io.h fasta-genome-io.c
sam-parse.o: sam-parse.h sam-parse.c