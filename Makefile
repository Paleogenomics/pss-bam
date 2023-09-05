CC=gcc
#CFLAGS=-Og
CFLAGS=-gdwarf-2 -g
OBJS=fasta-genome-io.o sam-parse.o
FK_OBJS=fasta-genome-io.o sam-parse.o kmer.o
HOB_OBJS=sam-parse.o
GKC=fasta-genome-io.o kmer.o
LDFLAGS=-lz

all: main fragkon hangover-bam genome-kmer-count

main: $(OBJS)
	@echo Building pss-bam...
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) pss-bam.c -o pss-bam
	@echo Done.

fragkon: $(FK_OBJS) fragkon.c
	@echo Building fragkon...
	$(CC) $(CFLAGS) $(FK_OBJS) $(LDFLAGS) fragkon.c -o fragkon

hangover-bam: $(HOB_OBJS) hangover-bam.c
	@echo Building hangover-bam...
	$(CC) $(CFLAGS) $(HOB_OBJS) hangover-bam.c -o hangover-bam

genome-kmer-count: $(GKC) genome-kmer-count.c
	@echo Building genome-kmer-count...
	$(CC) $(CFLAGS) $(GKC) $(LDFLAGS) genome-kmer-count.c -o genome-kmer-count

fasta-genome-io.o: fasta-genome-io.h fasta-genome-io.c
sam-parse.o: sam-parse.h sam-parse.c
kmer.o: kmer.h kmer.c
