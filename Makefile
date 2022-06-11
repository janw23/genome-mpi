SEQCC=g++
PARCC=mpic++
CFLAGS=-Wall -std=c++17 -O2

genome_index: genome.cpp data_source.cpp data_source.h mpi_utils.h mpi_vector.h
	$(PARCC) $(CFLAGS) -o genome_index genome.cpp data_source.cpp

genome_index_seq: genome-seq.cpp
	$(SEQCC) $(CFLAGS) -o genome_index_seq genome-seq.cpp

clean:
	rm -f genome_index_seq genome_index