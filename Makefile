CC=g++
CFLAGS=-Wall -Wextra -std=c++14

genome_index: genome-seq.cpp
	$(CC) $(CFLAGS) -o genome_index genome-seq.cpp

clean:
	rm -f genome_index