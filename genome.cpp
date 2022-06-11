// Author: Jan Wawszczak 418479

#include <mpi.h>
#include <iostream>
#include <vector>
#include <exception>
#include <cassert>
#include "data_source.h"
#include "mpi_utils.h"
#include "mpi_vector.h"

// DEBUG HELPERS
template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec) {
    for (size_t i = 0; i < vec.size(); i++) {
        os << vec[i];
        if (i < vec.size() - 1) os << " ";
    }
    return os;
}
template <typename T>
std::ostream &operator<<(std::ostream &os, mpi_vector<T> &vec) {
    for (size_t i = 0; i < vec.size(); i++) {
        os << vec[i];
        if (i < vec.size() - 1) os << " ";
    }
    return os;
}


static std::tuple<uint64_t, uint64_t, std::string, std::string, std::string>
parseArgs(int argc, char *argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " :n :m :genome_in :queries_in :queries_out" << std::endl;
        MPI_Finalize();
        exit(1);
    }
    
    try {
        uint64_t n = std::stoull(argv[1]);
        uint64_t m = std::stoull(argv[2]);
        std::string genome_in = argv[3];
        std::string queries_in = argv[4];
        std::string genome_out = argv[5];
        return {n, m, genome_in, queries_in, genome_out};
    } catch (std::exception &e) {
        std::cerr << "parsing input args failed, reason: " << e.what() << std::endl;
        MPI_Finalize();
        exit(1);
    }
}

// TODO in the end make it work on vecot rof mpi_vectors
static mpi_vector<size_t>
suffixArray(mpi_vector<char> &genome) {
    assert(genome.rank != genome.nprocs - 1 || genome[genome.size() - 1] == '$');
    mpi_vector<size_t> B(genome);
    mpi_vector<uint64_t> SA(B.comm, B.size());

    { // Sort
        auto comp = [](const std::pair<size_t, size_t> &a, const std::pair<size_t, size_t> &b){ return a.first < b.first; };
        auto tmp = indexed(B);
        mpi_sort(tmp, comp);
        unpack_indexed(tmp, B, SA);
    }

    { // Rebucket

    }

    return B;
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    //MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN); // TODO just for debug

    auto [num_genomes, num_queries, genome_in, queries_in, queries_out] = parseArgs(argc, argv);
    DataSource data_source(genome_in.data());

    mpi_vector<char> genome(MPI_COMM_WORLD, data_source, 0);
    genome.push_back('$');
    
    std::cout << "Node[" << genome.rank << "]: initial genome: " << genome << "\n";

    auto B = suffixArray(genome);
    auto sorted_genome = mpi_vector<char>(B);

    std::cout << "Node[" << genome.rank << "]: sorted genome: " << sorted_genome << "\n";
    std::cout << "Node[" << genome.rank << "] B: " << B << "\n";

    MPI_Finalize();
    return 0;
}