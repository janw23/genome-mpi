// Author: Jan Wawszczak 418479

#include <mpi.h>
#include <iostream>
#include <vector>
#include <exception>
#include <cassert>
#include "data_source.h"
#include "mpi_utils.h"

// DEBUG HELPERS
template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec) {
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


static std::vector<size_t>
suffixArray(std::vector<char> genome_chunk, const MPIContext &mpi_context) {
    auto genome_offset = mpi_context.getNodeGenomeOffset(mpi_context.rank);
    if (mpi_context.rank == mpi_context.nprocs - 1) {
        //genome_chunk.push_back('$');
    }
    std::vector<size_t> B(genome_chunk.size());
    std::copy(genome_chunk.begin(), genome_chunk.end(), B.begin());
    std::vector<size_t> SA(B.size());

    { // Sort
        std::vector<std::pair<size_t, size_t>> tmp(B.size());
        for (size_t i = 0; i < tmp.size(); i++) {
            tmp[i].first = B[i];
            tmp[i].second = genome_offset + i;
        }

        auto comp = [](const std::pair<size_t, size_t> &a, const std::pair<size_t, size_t> &b){ return a.first < b.first; };
        mpiSort(tmp, comp, mpi_context);

        for (size_t i = 0; i < tmp.size(); i++) {
            B[i] = tmp[i].first;
            SA[i] = tmp[i].second;
        }
    }

    return SA;


    // Rebucket

    // Loop h=k, 2k, 4k, 8k, ...:
    //  reorder to string order
    //  check done
    //  B2 := shift B by h
    //  (B, B2, SA) := Sort (B, B2, idx)
    //  rebucket(B, B2)
    //  done := check-all-singleton(B)

    // return SA
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    //MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN); // TODO just for debug

    auto [num_genomes, num_queries, genome_in, queries_in, queries_out] = parseArgs(argc, argv);
    DataSource data_source(genome_in.data());
    MPIContext mpi_context(data_source, MPI_COMM_WORLD);

    // TODO just for tests read a chunk of first genome
    std::vector<char> genome(data_source.getNodeGenomeSize(0));
    data_source.getNodeGenomeValues(0, genome.data());
    
    std::cout << "I am node " << mpi_context.rank << " and got initial genome: " << genome.data() << "\n";

    auto SA = suffixArray(genome, mpi_context);

    std::cout << "I am node " << mpi_context.rank << ", SA: " << SA << "\n";

    MPI_Finalize();
    return 0;
}