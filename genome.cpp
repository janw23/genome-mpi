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
suffixArray_OBSOLETE(std::vector<char> genome_chunk, const MPIContext &mpi_context) {
    assert(!mpi_context.rankLast() || genome_chunk[genome_chunk.size() - 1] == '$');
    assert(!mpi_context.rankLast() || mpi_context.getNodeGenomeSize(mpi_context.rank) == genome_chunk.size());

    auto genome_offset = mpi_context.getNodeGenomeOffset(mpi_context.rank);
    std::vector<size_t> B(genome_chunk.size()); // TODO fit more chars in one machine word
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

    { // Rebucket
        // TODO I could add data overlap by default to optimize communication. So that boundaries would not need to be changed.
        // Receive the boundary.
        size_t prevB = 0, maxB = 0;
        if (mpi_context.rank > 0) {
            size_t tmp[2];
            MPI_Recv(&tmp, 2, MPI_UINT64_T, mpi_context.rank - 1, MPI_ANY_TAG, mpi_context.comm, MPI_STATUS_IGNORE);
            prevB = tmp[0];
            maxB = tmp[1];
        }

        // Perform local rebucketing.
        for (size_t i = 0; i < B.size(); i++) {
            if (prevB != B[i]) {
                prevB = B[i];
                B[i] = genome_offset + i;
                maxB = B[i];
            } else {
                prevB = B[i];
                B[i] = maxB;
            }
        }

        // Send the boundary.
        if (!mpi_context.rankLast()) {
            size_t tmp[2];
            tmp[0] = prevB;
            tmp[1] = maxB;
            MPI_Send(&tmp, 2, MPI_UINT64_T, mpi_context.rank+1, 0, mpi_context.comm);
        }

        // if (mpi_context.rankFirst()) MPI_Send(&B[B.size()-1], 1, MPI_UINT64_T, mpi_context.rank+1, 0, mpi_context.comm);
        // else if (mpi_context.rankLast()) MPI_Recv(&prevB, 1, MPI_UINT64_T, mpi_context.rank-1, MPI_ANY_TAG, mpi_context.comm, MPI_STATUS_IGNORE);
        // else MPI_Sendrecv(&B[B.size()-1], 1, MPI_UINT64_T, mpi_context.rank+1, 0, &prevB, 1, MPI_UINT64_T, mpi_context.rank-1, MPI_ANY_TAG, mpi_context.comm, MPI_STATUSES_IGNORE);
        // MPI_Scan(MPI_IN_PLACE, B.data(), B.size(), MPI_UINT64_T, MPI_MAX, mpi_context.comm);
        return B;
    }

    // Loop h=k, 2k, 4k, 8k, ...:
    for (size_t h = 1;; h *= 2) {
        { // Reorder to string order
            // This means sending B[i] to node resp. for SA[i].
            // TODO I need SA[i] -> rank translation 
            // bucket elems according to the target node
            // exchange bucket sizes
            // alltoallv in place
            
        }


        //  check done
        //  B2 := shift B by h
        //  (B, B2, SA) := Sort (B, B2, idx)
        //  rebucket(B, B2)
        //  done := check-all-singleton(B)
    }

    return SA;
}

// TODO in the end make it work on vecot rof mpi_vectors
static mpi_vector<size_t>
suffixArray(mpi_vector<char> &genome) {

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
    if (mpi_context.rankLast()) genome.push_back('$'); // put guard at the end
    
    std::cout << "I am node " << mpi_context.rank << " and got initial genome: " << genome.data() << "\n";

    auto B = suffixArray_OBSOLETE(genome, mpi_context);

    std::cout << "I am node " << mpi_context.rank << ", B: " << B << "\n";

    MPI_Finalize();
    return 0;
}