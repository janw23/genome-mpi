// Author: Jan Wawszczak 418479

#include <mpi.h>
#include <iostream>
#include <vector>
#include <exception>
#include <cassert>
#include "data_source.h"
#include "mpi_utils.h"
#include "mpi_vector.h"
#include <any> // TODO remove

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

template <typename T>
static void snapshot(std::vector<std::any> &dbg, const mpi_vector<T> &vec) {
    dbg.push_back(vec.gather(0));
}
template <typename T>
static void snapshot(std::vector<std::any> &dbg, const std::vector<T> &vec) {
    dbg.push_back(vec);
}
// DEBUG HELPERS END

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

// TODO remove
std::vector<std::size_t> gspf(std::string const &genome, std::vector<std::any> &dbg) {
    // Put k-mers into B.
    // TODO Squeeze k>1 - mers into B here to optimize num iterations.
    // TODO starting with k>1 would optimize runtime, but not mem usage.
    std::cerr << "seq genome " << genome << "\n";
    std::vector<std::size_t> B(genome.begin(), genome.end());

    std::vector<std::size_t> SA(B.size());

    { // Sort B, keeping original indices.
        std::vector<std::pair<char, std::size_t>> tmp(B.size());
        for (std::size_t i = 0; i < B.size(); i++) {
            tmp[i] = std::make_pair(B[i], i);
        }

        auto comp = [](const std::pair<size_t, size_t> &a, const std::pair<size_t, size_t> &b){ return a.first < b.first; };
        std::sort(tmp.begin(), tmp.end(), comp); // TODO no need to compare by second element of pair
        
        for (std::size_t i = 0; i < tmp.size(); i++) {
            B[i] = tmp[i].first;
            SA[i] = tmp[i].second;
        }
    }

    snapshot(dbg, B);
    snapshot(dbg, SA);

    { // Rebucket
        std::size_t g = 0;
        std::size_t prev = B[0];
        B[0] = 0;
        for (std::size_t i = 1; i < B.size(); i++) {
            if (prev != B[i]) {
                g = i;
            }
            prev = B[i];
            B[i] = g;
        }
    }

    snapshot(dbg, B);

    // Iterate algorithm, starting with k chosen for initial k-mers TODO here k=1
    for (std::size_t h = 1;; h *= 2) {
        { // Reorder to string-order. TODO in place?
            auto tmp = B;
            for (std::size_t i = 0; i < B.size(); i++) {
                tmp[SA[i]] = B[i];
            }
            std::swap(tmp, B);
        }

        snapshot(dbg, B);

        // Shift B.
        std::vector<std::size_t> B2(B.size());
        std::copy(B.begin() + h, B.end(), B2.begin()); // copy values shifted by h

        snapshot(dbg, B2);
        break;

        {   // Sort according to B, B2 keeping original indices.
            std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> tmp(B.size());
            for (std::size_t i = 0; i < B.size(); i++) {
                tmp[i] = std::make_tuple(B[i], B2[i], i);
            }
            std::sort(tmp.begin(), tmp.end()); // TODO no need to compare by third element of tuple

            for (std::size_t i = 0; i < B.size(); i++) {
                B[i] = std::get<0>(tmp[i]);
                B2[i] = std::get<1>(tmp[i]);
                SA[i] = std::get<2>(tmp[i]);
            }
        }

        // Rebucket.
        {
            std::size_t g = 0;
            std::size_t prev = B[0];
            B[0] = 0;
            std::size_t singletons = 1;
            for (std::size_t i = 1; i < B.size(); i++) {
                if (prev != B[i] || B2[i-1] != B2[i]) { // TODO this seems correct but is different than paper?
                    g = i;
                    singletons++;
                }
                prev = B[i];
                B[i] = g;
            }
            if(singletons == B.size()) break; // Check algorithm-done condition.
        }
    }

    return SA;
}

// TODO in the end make it work on vecot rof mpi_vectors
static mpi_vector<size_t>
suffixArray(mpi_vector<char> &genome, std::vector<std::any> &dbg) {
    assert(genome.rank != genome.nprocs - 1 || genome[genome.size() - 1] == '$');
    mpi_vector<size_t> B(genome);
    mpi_vector<uint64_t> SA(B.comm, B.size());

    { // Sort
        auto comp = [](const std::pair<size_t, size_t> &a, const std::pair<size_t, size_t> &b){ return a.first < b.first; };
        auto tmp = indexed(B);
        mpi_sort(tmp, comp);
        unpack_indexed(tmp, B, SA);
    }

    snapshot(dbg, B);
    snapshot(dbg, SA);

    { // Rebucket
        std::pair<size_t, size_t> context = {0, 0}; // prevB, maxB
        B.iter_apply(context, [&B](decltype(context) &ctx, size_t idx) {
            auto [prev, max] = ctx;
            if (prev != B[idx]) {
                prev = B[idx];
                B[idx] = B.global_index(idx);
                max = B[idx];
            } else {
                prev = B[idx];
                B[idx] = max;
            }
            ctx = {prev, max};
        });
    }

    snapshot(dbg, B);

    for (size_t h = 1;; h *= 2) {
        B.reorder(SA);
        snapshot(dbg, B);

        auto B2 = B;
        B2.shift_left(h, 0);
        snapshot(dbg, B2);

        break; // TODO remove
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
    
    std::vector<std::any> dbg_seq, dbg_par;
    auto full_genome = genome.gather(0);
    if (genome.rank == 0) std::cerr << "full_genome.size() = " << full_genome.size() << "\n";

    if (genome.rank == 0) gspf(std::string(full_genome.data(), full_genome.size()), dbg_seq);
    // std::cout << "Node[" << genome.rank << "]: initial genome: " << genome << "\n";
    auto B = suffixArray(genome, dbg_par);
    // std::cout << "Node[" << genome.rank << "] post B: " << B << "\n";

    if (genome.rank == 0) {
        bool same = true;
        assert(dbg_seq.size() == dbg_par.size());
        for (size_t i = 0; i < dbg_par.size(); i++) {
            const auto &par = std::any_cast<std::vector<size_t>>(dbg_par[i]);
            const auto &seq = std::any_cast<std::vector<size_t>>(dbg_seq[i]);
            if (par != seq) {
                std::cerr << "difference at " << i << "\n";
                std::cerr << "par: " << par << "\n";
                std::cerr << "seq: " << seq << "\n";
                same = false;
            }
        }
        std::cerr << "par and seq are the same: " << same << "\n";
    }

    MPI_Finalize();
    return 0;
}