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
    //dbg.push_back(vec.gather(0));
}
template <typename T>
static void snapshot(std::vector<std::any> &dbg, const std::vector<T> &vec) {
    //dbg.push_back(vec);
}

uint64_t num_occurences_naive(std::string const &genome, std::string const &query) {
    uint64_t num_occurences_naive = 0;
    std::size_t pos = 0;
    while ((pos = genome.find(query, pos)) != std::string::npos) {
        num_occurences_naive++;
        pos++;
    }
    return num_occurences_naive;
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

        {   // Sort according to B, B2 keeping original indices.
            std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> tmp(B.size());
            for (std::size_t i = 0; i < B.size(); i++) {
                tmp[i] = std::make_tuple(B[i], B2[i], i);
            }
            auto comp = [](
                const std::tuple<size_t, size_t, size_t> &a,
                const std::tuple<size_t, size_t, size_t> &b){
                    return 
                        std::tie(std::get<0>(a), std::get<1>(a))
                        <
                        std::tie(std::get<0>(b), std::get<1>(b));
            };
            std::sort(tmp.begin(), tmp.end(), comp); // TODO no need to compare by third element of tuple

            for (std::size_t i = 0; i < B.size(); i++) {
                B[i] = std::get<0>(tmp[i]);
                B2[i] = std::get<1>(tmp[i]);
                SA[i] = std::get<2>(tmp[i]);
            }
        }

        snapshot(dbg, B);
        snapshot(dbg, B2);
        snapshot(dbg, SA);

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

        snapshot(dbg, B);
    }

    return SA;
}

// TODO in the end make it work on vecot rof mpi_vectors
static mpi_vector<size_t>
suffix_array(mpi_vector<char> &genome, std::vector<std::any> &dbg) {
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
        auto B2_local = B2.gather(0);
        B2.shift_left(h, 0);
        B2_local = B2.gather(0);
        snapshot(dbg, B2);

        { // Sort
            auto comp = [](
                const std::tuple<size_t, size_t, size_t> &a,
                const std::tuple<size_t, size_t, size_t> &b){
                    return 
                        std::tie(std::get<0>(a), std::get<1>(a))
                        <
                        std::tie(std::get<0>(b), std::get<1>(b));
            };
            auto tmp = indexed(B, B2);
            mpi_sort(tmp, comp);
            unpack_indexed(tmp, B, B2, SA);
        }

        snapshot(dbg, B);
        snapshot(dbg, B2);
        snapshot(dbg, SA);

        { // Rebucket
            std::tuple<size_t, size_t, size_t, uint64_t> context = {0, 0, 0, 0}; // prevB, prevB2, max, singletons
            B.iter_apply(context, [&B, &B2](decltype(context) &ctx, size_t idx) {
                auto [prevB, prevB2, max, singletons] = ctx;
                if (B.global_index(idx) == 0 || prevB != B[idx] || prevB2 != B2[idx]) {
                    max = B.global_index(idx);
                    singletons++;
                }
                prevB = B[idx];
                prevB2 = B2[idx];
                B[idx] = max;
                ctx = {prevB, prevB2, max, singletons};
            });

            // Broadcast info whether algorithm is done.
            char done = (std::get<3>(context) == B.global_size());
            MPI_Bcast(&done, 1, MPI_CHAR, B.nprocs - 1, B.comm);
            if (done) break;
        }

        snapshot(dbg, B);
    }

    return SA;
}

static std::pair<uint64_t, bool> find_lower(const mpi_vector<size_t> &SA,
                                            const mpi_vector<char> &genome,
                                            const std::string &query) {
    auto comm = genome.comm;
    const int done_tag = 0, cmp_tag = 1, lsa_tag = 2;
    MPI_Request requests[3];

    std::pair<uint64_t, bool> done_buf, done_msg; // received when algorithm stops and holds the result
    std::tuple<uint64_t, uint64_t, uint64_t> cmp_buf, cmp_msg; // left, right, SA[(l+r)/2]
    std::pair<uint64_t, uint64_t> lsa_buf, lsa_msg; // left, right

    // This will determine whether we are the starting process.
    lsa_msg = {0, genome.global_size()};
    bool starter_node = genome.rank == genome.node_with_global_index((lsa_msg.first + lsa_msg.second) / 2);

    bool done = false;
    while (!done) {
        MPI_Irecv(&done_buf, sizeof(done_buf), MPI_BYTE, MPI_ANY_SOURCE, done_tag, comm, &requests[done_tag]);
        MPI_Irecv(&cmp_buf, sizeof(cmp_buf), MPI_BYTE, MPI_ANY_SOURCE, cmp_tag, comm, &requests[cmp_tag]);
        MPI_Irecv(&lsa_buf, sizeof(lsa_buf), MPI_BYTE, MPI_ANY_SOURCE, lsa_tag, comm, &requests[lsa_tag]);

        int msg_tag = lsa_tag; // alos an initial message for starter node
        if (!starter_node) MPI_Waitany(3, requests, &msg_tag, MPI_STATUS_IGNORE);
        else starter_node = false;

        if (msg_tag == done_tag) done_msg = done_buf;
        else if (msg_tag == cmp_tag) cmp_msg = cmp_buf;
        else if (msg_tag == lsa_tag) lsa_msg = lsa_buf;

        bool repeat = true; // allows to repeat the request locally
        while (repeat) {
            repeat = false;

            if (msg_tag == done_tag) {
                done = true;
                break; // TODO
            } else if (msg_tag == cmp_tag) {
                auto [left, right, sa] = cmp_msg;
                size_t lidx = genome.local_index(sa);
                int cmp = strcmp(query.data(), &genome.local_data()[lidx]);

                assert(left <= right);
                if (left == right) {
                    // Broadcast finding of the results
                    for (int dst = 0; dst < genome.nprocs; dst++) {
                        if (dst != genome.rank) {
                            done_msg = {left, cmp == 0};
                            MPI_Send(&done_msg, sizeof(done_msg), MPI_BYTE, dst, done_tag, comm);
                        }
                    }
                    repeat = true;
                    msg_tag = done_tag;
                } else { // left != right
                    uint64_t center = (left + right) / 2;
                    if (cmp < 0) {
                        right = center - 1;
                    } else if (cmp == 0) {
                        right = center;
                    } else {
                        left = center + 1;
                    }

                    auto node = genome.node_with_global_index(center);
                    lsa_msg = {left, right};
                    msg_tag = lsa_tag;
                    repeat = node == genome.rank;
                    if (!repeat) MPI_Send(&lsa_msg, sizeof(lsa_msg), MPI_BYTE, node, lsa_tag, comm);
                }
            } else { // lsa_tag
                auto [left, right] = lsa_msg;
                uint64_t gidx = (left + right) / 2;
                size_t lidx = SA.local_index(gidx);
                auto sa = SA[lidx];

                auto node = SA.node_with_global_index(sa);
                cmp_msg = {left, right, sa};
                msg_tag = cmp_tag;
                repeat = node == SA.rank;
                if (!repeat) MPI_Send(&cmp_msg, sizeof(cmp_msg), MPI_BYTE, node, cmp_tag, comm);
            }
        }
    }

    MPI_Request_free(&requests[done_tag]);
    MPI_Request_free(&requests[cmp_tag]);
    MPI_Request_free(&requests[lsa_tag]);

    return done_msg;
}

// static uint64_t find_lower(const mpi_vector<size_t> &SA,
//                            const mpi_vector<char> &genome,
//                            const std::string &query) {
    
//     const int msg_tag = 0;
//     const int done_tag = 1;
//     char done = 0;
//     MPI_request requests[2];
//     MPI_Irecv(&done, 1, MPI_CHAR, MPI_ANY_SOURCE, done_tag, genome.comm, &requests[done_tag]);
    
//     while (!done) {
//         std::tuple<uint64_t, uint64_t, uint64_t> msg;
//         MPI_Irecv(&interval, sizeof(interval), MPI_BYTE, MPI_ANY_SOURCE, msg_tag, genome.comm, &requests[msg_tag]);

//         int index;
//         MPI_Waitany(2, requests, &index, MPI_STATUS_IGNORE);

//         if (index == done_tag) {
//             done = true;
//             MPI_Request_free(&requests[msg_tag]);
//         } else {
//             // Actually got a requests to process.
//             auto [left, right, gidx] = interval;
//             assert(genome.node_with_global_index(gidx) == genome.rank()); // make sure we are responsible for this index
//             assert (left <= right);
            
//             uint64_t center = (left + right) / 2;
//             auto local_index = genome.local_index(gidx);

//             int cmp = strcmp(genome.local_data().data(), query.data());

//             if (left == right) {
//                 if (cmp == 0) {
//                     // found!
//                 } else {
//                     // not found, stop!
//                 }
//             } else {
//                 if (cmp < 0) {
//                     right = center - 1;
//                 } else if (cmp == 0) {
//                     right = center;
//                 } else {
//                     left = center + 1;
//                 }
//             }

//         }
//     }
// }

static uint64_t num_occurences(const mpi_vector<size_t> &SA,
                               const mpi_vector<char> &genome,
                               const std::string &query) {

    assert(SA.rank == genome.rank);
   auto lower = find_lower(SA, genome, query);
   if (genome.rank == 0) std::cerr << "lower: (" << lower.first << ", " << lower.second << ")\n";
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    auto [num_genomes, num_queries, genome_in, queries_in, queries_out] = parseArgs(argc, argv);
    DataSource data_source(genome_in.data());

    mpi_vector<char> genome(MPI_COMM_WORLD, data_source, 0);
    genome.push_back('$');
    
    std::vector<std::any> dbg_seq, dbg_par;
    auto full_genome_vec = genome.gather(0);
    auto full_genome = std::string(full_genome_vec.data(), full_genome_vec.size());
    if (genome.rank == 0) auto SA_seq = gspf(full_genome, dbg_seq);

    auto SA = suffix_array(genome, dbg_par);

    if (genome.rank == 0) {
        auto num_occurs_seq = num_occurences_naive(full_genome, "CGTC");
        std::cerr << "num occurs seq: " << num_occurs_seq << "\n"; 
    }

    auto num_occurs = num_occurences(SA, genome, "CGTC");

    MPI_Finalize();
    return 0;
}