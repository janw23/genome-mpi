#ifndef _MPI_VECTOR_H_
#define _MPI_VECTOR_H_

#include <mpi.h>
#include <cinttypes>
#include <vector>
#include <algorithm>
#include <cassert>
#include "data_source.h"

template <typename T>
class mpi_vector {
public:
    const MPI_Comm comm;
    const int rank;
    const int nprocs;

    mpi_vector(MPI_Comm comm);
    mpi_vector(MPI_Comm comm, size_t size);
    mpi_vector(MPI_Comm comm, DataSource &data_source, size_t genome_index);
    template <typename T2>
    mpi_vector(const mpi_vector<T2> &vec);

    T &operator[](size_t index);
    const T &operator[](size_t index) const;

    void push_back(T elem);

    const std::vector<T> &local_data() const;

    typename std::vector<T>::const_iterator cbegin() const;
    typename std::vector<T>::const_iterator cend() const;

    size_t size() const;

    uint64_t global_index(size_t index) const;
    uint64_t global_offset() const;

    size_t node_with_global_index(uint64_t index) const;

    // Iterate over whole vector, at each step applying
    // a function (context, local_index) -> void
    // and passing the context between the nodes.
    template <typename Context, typename Func>
    void iter_apply(Context ctx, Func f);

    // Element [local_idx] goes to position global_indices[local_idx] in global array.
    void reorder(const mpi_vector<uint64_t> &global_indices);

    // Wrapped standard MPI operations

    // Gathers all data at process with [root]. Other processes get empty vector.
    std::vector<T> gather(int root) const;

    // Distributes [data] to processes according to chunk_sizes.
    // [data] is insignificant at non-root processes.
    // [data] must have the same size as global_size.
    void scatter(const std::vector<T> &data, int root);
    


private:
    std::vector<T> datavec;

    std::vector<size_t> chunk_sizes;
    uint64_t chunk_offset;
    uint64_t global_size;

    void update_chunks();
};

static int get_rank(MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);
    return rank;
}

static int get_nprocs(MPI_Comm comm) {
    int nprocs;
    MPI_Comm_size(comm, &nprocs);
    return nprocs;
}

template <typename T>
mpi_vector<T>::mpi_vector(MPI_Comm comm)
: comm(comm), rank(get_rank(comm)), nprocs(get_nprocs(comm)), chunk_sizes(nprocs) {}

template <typename T>
mpi_vector<T>::mpi_vector(MPI_Comm comm, size_t local_size) : mpi_vector(comm) {
    datavec = std::vector<T>(local_size);
    update_chunks();
}

template <>
mpi_vector<char>::mpi_vector(MPI_Comm comm, DataSource &data_source, size_t genome_index) 
: mpi_vector(comm) {
    datavec = std::vector<char>(data_source.getNodeGenomeSize(genome_index));
    data_source.getNodeGenomeValues(genome_index, datavec.data());
    update_chunks();
}

template <typename T>
template <typename T2>
mpi_vector<T>::mpi_vector(const mpi_vector<T2> &vec) : mpi_vector(vec.comm) {
    datavec = std::vector<T>(vec.local_data().size());
    // We copy this way to allow type conversion.
    std::copy(vec.local_data().cbegin(), vec.local_data().cend(), datavec.begin());
    update_chunks();
}

template <typename T>
T &mpi_vector<T>::operator[](size_t index) {
    return datavec[index];
}

template <typename T>
const T &mpi_vector<T>::operator[](size_t index) const {
    return datavec[index];
}

template <typename T>
const std::vector<T> &mpi_vector<T>::local_data() const {
    return datavec;
}

template <typename T>
typename std::vector<T>::const_iterator mpi_vector<T>::cbegin() const {
    return datavec.cbegin();
}

template <typename T>
typename std::vector<T>::const_iterator mpi_vector<T>::cend() const {
    return datavec.cend();
}

template <typename T>
size_t mpi_vector<T>::size() const {
    return datavec.size();
}

template <typename T>
void mpi_vector<T>::push_back(T elem) {
    if (rank == nprocs - 1) {
        datavec.push_back(std::move(elem));
    }
    update_chunks();
}

template <typename T>
uint64_t mpi_vector<T>::global_index(size_t index) const {
    return chunk_offset + index;
}

template <typename T>
uint64_t mpi_vector<T>::global_offset() const {
    return global_index(0);
}

template <typename T>
uint64_t mpi_vector<T>::node_with_global_index(uint64_t global_index) const {
    assert(global_index < global_size);

    // This brutal solution is okay because there are not that many nodes.
    uint64_t offset = 0;
    for (size_t i = 0; i < chunk_sizes.size(); i++) {
        if (offset <= global_index && global_index < offset + chunk_sizes[i]) {
            return i;
        }
        offset += chunk_sizes[i];
    }

    assert(false);
    return 0;
}

template <typename T>
void mpi_vector<T>::update_chunks() {
    assert(rank >= 0);
    assert(nprocs > 0);
    assert(chunk_sizes.size() == static_cast<size_t>(nprocs));

    // Gather all chunk sizes from other nodes.
    size_t chunk_size = size();
    MPI_Allgather(&chunk_size, 1, MPI_UINT64_T, chunk_sizes.data(), 1, MPI_INT64_T, comm);

    uint64_t offset = 0;
    for (size_t i = 0; i < chunk_sizes.size(); i++) {
        if (i == static_cast<size_t>(rank)) chunk_offset = offset;
        offset += chunk_sizes[i];
    }
    global_size = offset;
}

template <typename T>
template <typename Context, typename Func>
void mpi_vector<T>::iter_apply(Context ctx, Func f) {
    if (rank > 0) {
        MPI_Recv(&ctx, sizeof(ctx), MPI_BYTE, rank - 1, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
    }
    for (size_t i = 0; i < size(); i++) {
        f(ctx, i);
    }
    if (rank < nprocs - 1) {
        MPI_Send(&ctx, sizeof(ctx), MPI_BYTE, rank + 1, 0, comm);
    }
}

template <typename T>
static std::pair<std::vector<int>, std::vector<int>>
counts_displacements(const std::vector<size_t> &chunk_sizes) {
    std::vector<int> counts(chunk_sizes.size());
    std::vector<int> displs(chunk_sizes.size());

    int offset = 0;
    for (size_t i = 0; i < chunk_sizes.size(); i++){
        counts[i] = sizeof(T) * chunk_sizes[i];
        displs[i] = offset;
        offset += counts[i];
    }

    return {counts, displs};
}

// TODO remove
template <typename T>
static void printvec(std::string name, std::vector<T> const &vec) {
    return;
    std::cerr << name << ":";
    for (T e : vec) std::cerr << " " << e;
    std::cerr << "\n";
}

template <typename T>
void mpi_vector<T>::reorder(const mpi_vector<uint64_t> &global_indices) {
    assert(size() == global_indices.size());

    // Count the number of elements to send to each process from this one.
    std::vector<size_t> sendcounts(nprocs);
    for (size_t i = 0; i < size(); i++) {
        sendcounts[node_with_global_index(global_indices[i])]++; // TODO this might be slow due to brutal impl
    }

    std::vector<size_t> bucket_displacements(nprocs);
    for (size_t i = 1; i < sendcounts.size(); i++) {
        bucket_displacements[i] = bucket_displacements[i-1] + sendcounts[i-1];
    }

    // Put data in buckets depending on their destination process.
    std::vector<std::pair<T, uint64_t>> sendbuf(size()); // contains data and global index
    for (size_t i = 0; i < size(); i++) {
        auto dst = node_with_global_index(global_indices[i]);
        sendbuf[bucket_displacements[dst]] = {datavec[i], global_indices[i]};
        bucket_displacements[dst]++;
    }
    
    // Determine how much data each process will receive from each one.
    std::vector<size_t> recvcounts(sendcounts.size());
    MPI_Alltoall(sendcounts.data(), 1, MPI_UINT64_T, recvcounts.data(), 1, MPI_UINT64_T, comm);

    // Exchange the actual data.
    auto [send_bytes, send_displs] = counts_displacements<std::pair<T, uint64_t>>(sendcounts);
    auto [recv_bytes, recv_displs] = counts_displacements<std::pair<T, uint64_t>>(recvcounts);

    std::vector<std::pair<T, uint64_t>> recvbuf(sendbuf.size());
    MPI_Alltoallv(
        sendbuf.data(), send_bytes.data(), send_displs.data(), MPI_BYTE,
        recvbuf.data(), recv_bytes.data(), recv_displs.data(), MPI_BYTE, comm
    );

    // Place the received data in the correct index locally.
    for (const auto &data : recvbuf) {
        size_t local_idx = data.second - global_offset();
        datavec[local_idx] = data.first;
    }
}

template <typename T>
std::vector<T> mpi_vector<T>::gather(int root) const {
    if (rank == root) {
        std::vector<T> buff(global_size);
        auto [counts, displs] = counts_displacements<T>(chunk_sizes);

        MPI_Gatherv(
            datavec.data(), sizeof(T) * datavec.size(), MPI_BYTE,
            buff.data(), counts.data(), displs.data(), MPI_BYTE, root, comm
        );
        return buff;
    } else {
        MPI_Gatherv(
            datavec.data(), sizeof(T) * datavec.size(), MPI_BYTE,
            nullptr, nullptr, nullptr, nullptr, root, comm
        );
        return std::vector<T>();
    }
}

template <typename T>
void mpi_vector<T>::scatter(const std::vector<T> &data, int root) {
    if (rank == root) {
        auto [counts, displs] = counts_displacements<T>(chunk_sizes);
        MPI_Scatterv(
            data.data(), counts.data(), displs.data(), MPI_BYTE,
            datavec.data(), sizeof(T) * datavec.size(), MPI_BYTE, root, comm
        );
    } else {
        MPI_Scatterv(
            nullptr, nullptr, nullptr, nullptr,
            datavec.data(), sizeof(T) * datavec.size(), MPI_BYTE, root, comm
        );
    }
}


#endif // _MPI_VECTOR_H_