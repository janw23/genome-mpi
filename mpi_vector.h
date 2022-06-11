#ifndef _MPI_VECTOR_H_
#define _MPI_VECTOR_H_

#include <mpi.h>
#include <cinttypes>
#include <vector>
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

    void push_back(T elem);

    const std::vector<T> &local_data() const;

    size_t size() const;

    uint64_t global_index(size_t index) const;

    size_t node_with_global_index(uint64_t index) const;

    // Wrapped MPI operations

    // Gathers all data at process with [rank]. Other processes get empty vector.
    std::vector<T> gather(int rank);
    

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
mpi_vector<T>::mpi_vector(MPI_Comm comm, size_t size) : mpi_vector(comm) {
    datavec = std::vector<T>(size);
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
const std::vector<T> &mpi_vector<T>::local_data() const {
    return datavec;
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
std::vector<T> mpi_vector<T>::gather(int root) {
    if (rank == root) {
        std::vector<T> buff(global_size);
        std::vector<int> recvcounts(chunk_sizes.size());
        std::vector<int> displs(chunk_sizes.size());

        int offset = 0;
        for (size_t i = 0; i < chunk_sizes.size(); i++){
            recvcounts[i] = sizeof(T) * chunk_sizes[i];
            displs[i] = offset;
            offset += recvcounts[i];
        }

        MPI_Gatherv(datavec.data(), sizeof(T) * datavec.size(), MPI_BYTE, buff.data(), recvcounts.data(), displs.data(), MPI_BYTE, root, comm);
        return buff;
    } else {
        MPI_Gatherv(datavec.data(), sizeof(T) * datavec.size(), MPI_BYTE, nullptr, nullptr, nullptr, nullptr, root, comm);
        return std::vector<T>();
    }
}


#endif // _MPI_VECTOR_H_