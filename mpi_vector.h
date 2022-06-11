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
    mpi_vector(MPI_Comm comm);
    mpi_vector(MPI_Comm comm, size_t size);
    mpi_vector(MPI_Comm comm, DataSource &data_source, size_t genome_index);

    T &operator[](size_t index) const;

    void push_back(T elem);

    size_t size() const;

    uint64_t global_index(size_t index) const;

    size_t node_with_global_index(uint64_t index) const;
    

private:
    std::vector<T> datavec;
    int rank;
    int nprocs;
    MPI_Comm comm;

    std::vector<size_t> chunk_sizes;
    uint64_t chunk_offset;
    uint64_t global_size;

    void update_chunks();
};

template <typename T>
mpi_vector<T>::mpi_vector(MPI_Comm comm) : comm(comm) {
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);
    chunk_sizes = std::vector<size_t>(nprocs);
}

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
T &mpi_vector<T>::operator[](size_t index) const {
    return datavec[index];
}

template <typename T>
size_t mpi_vector<T>::size() const {
    return datavec.size();
}

template <typename T>
void mpi_vector<T>::push_back(T elem) {
    datavec.push_back(std::move(elem));
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
    assert(rank > 0);
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

#endif // _MPI_VECTOR_H_