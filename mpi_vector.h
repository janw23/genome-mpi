#ifndef _MPI_VECTOR_H_
#define _MPI_VECTOR_H_

#include <mpi.h>
#include <cinttypes>
#include <vector>
#include <cassert>

template <typename T>
class mpi_vector {
public:
    mpi_vector(MPI_COMM comm);
    mpi_vector(MPI_COMM comm, size_t size);

    T &operator[](size_t index) const;

    void push_back(T elem);

    size_t size() const;

    uint64_t global_index(size_t index) const;

    size_t node_with_global_index(uint64_t index) const;
    

private:
    std::vector<T> datavec;
    int rank;
    int nprocs;
    MPI_COMM comm;

    std::vector<size_t> chunk_sizes;
    uint64_t chunk_offset;
    uint64_t global_size;

    void update_chunks();
}

mpi_vector::mpi_vector(MPI_COMM comm) : comm(comm) {
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);
    update_chunks();
}

mpi_vector::mpi_vector(MPI_COMM comm, size_t size) : mpi_vector(comm), datavec(size) {
    update_chunks();
}

template <typename T>
T &mpi_vector::operator[](size_t index) const {
    return datavec[index];
}

size_t mpi_vector::size() const {
    return datavec.size();
}

template <typename T>
void mpi_vector::push_back(T elem) {
    datavec.push_back(std::move(elem));
    update_chunks();
}

uint64_t mpi_vector::global_index(size_t index) const {
    return chunk_offset + index;
}

uint64_t mpi_vector::node_with_global_index(uint64_t global_index) const {
    assert(global_index < global_size);

    // This brutal solution is okay because there are not that many nodes.
    uint64_t offset = 0;
    for (size_t i = 0; i < chunk_sizes; i++) {
        if (offset <= global_index && global_index < offset + chunk_sizes[i]) {
            return i;
        }
        offset += chunk_sizes[i];
    }

    assert(false);
    return 0;
}

#endif // _MPI_VECTOR_H_