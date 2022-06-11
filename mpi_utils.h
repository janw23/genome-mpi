#ifndef _MPI_UTILS_H_
#define _MPI_UTILS_H_

#include <vector>
#include <algorithm>
#include "mpi_vector.h"


template<typename T, typename Compare>
mpi_vector<T> &mpi_sort(mpi_vector<T> &data, Compare comp) {
    auto full_data = data.gather(0);
    if (data.rank == 0) std::sort(full_data.begin(), full_data.end(), comp);
    data.scatter(full_data, 0);
    return data; // allows chaining functions
}

template <typename T>
mpi_vector<std::pair<T, uint64_t>> indexed(const mpi_vector<T> &vec) {
    mpi_vector<std::pair<T, uint64_t>> tmp(vec.comm, vec.size());
    for (size_t i = 0; i < vec.size(); i++) {
        tmp[i].first = vec[i];
        tmp[i].second = vec.global_offset() + i;
    }
    return tmp;
}

template <typename T>
void unpack_indexed(const mpi_vector<std::pair<T, uint64_t>> &vec,
                           mpi_vector<T> &vals, mpi_vector<uint64_t> &idxs) {
    assert(vals.size() == idxs.size());
    assert(vals.size() == vec.size());

    for (size_t i = 0; i < vec.size(); i++) {
        vals[i] = vec[i].first;
        idxs[i] = vec[i].second;
    }
}


#endif //_MPI_UTILS_H_