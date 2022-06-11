#ifndef _MPI_UTILS_H_
#define _MPI_UTILS_H_

#include <vector>
#include <algorithm>
#include "mpi_vector.h"


template<typename T, typename Compare>
void mpi_sort(mpi_vector<T> &data, Compare comp) {
    auto full_data = data.gather(0);
    if (data.rank == 0) std::sort(full_data.begin(), full_data.end(), comp);
    data.scatter(full_data, 0);
}

#endif //_MPI_UTILS_H_