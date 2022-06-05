#ifndef _MPI_UTILS_HPP_
#define _MPI_UTILS_HPP_

#include <vector>
#include <mpi.h>
#include <algorithm>
#include "data_source.h"

class MPIContext {
public:
    const MPI_Comm comm;
    const int rank;
    const int nprocs;

    MPIContext(DataSource &data_source, MPI_Comm comm);
    uint64_t getNodeDataSize(uint64_t node_rank) const;

private:
    std::vector<uint64_t> data_chunks_size;

    int commRank(MPI_Comm comm) const;
    int commSize(MPI_Comm comm) const;
};

template<typename T, typename Compare>
void mpi_sort(std::vector<T> &data_chunk, Compare comp, const MPIContext &mpi_context) {
    if (mpi_context.rank == 0) {
        std::vector<char> full_genome(data_chunk);

        // Receive full genome from other processes.
        for (int proc = 1; proc < mpi_context.nprocs; proc++) {
            std::vector<char> proc_data(mpi_context.getNodeDataSize(proc));
            MPI_Recv(proc_data.data(), proc_data.size(), MPI_CHAR, proc, MPI_ANY_TAG, mpi_context.comm, MPI_STATUS_IGNORE);
            std::copy(proc_data.begin(), proc_data.end(), std::back_inserter(full_genome));
        }

        // Sort full genome.
        std::sort(full_genome.begin(), full_genome.end(), comp);

        std::copy(full_genome.begin(), full_genome.begin() + data_chunk.size(), data_chunk.begin());
        size_t offset = data_chunk.size();

        // Send full genome chunks to other processes.
        for (int proc = 1; proc < mpi_context.nprocs; proc++) {
            uint64_t target_node_chunk_size = mpi_context.getNodeDataSize(proc);
            MPI_Send(full_genome.data() + offset, target_node_chunk_size, MPI_CHAR, proc, 0, mpi_context.comm);
            offset += target_node_chunk_size;
        }
    } else {
        MPI_Send(data_chunk.data(), data_chunk.size(), MPI_CHAR, 0, 0, mpi_context.comm);
        MPI_Recv(data_chunk.data(), data_chunk.size(), MPI_CHAR, 0, MPI_ANY_TAG, mpi_context.comm, MPI_STATUS_IGNORE);
    }
}

#endif //_MPI_UTILS_HPP_