#ifndef _MPI_UTILS_HPP_
#define _MPI_UTILS_HPP_

#include <vector>
#include <mpi.h>
#include <algorithm>
#include "data_source.h"

// https://stackoverflow.com/questions/60308967/an-error-occured-in-mpi-recv-while-sending-an-array
#define MPI_Error_Check(x) {const int err=x; if(x!=MPI_SUCCESS) { fprintf(stderr, "MPI ERROR %d at %d.", err, __LINE__);}}

class MPIContext {
public:
    const MPI_Comm comm;
    const int rank;
    const int nprocs;

    MPIContext(DataSource &data_source, MPI_Comm comm);
    uint64_t getNodeGenomeSize(int node_rank) const;
    uint64_t getNodeGenomeOffset(int node_rank) const;
    uint64_t getNodeWithIndexOfData(uint64_t idx) const;

    bool rankFirst() const;
    bool rankLast() const;

private:
    std::vector<uint64_t> genome_sizes_prefix_sums;

    int commRank(MPI_Comm comm) const;
    int commSize(MPI_Comm comm) const;
};

template<typename T, typename Compare>
void mpiSort(std::vector<T> &data_chunk, Compare comp, const MPIContext &mpi_context) {
    const size_t chunk_size = data_chunk.size();

    if (mpi_context.rank == 0) {
        // Gather chunk sizes at root.
        std::vector<uint64_t> node_chunk_size(mpi_context.nprocs);
        MPI_Gather(&chunk_size, 1, MPI_UINT64_T, node_chunk_size.data(), 1, MPI_UINT64_T, 0, mpi_context.comm);

        std::vector<T> full_data(data_chunk);

        // Receive full data from other processes.
        for (int proc = 1; proc < mpi_context.nprocs; proc++) {
            std::vector<T> proc_data(node_chunk_size[proc]);
            MPI_Error_Check(MPI_Recv(static_cast<void*>(proc_data.data()), sizeof(T) * proc_data.size(), MPI_BYTE, proc, MPI_ANY_TAG, mpi_context.comm, MPI_STATUS_IGNORE));
            std::copy(proc_data.begin(), proc_data.end(), std::back_inserter(full_data));
        }

        // Sort full data.
        std::sort(full_data.begin(), full_data.end(), comp);

        std::copy(full_data.begin(), full_data.begin() + data_chunk.size(), data_chunk.begin());
        size_t offset = data_chunk.size();

        // Send data chunks to other processes.
        for (int proc = 1; proc < mpi_context.nprocs; proc++) {
            MPI_Error_Check(MPI_Send(static_cast<void*>(full_data.data() + offset), sizeof(T) * node_chunk_size[proc], MPI_BYTE, proc, 0, mpi_context.comm));
            offset += node_chunk_size[proc];
        }
    } else {
        MPI_Gather(&chunk_size, 1, MPI_UINT64_T, nullptr, 1, MPI_UINT64_T, 0, mpi_context.comm); // Send chunk sizes to root.
        MPI_Error_Check(MPI_Send(static_cast<void*>(data_chunk.data()), sizeof(T) * data_chunk.size(), MPI_BYTE, 0, 0, mpi_context.comm));
        MPI_Error_Check(MPI_Recv(static_cast<void*>(data_chunk.data()), sizeof(T) * data_chunk.size(), MPI_BYTE, 0, MPI_ANY_TAG, mpi_context.comm, MPI_STATUS_IGNORE));
    }
}

#endif //_MPI_UTILS_HPP_