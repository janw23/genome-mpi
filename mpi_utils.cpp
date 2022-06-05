#include "mpi_utils.h"

// ################################ MPIContext ################################

MPIContext::MPIContext(DataSource &data_source, MPI_Comm comm) 
: comm(comm), rank(commRank(comm)), nprocs(commSize(comm)), data_chunks_size(nprocs) {
    const size_t genome_idx = 0; // TODO impl should handle all genomes
    const uint64_t local_chunk_size = data_source.getNodeGenomeSize(genome_idx);
    MPI_Alltoall(&local_chunk_size, 1, MPI_UINT64_T, data_chunks_size.data(), 1, MPI_UINT64_T, comm);
}

uint64_t MPIContext::getNodeDataSize(uint64_t node_rank) const {
    return data_chunks_size[node_rank];
}

int MPIContext::commRank(MPI_Comm comm) const {
    int rank;
    MPI_Comm_rank(comm, &rank);
    return rank;
}

int MPIContext::commSize(MPI_Comm comm) const {
    int size;
    MPI_Comm_size(comm, &size);
    return size;
}
