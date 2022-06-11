#include "mpi_utils.h"

// ################################ MPIContext ################################

MPIContext::MPIContext(DataSource &data_source, MPI_Comm comm) 
: comm(comm), rank(commRank(comm)), nprocs(commSize(comm)), genome_sizes_prefix_sums(nprocs + 1, 0) {
    const size_t genome_idx = 0; // TODO impl should handle all genomes
    const uint64_t local_chunk_size = data_source.getNodeGenomeSize(genome_idx) + rankLast(); // last node gets additional field for guard character
    MPI_Error_Check(MPI_Allgather(&local_chunk_size, 1, MPI_UINT64_T, genome_sizes_prefix_sums.data() + 1, 1, MPI_UINT64_T, comm));

    for (size_t i = 1; i < genome_sizes_prefix_sums.size(); i++) {
        genome_sizes_prefix_sums[i] += genome_sizes_prefix_sums[i-1];
    }
}

uint64_t MPIContext::getNodeGenomeSize(int node_rank) const {
    return genome_sizes_prefix_sums[node_rank + 1] - genome_sizes_prefix_sums[node_rank];
}

uint64_t MPIContext::getNodeGenomeOffset(int node_rank) const {
    return genome_sizes_prefix_sums[node_rank];
}

uint64_t MPIContext::getNodeWithIndexOfData(uint64_t idx) const {
    
}

bool MPIContext::rankFirst() const {
    return rank == 0;
}

bool MPIContext::rankLast() const {
    return rank == nprocs - 1;
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
