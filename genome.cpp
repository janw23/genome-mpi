// Author: Jan Wawszczak 418479

#include <mpi.h>
#include <iostream>
#include <vector>
#include <exception>
#include "data_source.h"
#include "mpi_utils.h"

static std::tuple<uint64_t, uint64_t, std::string, std::string, std::string>
parseArgs(int argc, char *argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " :n :m :genome_in :queries_in :queries_out" << std::endl;
        MPI_Finalize();
        exit(1);
    }
    
    try {
        uint64_t n = std::stoull(argv[1]);
        uint64_t m = std::stoull(argv[2]);
        std::string genome_in = argv[3];
        std::string queries_in = argv[4];
        std::string genome_out = argv[5];
        return {n, m, genome_in, queries_in, genome_out};
    } catch (std::exception &e) {
        std::cerr << "parsing input args failed, reason: " << e.what() << std::endl;
        MPI_Finalize();
        exit(1);
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    auto [num_genomes, num_queries, genome_in, queries_in, queries_out] = parseArgs(argc, argv);
    DataSource data_source(genome_in.data());

    // TODO just for tests read a chunk of first genome
    std::vector<char> genome(data_source.getNodeGenomeSize(0));
    data_source.getNodeGenomeValues(0, genome.data());
    
    MPIContext mpi_context(data_source, MPI_COMM_WORLD);
    std::cout << "I am node " << mpi_context.rank << " and got initial genome: " << genome.data() << "\n";

    mpi_sort(genome, std::less<char>(), mpi_context);

    std::cout << "I am node " << mpi_context.rank << " and got sorted genome: " << genome.data() << "\n";

    MPI_Finalize();
    return 0;
}