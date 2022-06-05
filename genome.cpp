// Author: Jan Wawszczak 418479

#include <mpi.h>
#include <iostream>

static std::tuple<uint64_t, uint64_t, std::string, std::string, std::string>
parse_args(int argc, char *argv[]) {
    return {0, 0, "xd", "xd", "dx"};
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    auto [num_genomes, num_queries, genome_in, queries_in, queries_out] = parse_args(argc, argv);

    std::cout << num_genomes << "\n";

    for (int i = 0; i < argc; i++) {
        std::cout << "argv[" << i << "] = " << argv[i] << "\n";
    }
    
    MPI_Finalize();
    return 0;
}