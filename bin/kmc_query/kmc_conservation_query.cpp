#include <iostream>
#include <string>
#include <getopt.h>
#include "kmc_api/kmc_file.h"
#include "kmc_api/kmer_defs.h"

int main(int argc, char *argv[]) {
//------------------------------------------------------------------------------
// Parse input parameters for kmer database and fasta file 
//------------------------------------------------------------------------------
    std::string kmc_database_prefix;
    std::string fasta_file;
    int c;
    while (1) {
        static struct option long_options[] = {
            {"d", required_argument, 0, '1'},
            {"f", required_argument, 0, '2'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "1:2:", long_options, &option_index);
        if (c == -1)
            break;
        switch (c) {
            case '1':
                kmc_database_prefix = optarg;
                break;
            case '2':
                fasta_file = optarg;
                break;
            default:
                break;
        }
    }
    if (kmc_database_prefix.empty() || fasta_file.empty()) {
        std::cerr << "Usage: " << argv[0] << " --d <kmc_database_prefix> --f <fasta>" << std::endl;
        return 1;
    }
    // std::cout << "KMC database prefix: " << kmc_database_prefix << std::endl;
    // std::cout << "Fasta file: " << fasta_file << std::endl;

//------------------------------------------------------------------------------
// Extract read from fasta file
//------------------------------------------------------------------------------
// TODO: parse single line fasta_file


//------------------------------------------------------------------------------
// Query kmer counts of a read to a kmer database
//------------------------------------------------------------------------------
    // Open the KMC database file for reading
    CKMCFile kmer_data_base; 
    if (!kmer_data_base.OpenForRA(kmc_database_prefix)) {
        std::cerr << "Error: Unable to open KMC database." << std::endl;
        return 1;
    }
    else
    {
        std::vector<uint32_t> read_count_vector;
        kmer_data_base.GetCountersForRead(
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            read_count_vector);

        // print counts of canonical kmers
        for (auto counter : read_count_vector)
          std::cout << counter << std::endl;

        // Close the KMC database
        kmer_data_base.Close();
    }

    return 0;
}


