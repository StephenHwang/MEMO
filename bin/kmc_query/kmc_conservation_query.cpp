#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include "kmc_api/kmc_file.h"
#include "kmc_api/kmer_defs.h"


int query_kmc_database(CKMCFile& kmer_data_base, std::string fasta_record) {
    // std::cout << "Fasta record: " << fasta_record << std::endl;
    std::vector<uint32_t> read_count_vector;
    kmer_data_base.GetCountersForRead(fasta_record, read_count_vector);
    // print counts of canonical kmers
    for (auto counter : read_count_vector)
      std::cout << counter << std::endl;
    return 0;
}

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
// Query kmc database for kmer counts of a each record in a fasta 
//------------------------------------------------------------------------------
    // open kmer_data_base
    CKMCFile kmer_data_base; 
    if (!kmer_data_base.OpenForRA(kmc_database_prefix)) {
        std::cerr << "Error: Unable to open KMC database." << std::endl;
        return 1;
    }

    // read each fasta record, counting k-mers
    std::string line, id, DNA_sequence;
    std::ifstream fasta_file_stream;
    fasta_file_stream.open(fasta_file);
    if (fasta_file_stream.is_open()) {
        while (std::getline(fasta_file_stream, line)) {
            if(line.empty())
                continue;
            if (line[0] == '>') {
                if(!id.empty())
                    query_kmc_database(kmer_data_base, DNA_sequence);
                id = line.substr(1);
                DNA_sequence.clear();
            }
            else {
                DNA_sequence += line;
            }
        }
        if(!id.empty())
            query_kmc_database(kmer_data_base, DNA_sequence);
        fasta_file_stream.close();
    }

    // close kmer_data_base
    kmer_data_base.Close();
    return 0;
}


