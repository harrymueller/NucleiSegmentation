#include "bin_genes.hpp"

void BinGenes::bin_genes(Mat im, std::string gem_file)
{
    // init array of size <max_ids+1>
    int max_bin_id = IO::find_max_bin_id(im);
    int counts[max_bin_id+1];
    BinGenes::empty_array(counts);

    // 


    std::ofstream file;
    file.open(filename);

    file
}

void BinGenes::empty_array(int* counts, int length)
{
    for ( int i = 0; i < length; i++ ) counts[i] = 0;
}

void BinGenes::save_genes(std::string gene, std::string filename, int* counts, int length) 
{
    std::ofstream file;
    file.open(filename);

    for ( int i = 0; i < length; i++ ) {
        if (counts[i] > 0)
            file << gene << "\t" << i << "\t" << counts[i] << "\n";
    }

    file.close();
}