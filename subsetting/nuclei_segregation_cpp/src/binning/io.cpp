#include "io.hpp"

/*
    Given an image csv containing segments, creates an image object
*/
Image IO::read_markers(std::string filename)
{
    // template variables
    std::string l, val;
    std::vector<std::vector<int>> data;
    std::vector<int> row;
    int i, j;

    // open file and check that it is valid
    std::fstream file(filename);
    if (!file.is_open()) {
        std::cout << "Could not open file '" << filename << "'\n";
        exit(EXIT_FAILURE);
    }

    // read line by line
    while (std::getline(file, l)) {
        row.clear();

        std::stringstream temp(l);
        while (getline(temp, val, ',')) {
            row.push_back(atoi(val.c_str()));
        }
        data.push_back(row);
    }

    // Now add all the data into a Mat element
    Mat im = Mat::zeros((int) data.size(), (int) data[0].size(), CV_32S);
    // Loop over vectors and add the data
    for ( i = 0; i < (int) data.size(); i++ ){
        for ( j= 0; j< (int) data[0].size(); j++ ){
            im.at<int>(i,j) = data[i][j];
        }
    }

    return Image(im);
}

/*
    Creates a sparse csv containing <bin_id  x  y>
*/
void IO::sparse_segments(Mat im, std::string filename)
{
    // var declaration
    int val, i, j;

    // output file stream
    std::ofstream file;
    file.open(filename);

    // header
    file << "bin_id\tx\ty\n";
    // loop through image matrix saving labels
    for ( i = 0; i < im.rows; i++ ) {
        for ( j = 0; j < im.cols; j++ ) {
            val = im.at<int>(i,j);
            if (val > 1)
                file << val << "\t" << j << "\t" << i << "\n";
        }
    }

    file.close();
}

/*
    Finds the mean coordinates for each bin_id
*/
void IO::bin_coords(Mat im, std::string output)
{
    // var declaration
    int i, j, val;

    // get max bin id
    int max_id = find_max_bin_id(im);

    // create an array to store <sum x, sum y, n>
    long data[max_id+1][3];
    for ( i = 0; i < max_id+1; i++) {
        data[i][0] = 0; data[i][1] = 0; data[i][2] = 0; // set array vals to 0
    };

    // fill array
    for ( i = 0; i < im.rows; i++ ) { // y
        for ( j = 0; j < im.cols; j++ ) { // x
            val = im.at<int>(i,j);
            if (val > 1) {
                data[val][0] += j;
                data[val][1] += i;
                data[val][2] += 1;
            }
        }
    }

    // output file stream
    std::ofstream file;
    file.open(output);

    // header
    file << "bin_id,x,y\n";
    for ( i = 2; i < max_id + 1; i++ ) {
        if (data[i][2] == 0) std::cout << i << " "; // not bins found with that id
        else file << i << "," << data[i][0] * 1.0 / data[i][2] << "," << data[i][1] * 1.0 / data[i][2] << "\n";
    }
    std::cout << "\n";
    file.close();
}

/*
    Finds the largest bin id, searches from bottom right to top right, rows first
*/
int IO::find_max_bin_id(Mat im) {
    int i, j, val;
    for ( i = im.rows - 1; i >= 0; i-- ) { // y
        for ( j = im.cols - 1; j >= 0; j-- ) { // x
            val = im.at<int>(i,j);
            if (val > 1)
                return val;
        }
    }
    return -1;
}