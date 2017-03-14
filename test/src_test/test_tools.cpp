// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "include/test_tools.hpp"

MatXi read_matxi(const char *filename, int max_buffer_size)
{
    MatXd temp = read_matxd(filename, max_buffer_size);
    MatXi output = temp.cast <int> ();
    return output;
}
MatXd read_matxd(const char *filename, int max_buffer_size)
{
    using namespace std;

    int cols = 0, rows = 0;
    double* buff = new double[max_buffer_size];

    // Read numbers from file into buffer.
    ifstream infile;
    infile.open(filename);
    while (! infile.eof())
    {
        string line;
        getline(infile, line);

        int temp_cols = 0;
        stringstream stream(line);
        while(! stream.eof())
            stream >> buff[cols*rows+temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
    }

    infile.close();

    rows--;

    // Populate matrix with numbers.
    MatXd result(rows,cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[ cols*i+j ];

    delete [] buff;
    return result;
};
