#ifndef readArray_h
#define readArray_h

#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <stdexcept>

using namespace std;
using namespace Eigen;

#define MAXBUFSIZE  ((int) 1e6)

ArrayXXd readArray(const char *filename)
    {
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];

    // Read numbers from file into buffer.
    ifstream infile;
    infile.open(filename);

    if(infile.is_open()==false){//Check if file opens correctly - most likely problem is it not existing, but could be others
        std::string error_message("\nFile not found or file error!\nFilename: \n ");
        error_message.append( filename ) ;
        throw std::invalid_argument(error_message );
        
        ArrayXXd error(1,1);
        return error;
        }
    while (! infile.eof())
        {
        string line;
        getline(infile, line);

        int temp_cols = 0;
        stringstream stream(line);
        while(stream >> buff[cols*rows+temp_cols]) //! stream.eof()
            temp_cols++;
        
        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;
        rows++;
        }

    infile.close();

    //cols--;
    //rows--;
    // Populate matrix with numbers.
    ArrayXXd result(rows,cols);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[ cols*i+j ];
    
    return result;
    };


#endif