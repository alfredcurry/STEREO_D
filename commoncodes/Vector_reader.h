#ifndef readVector_h
#define readVector_h

#include <stdexcept>
#include <vector>
#include <string>

template <typename Real>
void readVector(std::vector<Real> &vec , const int &ROWS , std::string &filename){
    std::ifstream fin;
    Real item;
    fin.open(filename);
    if(!fin.is_open()){
            std::string error_message("\nCoefficient file not found or file error!\nFilename: \n ");
            error_message.append( filename ) ;
            throw std::invalid_argument(error_message );
    }

    for (int row = 0; row < ROWS; row++)
    {
        fin >> item;
        vec[row] = item ; 
    }
    fin.close();

}

#endif