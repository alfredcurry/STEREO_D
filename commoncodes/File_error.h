#ifndef File_Error_H
#define File_Error_H

static void file_error(std::string filename){
    std::string error_message("\nCoefficient file not found or file error!\nFilename: \n ");
    error_message.append( filename ) ;
    throw std::invalid_argument(error_message );   
}
#endif