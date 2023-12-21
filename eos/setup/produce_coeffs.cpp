#include "2Dinterpolator.cpp"
#include "EquationOfStateIdeal.hpp"
#include <string>
#include <fstream>
#include "Matrix_reader.h"

int main(int arg , char** argv){

    std::string x_name = "P", y_name = "T", f_name(argv[1]) , eos_name(argv[2]) ;
    std::string filebase , P_string , T_string;
    char mu_string[10];
    double mu_m = 20;
    int Pgrid_size = std::strtod(argv[3] , nullptr ) , Tgrid_size = std::strtod(argv[4] , nullptr ) ;
    //snprintf(mu_string,10 , "mu%.1f", mu_m);
    filebase = "eos/"+eos_name;
    //filebase.append(mu_string);
    std::cout << filebase << std::endl;
    std::ofstream fout;

    Interpolators::BicubicInterpolator<double> interp ;

    Eigen::ArrayXd P , T ;

    P_string = (filebase+"/grid/"+x_name+std::to_string(Pgrid_size)+".txt");
    T_string = (filebase+"/grid/"+y_name+std::to_string(Tgrid_size)+".txt");
    P = readArray(&P_string[0]); //Eigen::ArrayXd::LinSpaced(1000, 0, 11.5 );
    T =  readArray(&T_string[0]); //Eigen::ArrayXd::LinSpaced(500 , 0 , 4);

    Eigen::ArrayXXd log_rho , rho_dP , rho_dT , rho_dPdT , P_grid(P.size(),T.size()), T_grid(P.size(),T.size());

   // P = pow(10,P);
    //T = pow(10,T);

    for(int i=0 ; i < T.size() ; i++){
        P_grid.col(i) = P;
    }
    for(int i=0 ; i < P.size() ; i++){
        T_grid.row(i) = T;
    }

    if( mu_m ==1 || mu_m ==2){
    //For the ideal case
        class Ideal<Eigen::ArrayXXd> ideal;
        std::vector<double> Pv , Tv ;
        for(int j=0; j< P.size(); j++){
            Pv.push_back(P(j));
        }
        for(int j=0; j< T.size(); j++){
            Tv.push_back(T(j));
        }

        if (f_name == "del"){
            log_rho = log10(ideal.del(pow(10,P_grid),pow(10,T_grid) , mu_m));
            rho_dP = pow(10,P_grid-log_rho)*ideal.DiffdelP(pow(10,P_grid),pow(10,T_grid), mu_m);
            rho_dT = pow(10,T_grid-log_rho)*ideal.DiffdelT(pow(10,P_grid),pow(10,T_grid), mu_m);
            rho_dPdT = pow(10,T_grid+P_grid-log_rho)*log(10)*ideal.Diff2delPT(pow(10,P_grid), pow(10,T_grid) , mu_m);
        }else
        if (f_name == "logrho"){
            std::cout << "it's rho"<<std::endl;

            log_rho = log10(ideal.EoSrho(pow(10,P_grid), pow(10,T_grid) , mu_m));
            std::cout << "rho done "<<std::endl;
            rho_dP = pow(10,P_grid-log_rho)*ideal.DiffrhoP(pow(10,P_grid), pow(10,T_grid) , mu_m);
            rho_dT = pow(10,T_grid-log_rho)*ideal.DiffrhoT(pow(10,P_grid), pow(10,T_grid), mu_m);
            rho_dPdT = pow(10,T_grid+P_grid-log_rho)*log(10)*ideal.Diff2rhoPT(pow(10,P_grid),pow(10,T_grid) , mu_m);  
            std::cout << "all done "<<std::endl;

        }
        fout.open(filebase+"/grid/"+x_name+std::to_string(P.size())+".txt");
            fout << P;
        fout.close();
        fout.open(filebase+"/grid/"+y_name+std::to_string(T.size())+".txt");
            fout << T;
        fout.close();
        std::cout << "PT saved "<<std::endl;
    }
    

    if(mu_m == 0){
    fout.open(filebase+"/data/"+f_name+std::to_string(P.size())+"x"+std::to_string(T.size())+".txt");
        fout << log_rho;
    fout.close();
    std::cout << "rho saved "<<std::endl;

    fout.open(filebase+"/data/"+f_name+"_d"+x_name+".txt");
        fout << rho_dP;
    fout.close();
    fout.open(filebase+"/data/"+f_name+"_d"+y_name+".txt");
        fout << rho_dT;
    fout.close();
    fout.open(filebase+"/data/"+f_name+"_d"+x_name+"d"+y_name+".txt");
        fout << rho_dPdT;
    fout.close();
    std::cout << "all saved "<<std::endl;
    }

    //interp.setupInterpolator(Pv,Tv,log_rho,rho_dP,rho_dT,rho_dPdT);
    int grid_size[2] = {Pgrid_size,Tgrid_size};
    interp.read_setup(filebase, x_name , y_name , f_name , grid_size);
    interp.save_coeff(filebase, f_name);
    interp.read_coeff(filebase);
    //interp.setup_from_coeff(filebase, x_name, y_name, f_name);
    double f_new , f_new_true;

    int P_value =100 , T_value =1873;
    int *indexes ;
    
    indexes = interp.return_index(P_value, T_value);
    f_new = interp.return_value(P_value, T_value, indexes[0],  indexes[1]);
    
    if( f_name == "log_rho"){
        f_new_true = T_value/Ideal<double>::EoSrho(P_value , T_value , mu_m)*Ideal<double>::DiffrhoT(P_value , T_value , mu_m);
    }else
    if(f_name == "del"){
        f_new_true = Ideal<double>::del(P_value , T_value , mu_m);

    }
    
    std::cout << f_new << std::endl;
    //<< " " << f_new[1] << " " << f_new_true << std::endl;
    //std::cout << rho_dT(indexes[0],indexes[1]) << " " << rho_dT(indexes[0]+1,indexes[1]) << std::endl;
    //std::cout << rho_dT(indexes[0],indexes[1]+1) << " " << rho_dT(indexes[0]+1,indexes[1]+1) << std::endl;

    return 0;
}