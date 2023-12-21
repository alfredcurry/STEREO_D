#ifndef _INTERPOLATOR_H
#define _INTERPOLATOR_H

#include <Eigen/Dense>


Eigen::ArrayXd interpolate(const Eigen::ArrayXd &m , const Eigen::ArrayXd &m_grid , const Eigen::ArrayXd &y_grid )
{   
    int J = m.rows() ; 
    Eigen::ArrayXd::Index ind, indup; 
    Eigen::ArrayXd y_interp(J);
    double a ;
    
    for( int j=0; j < J ; j++){

        a = ( m(j) - m_grid ).abs().minCoeff(&ind);

        if (a==0){

            y_interp(j) = y_grid(ind);

        }
        else {

            if ( m(j) > m_grid(ind) ){
                indup = ind+1;
            }
            else if ( m(j) < m_grid(ind) ){
                indup = ind;
                ind = ind-1;
            }
            if (ind < 0){
                y_interp.row(j) = (m(j))*y_grid(indup) /(m_grid(indup));
            }else if( indup >= m_grid.size()){
                y_interp.row(j) = (1-m(j))*y_grid(ind)/(1-m_grid(ind));

            }else{
                y_interp.row(j) = ((m(j)-m_grid(ind))*y_grid(indup) + (m_grid(indup)-m(j))*y_grid(ind))/(m_grid(indup)-m_grid(ind));
            }
        }
    }

    return y_interp;
}
#endif