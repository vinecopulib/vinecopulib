//
// Created by Vatter Thibault on 13/12/16.
//

#ifndef VINECOPLIB_BICOP_H_INCLUDED__
#define VINECOPLIB_BICOP_H_INCLUDED__

#include <Eigen/Dense>
typedef Eigen::VectorXd VecXd;
typedef Eigen::MatrixXd MatXd;

class Bicop
{

public:

    // hfunctions: the conditioning variable is put second
    virtual VecXd hfunc1(const MatXd *u) = 0;
    virtual VecXd hfunc2(const MatXd *u) = 0;

    // PDF
    virtual VecXd pdf(const MatXd *u) = 0;

    // fit statistics
    //virtual double loglik(MatXd u) = 0;

};

#endif // __BICOP_H_INCLUDED__
