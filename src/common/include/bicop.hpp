//
// Created by Vatter Thibault on 13/12/16.
//

#ifndef __BICOP_H_INCLUDED__
#define __BICOP_H_INCLUDED__

#include <Eigen/Dense>
using namespace Eigen;

class BiCop
{

public:
    BiCop();
    ~BiCop();

    // hfunctions: the conditioning variable is put second
    virtual VectorXd hfunc1(MatrixXd u) = 0;
    virtual VectorXd hfunc2(MatrixXd u) = 0;

    // PDF
    virtual VectorXd pdf(MatrixXd u) = 0;

    // fit statistics
    virtual double loglik(MatrixXd u) = 0;

};

#endif // __BICOP_H_INCLUDED__
