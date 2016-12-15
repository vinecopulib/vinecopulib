//
// Created by Vatter Thibault on 13/12/16.
//

#ifndef VINECOPLIB_BICOP_H_INCLUDED__
#define VINECOPLIB_BICOP_H_INCLUDED__

#include <Eigen/Dense>
#include <random>
#include <nlopt.hpp>
typedef Eigen::VectorXd VecXd;
typedef Eigen::MatrixXd MatXd;

class Bicop
{

public:

    // hfunctions: the conditioning variable is put second
    virtual VecXd hfunc1(const MatXd &u) = 0;
    virtual VecXd hfunc2(const MatXd &u) = 0;

    // PDF
    virtual VecXd pdf(const MatXd &u) = 0;

    // Inverse of the h-functions
    VecXd hinv1(const MatXd &u);
    VecXd hinv2(const MatXd &u);

    // Simulation
    MatXd simulate(int n);


    // fit statistics
    //virtual double loglik(MatXd u) = 0;

};

#endif // __BICOP_H_INCLUDED__
