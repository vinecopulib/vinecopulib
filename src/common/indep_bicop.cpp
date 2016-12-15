//
// Created by Vatter Thibault on 14/12/16.
//

#include "include/indep_bicop.h"

// constructor
IndepBicop::IndepBicop()
{
    family_ = 0;
    parameters_ = VecXd::Zero(1);
}

// hfunctions: the conditioning variable is put second
VecXd IndepBicop::hfunc1(const MatXd &u)
{
    VecXd v = u.col(1);
    return(v);
}

VecXd IndepBicop::hfunc2(const MatXd &u)
{
    VecXd v = u.col(0);
    return(v);
}

// PDF
VecXd IndepBicop::pdf(const MatXd &u)
{
    return(VecXd::Ones(u.rows()));
}