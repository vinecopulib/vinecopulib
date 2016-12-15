//
// Created by Vatter Thibault on 13/12/16.
//

#include "include/bicop.h"

VecXd Bicop::hinv1(const MatXd &u)
{
    MatXd v = u;
    VecXd u1 = u.col(1);

    int j=0, br = 0, in = 1, it = 0;
    double tol = 1e-12, xl = 1e-10, xh = 1-1e-10;

    v.col(1) = xl*VecXd::Ones(u1.size());
    VecXd fl = (hfunc1(v)-u1).cwiseAbs();
    v.col(1) = xh*VecXd::Ones(u1.size());
    VecXd fh = (hfunc1(v)-u1).cwiseAbs();
    VecXd fm = VecXd::Ones(u1.size());

    for (j=0;j<u.rows();j++)
    {

        if (fl(j) <= tol) {
            v(j,1) = xl;
            br = 1;
        }
        if (fh(j) <= tol) {
            v(j,1) = xh;
            br = 1;
        }

        while (!br){
            v(j,1) = (xh + xl) / 2.0;
            fm(j) = hfunc1(v.row(j))(0)-u1(j);

            //stop if values become too close (avoid infinite loop)
            if (fabs(fm(j)) <= tol) br = 1;
            if (fabs(xl-xh) <= tol) br = 1;

            if (fm(j) < 0.0) {
                xl = v(j,1);
                fh(j) = fm(j);
            } else {
                xh = v(j,1);
                fl(j) = fm(j);
            }

            //stop if too many iterations are required (avoid infinite loop)
            ++it;
            if (it > 50) br = 1;
        }
        br = 0;
        in = 1;
        it = 0;
        xl = 1e-10;
        xh = 1-1e-10;
    }

    return(v.col(1));
}


VecXd Bicop::hinv2(const MatXd &u)
{
    MatXd v = u;
    VecXd u1 = u.col(0);

    int j=0, br = 0, in = 1, it = 0;
    double tol = 1e-12, xl = 1e-10, xh = 1-1e-10;

    v.col(0) = xl*VecXd::Ones(u1.size());
    VecXd fl = (hfunc2(v)-u1).cwiseAbs();
    v.col(0) = xh*VecXd::Ones(u1.size());
    VecXd fh = (hfunc2(v)-u1).cwiseAbs();
    VecXd fm = VecXd::Ones(u1.size());

    for (j=0;j<u.rows();j++)
    {

        if (fl(j) <= tol) {
            v(j,0) = xl;
            br = 1;
        }
        if (fh(j) <= tol) {
            v(j,0) = xh;
            br = 1;
        }

        while (!br){
            v(j,0) = (xh + xl) / 2.0;
            fm(j) = hfunc2(v.row(j))(0)-u1(j);

            //stop if values become too close (avoid infinite loop)
            if (fabs(fm(j)) <= tol) br = 1;
            if (fabs(xl-xh) <= tol) br = 1;

            if (fm(j) < 0.0) {
                xl = v(j,0);
                fh(j) = fm(j);
            } else {
                xh = v(j,0);
                fl(j) = fm(j);
            }

            //stop if too many iterations are required (avoid infinite loop)
            ++it;
            if (it > 50) br = 1;
        }
        br = 0;
        in = 1;
        it = 0;
        xl = 1e-10;
        xh = 1-1e-10;
    }

    return(v.col(0));
}

MatXd Bicop::simulate(int n)
{
    std::default_random_engine generator(time(0));
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    MatXd U = MatXd::Zero(n,2);
    for (int i=0; i<n; ++i) {
        U(i,0) = distribution(generator);
        U(i,1) = distribution(generator);
    }
    U.col(1) = hinv1(U);

    return(U);
}