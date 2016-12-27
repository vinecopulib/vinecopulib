/*
    Copyright 2016 Thibault Vatter

    This file is part of vinecoplib.

    vinecoplib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    vinecoplib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with vinecoplib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "include/r_instance.hpp"

RInstance::RInstance() {
    // Default sample size
    int n = 1e3;
    n_ = n;
    family_ = 0;
    parameters_ = 4 * VecXd::Ones(2);
    U_ = MatXd::Zero(n, 2);
    tau_ = 0.5;

    // Set-up some R commands
    std::string u1 = "u1 <- runif()";
    u1.insert(12, std::to_string(n_));
    std::string u2 = "u2 <- runif()";
    u2.insert(12, std::to_string(n_));
    std::string eval_fct = "par <- BiCopTau2Par(1,)";
    eval_fct.insert(22, std::to_string(tau_));
    std::string eval_sim = "U <- BiCopSim(,1,par)";
    eval_sim.insert(14, std::to_string(n_));
    std::string load_VineCopula = "library(VineCopula)";

    // Declare the RInside instance and load the VineCopula package
    R.parseEval(load_VineCopula);

    // Evaluate the tests random matrix U and load it in C++
    // const Eigen::Map<Eigen::VectorXd> U2 = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(R.parseEval(u2));
    R.parseEval(eval_fct);
    R.parseEval(eval_sim);
    R.parseEval("u1 <- U[,1]");
    R.parseEval("u2 <- U[,2]");
    const Eigen::Map<VecXd> U1 = Rcpp::as<Eigen::Map<VecXd>>(R.parseEval(u1));
    const Eigen::Map<VecXd> U2 = Rcpp::as<Eigen::Map<VecXd>>(R.parseEval(u2));
    U_.col(0) = U1;
    U_.col(1) = U2;
}

VecXd RInstance::eval_in_R(std::string eval_fct, int start)
{

    std::string new_eval_fct = eval_fct;
    std::string par = std::to_string(parameters_(0));
    int family = get_family();

    // take care of the rotations
    int rotated_families_array[] = {3,4,6,7,8,9,10};
    std::vector<int> rotated_families;
    rotated_families.assign(rotated_families_array,rotated_families_array+7);
    if(std::find(rotated_families.begin(), rotated_families.end(), family) != rotated_families.end()) {
        int rotation = get_rotation();
        if (rotation == 180)
            family += 10;
        if (rotation == 90)
            family += 20;
        if (rotation == 270)
            family += 30;
    }

    // evaluate the function in R
    std::string fam = std::to_string(family);
    new_eval_fct.insert(start,fam);
    new_eval_fct.insert(start+fam.length()+1,par);
    if (parameters_.size() == 1) {
        new_eval_fct.insert(start+par.length()+fam.length()+2,std::to_string(0));
    } else {
        new_eval_fct.insert(start+par.length()+fam.length()+2,std::to_string(parameters_(1)));
    }
    const Eigen::Map<Eigen::VectorXd>  f = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(R.parseEval(new_eval_fct));
    VecXd result = f;
    return(result);
}

// getters --------------------------------

int RInstance::get_family()
{
    return this->family_;
}

VecXd RInstance::get_parameters()
{
    return this->parameters_;
}

int RInstance::get_rotation()
{
    return this->rotation_;
}

double RInstance::get_tau()
{
    return this->tau_;
}

// setters --------------------------------

void RInstance::set_family(int family)
{
    family_ = family;
}

void RInstance::set_parameters(VecXd parameters)
{
    parameters_ = parameters;
}

void RInstance::set_rotation(int rotation)
{
    rotation_ = rotation;
}

void RInstance::set_tau(double tau)
{
    tau_ = tau;
}

RInside RInstance::get_R()
{
    return this->R;
}
MatXd RInstance::get_U()
{
    return this->U_;
}
void RInstance::change_n(int n)
{
    this->n_ = n;

    // Set-up some R commands
    std::string u1 = "u1 <- runif()";
    u1.insert(12, std::to_string(n));
    std::string u2 = "u2 <- runif()";
    u2.insert(12, std::to_string(n));

    // Evaluate the tests random vectors u1 and u2 and load them in C++
    const Eigen::Map<Eigen::VectorXd>  U1 = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(R.parseEval(u1));
    const Eigen::Map<Eigen::VectorXd>  U2 = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(R.parseEval(u2));
    this->U_ = MatXd::Zero(n,2);
    this->U_.col(0) = U1;
    this->U_.col(1) = U2;
}
