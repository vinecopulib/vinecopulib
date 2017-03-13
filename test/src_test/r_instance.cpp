/*
* The MIT License (MIT)
*
* Copyright © 2017 Thibault Vatter and Thomas Nagler
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "include/r_instance.hpp"

RInstance::RInstance() {
    // Default sample size
    int n = 2e3;
    n_ = n;
    family_ = 0;
    rotation_ = 0;
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
    
    // Set the seed and load the VineCopula package
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
    int family = get_family();
    double parameter = 0;
    if (family != 0) {
        parameter = parameters_(0);
        // take care of the rotations
        std::vector<int> rotated_families = {3,4,6,7,8,9,10};
        if (is_member(family, rotated_families)) {
            int rotation = get_rotation();
            if (rotation == 90) {
                family += 20;
                parameter *= -1;
            } else if (rotation == 180) {
                family += 10;
            } else if (rotation == 270) {
                family += 30;
                parameter *= -1;
            }
        }
    }
        
    // evaluate the function in R
    std::string fam = std::to_string(family);
    std::string par = std::to_string(parameter);
    new_eval_fct.insert(start,fam);
    new_eval_fct.insert(start+fam.length()+1,par);
    if (parameters_.size() < 2) {
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

int RInstance::get_n()
{
    return this->n_;
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
