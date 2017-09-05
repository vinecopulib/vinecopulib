// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/misc/tools_os.hpp>
#include <stan/command/stanc_helper.hpp>
#include <stan/optimization/bfgs.hpp>
#include "rosenbrock.hpp"

using namespace vinecopulib;

namespace test_stan {

    typedef rosenbrock_namespace::rosenbrock Model;
    typedef stan::optimization::BFGSLineSearch<Model,stan::optimization::BFGSUpdate_HInv<> > Optimizer;

    TEST(test_stan, stanc_works) {
        std::stringstream out;
        std::stringstream err;
        int argc = 4;
        std::vector<const char*> argv_vec;
        argv_vec.push_back("main");
        argv_vec.push_back("--name=rosenbrock");
        argv_vec.push_back("--o=bin/rosenbrock.hpp");
        argv_vec.push_back("bin/rosenbrock.stan");
        const char** argv = &argv_vec[0];
        int rc = stanc_helper(argc, argv, &out, &err);
        EXPECT_TRUE(rc == 0) << "error=" << err.str() << std::endl;
    }

    TEST(test_stan, rosenbrock_bfgs_convergence) {
            // -1,1 is the standard initialization for the Rosenbrock function
            std::vector<double> cont_vector(2);
            cont_vector[0] = -1; cont_vector[1] = 1;
            std::vector<int> disc_vector;

            static const std::string DATA("");
            std::stringstream data_stream(DATA);
            stan::io::dump dummy_context(data_stream);

            Model rb_model(dummy_context);
            std::stringstream out;
            Optimizer bfgs(rb_model, cont_vector, disc_vector, &out);
            EXPECT_EQ("", out.str());

            int ret = 0;
            while (ret == 0) {
                    ret = bfgs.step();
            }
            bfgs.params_r(cont_vector);

            // Check that the return code is normal
            EXPECT_GE(ret,0);

            // Check the correct minimum was found
            EXPECT_NEAR(cont_vector[0],1.0,1e-6);
            EXPECT_NEAR(cont_vector[1],1.0,1e-6);

            // Check that it didn't take too long to get there
            EXPECT_LE(bfgs.iter_num(), 35);
            EXPECT_LE(bfgs.grad_evals(), 70);
    }

    TEST(test_stan, rosenbrock_bfgs_termconds) {
            // -1,1 is the standard initialization for the Rosenbrock function
            std::vector<double> cont_vector(2);
            cont_vector[0] = -1; cont_vector[1] = 1;
            std::vector<int> disc_vector;

            static const std::string DATA("");
            std::stringstream data_stream(DATA);
            stan::io::dump dummy_context(data_stream);

            Model rb_model(dummy_context);

            std::stringstream out;
            Optimizer bfgs(rb_model, cont_vector, disc_vector, &out);
            EXPECT_EQ("", out.str());
            int ret;

            bfgs._conv_opts.maxIts = 1e9;
            bfgs._conv_opts.tolAbsX = 0;
            bfgs._conv_opts.tolAbsF = 0;
            bfgs._conv_opts.tolRelF = 0;
            bfgs._conv_opts.tolAbsGrad = 0;
            bfgs._conv_opts.tolRelGrad = 0;

            bfgs._conv_opts.maxIts = 5;
            bfgs.initialize(cont_vector);
            while (0 == (ret = bfgs.step()));
            EXPECT_EQ(ret,stan::optimization::TERM_MAXIT);
            EXPECT_EQ(bfgs.iter_num(),bfgs._conv_opts.maxIts);
            bfgs._conv_opts.maxIts = 1e9;

            bfgs._conv_opts.tolAbsX = 1e-8;
            bfgs.initialize(cont_vector);
            while (0 == (ret = bfgs.step()));
            EXPECT_EQ(ret,stan::optimization::TERM_ABSX);
            bfgs._conv_opts.tolAbsX = 0;

            bfgs._conv_opts.tolAbsF = 1e-12;
            bfgs.initialize(cont_vector);
            while (0 == (ret = bfgs.step()));
            EXPECT_EQ(ret,stan::optimization::TERM_ABSF);
            bfgs._conv_opts.tolAbsF = 0;

            bfgs._conv_opts.tolRelF = 1e+4;
            bfgs.initialize(cont_vector);
            while (0 == (ret = bfgs.step()));
            EXPECT_EQ(ret,stan::optimization::TERM_RELF);
            bfgs._conv_opts.tolRelF = 0;

            bfgs._conv_opts.tolAbsGrad = 1e-8;
            bfgs.initialize(cont_vector);
            while (0 == (ret = bfgs.step()));
            EXPECT_EQ(ret,stan::optimization::TERM_ABSGRAD);
            bfgs._conv_opts.tolAbsGrad = 0;

            bfgs._conv_opts.tolRelGrad = 1e+3;
            bfgs.initialize(cont_vector);
            while (0 == (ret = bfgs.step()));
            EXPECT_EQ(ret,stan::optimization::TERM_RELGRAD);
            bfgs._conv_opts.tolRelGrad = 0;
    }

    TEST(test_stan, rosenbrock_lbfgs_convergence) {
            // -1,1 is the standard initialization for the Rosenbrock function
            std::vector<double> cont_vector(2);
            cont_vector[0] = -1; cont_vector[1] = 1;
            std::vector<int> disc_vector;

            static const std::string DATA("");
            std::stringstream data_stream(DATA);
            stan::io::dump dummy_context(data_stream);

            Model rb_model(dummy_context);
            std::stringstream out;
            Optimizer bfgs(rb_model, cont_vector, disc_vector, &out);
            EXPECT_EQ("", out.str());

            int ret = 0;
            while (ret == 0) {
                    ret = bfgs.step();
            }
            bfgs.params_r(cont_vector);

            // Check that the return code is normal
            EXPECT_GE(ret,0);

            // Check the correct minimum was found
            EXPECT_NEAR(cont_vector[0],1.0,1e-6);
            EXPECT_NEAR(cont_vector[1],1.0,1e-6);

            // Check that it didn't take too long to get there
            EXPECT_LE(bfgs.iter_num(), 35);
            EXPECT_LE(bfgs.grad_evals(), 70);
    }

    TEST(test_stan, rosenbrock_lbfgs_termconds) {
            // -1,1 is the standard initialization for the Rosenbrock function
            std::vector<double> cont_vector(2);
            cont_vector[0] = -1; cont_vector[1] = 1;
            std::vector<int> disc_vector;

            static const std::string DATA("");
            std::stringstream data_stream(DATA);
            stan::io::dump dummy_context(data_stream);

            Model rb_model(dummy_context);
            std::stringstream out;
            Optimizer bfgs(rb_model, cont_vector, disc_vector, &out);
            EXPECT_EQ("", out.str());
            int ret;

            bfgs._conv_opts.maxIts = 1e9;
            bfgs._conv_opts.tolAbsX = 0;
            bfgs._conv_opts.tolAbsF = 0;
            bfgs._conv_opts.tolRelF = 0;
            bfgs._conv_opts.tolAbsGrad = 0;
            bfgs._conv_opts.tolRelGrad = 0;

            bfgs._conv_opts.maxIts = 5;
            bfgs.initialize(cont_vector);
            while (0 == (ret = bfgs.step()));
            EXPECT_EQ(ret,stan::optimization::TERM_MAXIT);
            EXPECT_EQ(bfgs.iter_num(),bfgs._conv_opts.maxIts);
            bfgs._conv_opts.maxIts = 1e9;

            bfgs._conv_opts.tolAbsX = 1e-8;
            bfgs.initialize(cont_vector);
            while (0 == (ret = bfgs.step()));
            EXPECT_EQ(ret,stan::optimization::TERM_ABSX);
            bfgs._conv_opts.tolAbsX = 0;

            bfgs._conv_opts.tolAbsF = 1e-12;
            bfgs.initialize(cont_vector);
            while (0 == (ret = bfgs.step()));
            EXPECT_EQ(ret,stan::optimization::TERM_ABSF);
            bfgs._conv_opts.tolAbsF = 0;

            bfgs._conv_opts.tolRelF = 1e+4;
            bfgs.initialize(cont_vector);
            while (0 == (ret = bfgs.step()));
            EXPECT_EQ(ret,stan::optimization::TERM_RELF);
            bfgs._conv_opts.tolRelF = 0;

            bfgs._conv_opts.tolAbsGrad = 1e-8;
            bfgs.initialize(cont_vector);
            while (0 == (ret = bfgs.step()));
            EXPECT_EQ(ret,stan::optimization::TERM_ABSGRAD);
            bfgs._conv_opts.tolAbsGrad = 0;

            bfgs._conv_opts.tolRelGrad = 1e+3;
            bfgs.initialize(cont_vector);
            while (0 == (ret = bfgs.step()));
            EXPECT_EQ(ret,stan::optimization::TERM_RELGRAD);
            bfgs._conv_opts.tolRelGrad = 0;
    }

    TEST(test_stan, ConvergenceOptions) {
            stan::optimization::ConvergenceOptions<> a;

            EXPECT_FLOAT_EQ(a.maxIts, 10000);
            EXPECT_FLOAT_EQ(a.fScale, 1);
            EXPECT_FLOAT_EQ(a.tolAbsX, 1e-8);
            EXPECT_FLOAT_EQ(a.tolAbsF, 1e-12);
            EXPECT_FLOAT_EQ(a.tolAbsGrad, 1e-8);
            EXPECT_FLOAT_EQ(a.tolRelF, 1e+4);
            EXPECT_FLOAT_EQ(a.tolRelGrad, 1e+3);
    }

    TEST(test_stan, LsOptions) {
            stan::optimization::LSOptions<> a;

            EXPECT_FLOAT_EQ(a.c1, 1e-4);
            EXPECT_FLOAT_EQ(a.c2, 0.9);
            EXPECT_FLOAT_EQ(a.minAlpha, 1e-12);
            EXPECT_FLOAT_EQ(a.alpha0, 1e-3);
    }

    TEST(test_stan, ModelAdaptor) {
            Eigen::Matrix<double,Eigen::Dynamic,1> cont_vector(2);
            cont_vector[0] = -1; cont_vector[1] = 1;
            std::vector<int> disc_vector;

            static const std::string DATA("");
            std::stringstream data_stream(DATA);
            stan::io::dump dummy_context(data_stream);
            Model rb_model(dummy_context);
            std::stringstream out;
            stan::optimization::ModelAdaptor<Model> mod(rb_model, disc_vector, &out);
            EXPECT_EQ("", out.str());

            // test streams
            EXPECT_NO_THROW(stan::optimization::ModelAdaptor<Model> mod(rb_model, disc_vector, 0));
            EXPECT_NO_THROW(stan::optimization::ModelAdaptor<Model> mod(rb_model, disc_vector, &out));
    }

    TEST(test_stan, ModelAdaptor_fevals) {
            Eigen::Matrix<double,Eigen::Dynamic,1> cont_vector(2);
            cont_vector[0] = -1; cont_vector[1] = 1;
            std::vector<int> disc_vector;

            static const std::string DATA("");
            std::stringstream data_stream(DATA);
            stan::io::dump dummy_context(data_stream);
            Model rb_model(dummy_context);
            std::stringstream out;
            stan::optimization::ModelAdaptor<Model> mod(rb_model, disc_vector, &out);
            EXPECT_EQ("", out.str());

            EXPECT_FLOAT_EQ(mod.fevals(), 0);
    }

    TEST(test_stan, ModelAdaptor_operator_parens__matrix_double) {
            Eigen::Matrix<double,Eigen::Dynamic,1> cont_vector(2);
            cont_vector[0] = -1; cont_vector[1] = 1;
            std::vector<int> disc_vector;

            static const std::string DATA("");
            std::stringstream data_stream(DATA);
            stan::io::dump dummy_context(data_stream);
            Model rb_model(dummy_context);
            std::stringstream out;
            stan::optimization::ModelAdaptor<Model> mod(rb_model, disc_vector, &out);
            EXPECT_EQ("", out.str());

            Eigen::Matrix<double,Eigen::Dynamic,1> grad(2);
            grad[0] = 4;
            grad[1] = 0;
            double f;

            EXPECT_FLOAT_EQ(mod(cont_vector,f), 0);
    }

    TEST(test_stan, ModelAdaptor_operator_parens__matrix_double_matrix) {
            Eigen::Matrix<double,Eigen::Dynamic,1> cont_vector(2);
            cont_vector[0] = -1; cont_vector[1] = 1;
            std::vector<int> disc_vector;

            static const std::string DATA("");
            std::stringstream data_stream(DATA);
            stan::io::dump dummy_context(data_stream);
            Model rb_model(dummy_context);
            std::stringstream out;
            stan::optimization::ModelAdaptor<Model> mod(rb_model, disc_vector, &out);
            EXPECT_EQ("", out.str());

            Eigen::Matrix<double,Eigen::Dynamic,1> grad(2);
            grad[0] = 4;
            grad[1] = 0;
            double f;

            EXPECT_FLOAT_EQ(mod(cont_vector,f, grad), 0);
    }

    TEST(test_stan, ModelAdaptor_df) {
            Eigen::Matrix<double,Eigen::Dynamic,1> cont_vector(2);
            cont_vector[0] = -1; cont_vector[1] = 1;
            std::vector<int> disc_vector;

            static const std::string DATA("");
            std::stringstream data_stream(DATA);
            stan::io::dump dummy_context(data_stream);
            Model rb_model(dummy_context);
            std::stringstream out;
            stan::optimization::ModelAdaptor<Model> mod(rb_model, disc_vector, &out);
            EXPECT_EQ("", out.str());

            Eigen::Matrix<double,Eigen::Dynamic,1> grad(2);
            grad[0] = 4;
            grad[1] = 0;

            EXPECT_FLOAT_EQ(mod.df(cont_vector, grad), 0);
    }

    TEST(test_stan, BFGSLineSearch) {
            std::vector<double> cont_vector(2);
            cont_vector[0] = -1; cont_vector[1] = 1;
            std::vector<int> disc_vector;

            static const std::string DATA("");
            std::stringstream data_stream(DATA);
            stan::io::dump dummy_context(data_stream);
            Model rb_model(dummy_context);
            std::stringstream out;
            Optimizer bfgs(rb_model, cont_vector, disc_vector, &out);
            EXPECT_EQ("", out.str());

            // test streams
            EXPECT_NO_THROW(Optimizer bfgs(rb_model, cont_vector, disc_vector, 0));
            EXPECT_NO_THROW(Optimizer bfgs(rb_model, cont_vector, disc_vector, &out));
    }


    TEST(test_stan, BFGSLineSearch_initialize) {
            std::vector<double> cont_vector(2);
            cont_vector[0] = -1; cont_vector[1] = 1;
            std::vector<int> disc_vector;

            static const std::string DATA("");
            std::stringstream data_stream(DATA);
            stan::io::dump dummy_context(data_stream);
            Model rb_model(dummy_context);
            std::stringstream out;
            Optimizer bfgs(rb_model, cont_vector, disc_vector, &out);
            EXPECT_EQ("", out.str());

            bfgs.initialize(cont_vector);
            EXPECT_FLOAT_EQ(bfgs.curr_f(), 4);
            EXPECT_FLOAT_EQ(bfgs.curr_x()[0], -1);
            EXPECT_FLOAT_EQ(bfgs.curr_x()[1], 1);
            EXPECT_FLOAT_EQ(bfgs.curr_g()[0], -4);
            EXPECT_FLOAT_EQ(bfgs.curr_g()[1], 0);
            EXPECT_FLOAT_EQ(bfgs.curr_p()[0], 4);
            EXPECT_FLOAT_EQ(bfgs.curr_p()[1], 0);
    }

    TEST(test_stan, BFGSLineSearch_grad_evals) {
            std::vector<double> cont_vector(2);
            cont_vector[0] = -1; cont_vector[1] = 1;
            std::vector<int> disc_vector;

            static const std::string DATA("");
            std::stringstream data_stream(DATA);
            stan::io::dump dummy_context(data_stream);
            Model rb_model(dummy_context);
            std::stringstream out;
            Optimizer bfgs(rb_model, cont_vector, disc_vector, &out);
            EXPECT_EQ("", out.str());

            EXPECT_FLOAT_EQ(bfgs.grad_evals(), 1);
    }

    TEST(test_stan, BFGSLineSearch_logp) {
            std::vector<double> cont_vector(2);
            cont_vector[0] = -1; cont_vector[1] = 1;
            std::vector<int> disc_vector;

            static const std::string DATA("");
            std::stringstream data_stream(DATA);
            stan::io::dump dummy_context(data_stream);
            Model rb_model(dummy_context);
            std::stringstream out;
            Optimizer bfgs(rb_model, cont_vector, disc_vector, &out);
            EXPECT_EQ("", out.str());

            EXPECT_FLOAT_EQ(bfgs.logp(), -4);
    }

    TEST(test_stan, BFGSLineSearch_grad_norm) {
            std::vector<double> cont_vector(2);
            cont_vector[0] = -1; cont_vector[1] = 1;
            std::vector<int> disc_vector;

            static const std::string DATA("");
            std::stringstream data_stream(DATA);
            stan::io::dump dummy_context(data_stream);
            Model rb_model(dummy_context);
            std::stringstream out;
            Optimizer bfgs(rb_model, cont_vector, disc_vector, &out);
            EXPECT_EQ("", out.str());

            std::vector<double> grad;

            bfgs.grad(grad);
            EXPECT_FLOAT_EQ(bfgs.grad_norm(), 4);
    }

    TEST(test_stan, BFGSLineSearch_grad) {
            std::vector<double> cont_vector(2);
            cont_vector[0] = -1; cont_vector[1] = 1;
            std::vector<int> disc_vector;

            static const std::string DATA("");
            std::stringstream data_stream(DATA);
            stan::io::dump dummy_context(data_stream);
            Model rb_model(dummy_context);
            std::stringstream out;
            Optimizer bfgs(rb_model, cont_vector, disc_vector, &out);
            EXPECT_EQ("", out.str());
            std::vector<double> grad;

            bfgs.grad(grad);
            EXPECT_FLOAT_EQ(grad.size(), 2);
            EXPECT_FLOAT_EQ(grad[0], 4);
            EXPECT_FLOAT_EQ(grad[1], 0);
    }

    TEST(test_stan, BFGSLineSearch_params_r) {
            std::vector<double> cont_vector(2);
            cont_vector[0] = -1; cont_vector[1] = 1;
            std::vector<int> disc_vector;

            static const std::string DATA("");
            std::stringstream data_stream(DATA);
            stan::io::dump dummy_context(data_stream);
            Model rb_model(dummy_context);
            std::stringstream out;
            Optimizer bfgs(rb_model, cont_vector, disc_vector, &out);
            EXPECT_EQ("", out.str());
            std::vector<double> x;

            bfgs.params_r(x);
            EXPECT_FLOAT_EQ(x.size(), 2);
            EXPECT_FLOAT_EQ(x[0], -1);
            EXPECT_FLOAT_EQ(x[1], 1);
    }
}
