/*
Path integration test cases.
- Uniform medium
- Gaussian bump
*/

#include <cstdio>

#include "intergalactic_spec.h"

#include "catch.hpp"


TEST_CASE("uniform medium", "[spec]" ) {
    double odf = 1.0;
    double v_domain = 1e7;
    int n = 16;

    double *n_array = new double[n];
    double *v_array = new double[n];
    double *t_array = new double[n];
    double *tau0 = new double[n];
    double *tau1 = new double[n];
    double *tau2 = new double[n];

    // path data
    double n_0 = 1.0;
    double t_0 = 1e4;
    double b = sqrt(2.0 * k_cgs * t_0 / m_H_cgs);
    double vj = v_domain / n * (n/2 + 0.5);
    double tau_ans = odf * n_0 / 2 * (erf(vj / b) - erf((vj - v_domain) / b));

    SECTION("v = 0") {
        for (int i = 0; i < n; ++i) {
            n_array[i] = n_0;
            v_array[i] = 0.0;
            t_array[i] = t_0;
        }
    }

    SECTION("v = +const") {
        // path data
        for (int i = 0; i < n; ++i) {
            n_array[i] = n_0;
            v_array[i] = 0.5 * v_domain;
            t_array[i] = t_0;
        }
    }

    optical_depth_midpoint(odf, m_H_cgs, v_domain,
        n, n_array, v_array, t_array, n, tau0);
    optical_depth_midpoint_lss(odf, m_H_cgs, v_domain, 2 * n,
        n, n_array, v_array, t_array, n, tau1);
    optical_depth_pwc_exact(odf, m_H_cgs, v_domain,
        n, n_array, v_array, t_array, n, tau2);

    const double abs_tol = 1e-10;
    const double rel_tol = 1e-6;
    for (int i = 0; i < n; ++i) {
        REQUIRE( fabs(tau0[i] - tau_ans) <= abs_tol + rel_tol * tau_ans );
        REQUIRE( fabs(tau1[i] - tau_ans) <= abs_tol + rel_tol * tau_ans );
        REQUIRE( fabs(tau2[i] - tau_ans) <= abs_tol + rel_tol * tau_ans );
    }

    delete [] n_array;
    delete [] v_array;
    delete [] t_array;
    delete [] tau0;
    delete [] tau1;
    delete [] tau2;
}

TEST_CASE("gaussian bump", "[spec]") {
    double odf = 1.0;
    double v_domain = 1e7;
    int n = 256;

    double *n_array = new double[n];
    double *v_array = new double[n];
    double *t_array = new double[n];
    double *tau0 = new double[n];
    double *tau1 = new double[n];
    double *tau2 = new double[n];

    // path data
    double n_0 = 1.0;
    double v_0 = 0.5 * v_domain;
    double a = 0.1 * v_domain;
    for (int i = 0; i < n; ++i) {
        double x = ((double)i + 0.5) / n;
        double dx = (x * v_domain - v_0) / a;
        n_array[i] = n_0 * exp(-dx*dx);
        v_array[i] = 0.0;
        t_array[i] = 1e4;
    }

    optical_depth_midpoint(odf, m_H_cgs, v_domain,
        n, n_array, v_array, t_array, n, tau0);
    optical_depth_midpoint_lss(odf, m_H_cgs, v_domain, n,
        n, n_array, v_array, t_array, n, tau1);
    optical_depth_pwc_exact(odf, m_H_cgs, v_domain,
        n, n_array, v_array, t_array, n, tau2);

    double *tau_ans = new double[n];
    double b = sqrt(2.0 * k_cgs * t_array[0] / m_H_cgs);
    double c = sqrt(a*a + b*b);
    for (int i = 0; i < n; ++i) {
        double v = v_domain / n * (i + 0.5);

        tau_ans[i] = odf * n_0 * a / (2 * c) * exp( -pow(v - v_0, 2) / (c * c) )
            * ( erf( ( a*a * (v_domain - v) + b*b * (v_domain - v_0))
                     / (a * b * c) )
                - erf( (a*a * (-v) + b*b * (-v_0)) / (a * b * c) ) );

        // DEBUG
        /*
        double dx = fabs(v - v_0) / b;
        if (dx < 0.3) {
            printf("x %f   err0 %f   err1  %f   err2 %f\n",
               dx, tau0[i] / tau_ans[i] - 1, tau1[i] / tau_ans[i] - 1,
               tau2[i] / tau_ans[i] - 1);
        }
        */
    }

    // note: we should only check away from boundaries since the solution
    // is assuming it's isolated and the routines are periodic.
    const double abs_tol = 1e-10;
    const double rel_tol = 1e-2;

    for (int i = 0; i < n; ++i) {
        double x = ((double)i + 0.5) / n;
        double dx = (x * v_domain - v_0) / a;
        if (fabs(dx) < 2) {
            double t = tau_ans[i];
            REQUIRE( fabs(tau0[i] - t)
                     <= abs_tol + rel_tol * t );
            REQUIRE( fabs(tau1[i] - t)
                     <= abs_tol + rel_tol * t );
            REQUIRE( fabs(tau2[i] - t)
                     <= abs_tol + rel_tol * t );
        }
    }

    delete [] n_array;
    delete [] v_array;
    delete [] t_array;
    delete [] tau0;
    delete [] tau1;
    delete [] tau2;
}
