/*
Run a realistic problem with all methods.
*/

#include <cstdio>
#include <cmath>

#include <string>

#include "cosmo_spec.h"

#include "catch.hpp"

const std::string test_data_path = "test/sample_skewer.bin";

// cosmology
const double omega_m = 0.275;
const double h = 0.7;

// snapshot
const double z = 3.0;
const double e_z = sqrt(omega_m * pow(1.0 + z, 3) + (1.0 - omega_m));
const double da_dt = 100.0 * e_z / (1.0 + z);  // km/s h/Mpc

// box size
const double l = 20.0;  // Mpc/h
const double v_domain = 1.0e5 * da_dt * l;  // km/s -> cm/s

// resolution
const int n = 1024;

void
read_test_data(const std::string path, const int n, double * const n_array,
    double * const v_array, double * const t_array, double * const tau_array)
{
    FILE *f = fopen(path.c_str(), "r");
    fread(n_array, sizeof(double), n, f);
    fread(v_array, sizeof(double), n, f);
    fread(t_array, sizeof(double), n, f);
    fread(tau_array, sizeof(double), n, f);
    fclose(f);
}

TEST_CASE("regression on typical data", "[spec]") {
    double *n_array = new double[n];
    double *v_array = new double[n];
    double *t_array = new double[n];

    double *tau_ref = new double[n];
    double *tau0 = new double[n];
    double *tau1 = new double[n];
    double *tau2 = new double[n];

    read_test_data(test_data_path, n, n_array, v_array, t_array, tau_ref);

    double odf = od_prefactor_hilya(h, e_z);

    optical_depth_midpoint(odf, m_H_cgs, v_domain,
        n, n_array, v_array, t_array, n, tau0);
    optical_depth_midpoint_lss(odf, m_H_cgs, v_domain, 2 * n,
        n, n_array, v_array, t_array, n, tau1);
    optical_depth_pwc_exact(odf, m_H_cgs, v_domain,
        n, n_array, v_array, t_array, n, tau2);


    // DEBUG
    /*
    FILE *fs = fopen("reg_results.txt", "w");

    for (int i = 0; i < n; ++i) {
        double r = tau_ref[i];
        printf("i %i   tr %f   t0 %f   t1 %f   t2 %f\n",
            i, r, 100.0 * (tau0[i] / r - 1), 100.0 * (tau1[i] / r - 1),
            100.0 * (tau2[i] / r - 1));

        fprintf(fs, "%e %e %e %e\n", r, tau0[i], tau1[i], tau2[i]);
    }

    fclose(fs);
    */

    double abs_tol = 1e-10;
    double rel_tol = 3e-2; // 3 percent

    for (int i = 0; i < n; ++i) {
        REQUIRE( fabs(tau0[i] - tau_ref[i]) <= abs_tol + rel_tol * tau_ref[i] );
    }

    for (int i = 0; i < n; ++i) {
        REQUIRE( fabs(tau1[i] - tau_ref[i]) <= abs_tol + rel_tol * tau_ref[i] );
    }

    for (int i = 0; i < n; ++i) {
        REQUIRE( fabs(tau2[i] - tau_ref[i]) <= abs_tol + rel_tol * tau_ref[i] );
    }

    delete [] n_array;
    delete [] v_array;
    delete [] t_array;

    delete [] tau_ref;
    delete [] tau0;
    delete [] tau1;
    delete [] tau2;
}


// also test that the mp and mpss methods are the same
TEST_CASE("subsample test", "[spec]") {
    double *n_array = new double[n];
    double *v_array = new double[n];
    double *t_array = new double[n];

    double *tau_mp = new double[n];
    double *tau_ss1 = new double[n];
    double *tau_ss2 = new double[n];

    read_test_data(test_data_path, n, n_array, v_array, t_array, tau_mp);

    double odf = od_prefactor_hilya(h, e_z);

    optical_depth_midpoint(odf, m_H_cgs, v_domain,
        n, n_array, v_array, t_array, n, tau_mp);
    optical_depth_midpoint_lss(odf, m_H_cgs, v_domain, 1 * n,
        n, n_array, v_array, t_array, n, tau_ss1);
    optical_depth_midpoint_lss(odf, m_H_cgs, v_domain, 2 * n,
        n, n_array, v_array, t_array, n, tau_ss2);

    double abs_tol = 1e-10;
    double rel_tol = 1e-10;

    for (int i = 0; i < n; ++i) {
        REQUIRE( fabs(tau_ss1[i] - tau_mp[i]) <= abs_tol + rel_tol * tau_mp[i] );
    }

    rel_tol = 0.03;
    for (int i = 0; i < n; ++i) {
        REQUIRE( fabs(tau_ss2[i] - tau_mp[i]) <= abs_tol + rel_tol * tau_mp[i] );
    }

    delete [] n_array;
    delete [] v_array;
    delete [] t_array;

    delete [] tau_mp;
    delete [] tau_ss1;
    delete [] tau_ss2;
}

