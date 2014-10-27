/*
Copyright 2014, Casey W. Stark. See LICENSE for more information.

The spectrum header provides the functions for computing the optical depth.

*/

#ifndef __cs_spectrum_h__
#define __cs_spectrum_h__

#include "constants.h"

const double odf_const = pi * e_cgs * e_cgs / (m_e_cgs * c_cgs);

const double hilya_f = 0.41619671797998276;
const double hilya_lambda = 1.2150227340678867e-05; // in cm
const double odf_hilya = hilya_f * hilya_lambda;

double
odf_cosmo(const double h, const double e_z);

double
od_prefactor_hilya(const double h, const double e_z);

// Midpoint evaluated at the center of elements
void
optical_depth_midpoint(const double prefactor, const double m_x,
    const double v_domain,
    const int num_elements, const double * const n_array,
    const double * const vpara_array, const double * const t_array,
    const int num_pixels, double * const tau_array);

// Midpoint evaluated at evenly spaced points determined by `num_samples`.
// Performs linear interpolation of n, v, and t to subsample.
void
optical_depth_midpoint_lss(const double prefactor, const double m_x,
    const double v_domain, const int num_samples,
    const int num_elements, const double * const n_array,
    const double * const vpara_array, const double * const t_array,
    const int num_pixels, double * const tau_array);

// Analytic form for piecewise constant data.
void
optical_depth_pwc_exact(const double od_factor, const double m_x,
    const double v_domain,
    const int num_elements, const double * const n_array,
    const double * const vpara_array, const double * const t_array,
    const int num_pixels, double * const tau_array);

// Not sure if it would be useful, but we could add a method to subsample using
// the piecewise constant form by filling new arrays with linterp subsampled
// data.

#endif
