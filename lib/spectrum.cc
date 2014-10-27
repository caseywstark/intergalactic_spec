/*
Copyright 2014, Casey W. Stark. See LICENSE for more information.
*/

#include <cstdio>
#include <cmath>

#include "util.h"

#include "spectrum.h"

double
odf_cosmo(const double h, const double e_z)
{
    double H_z = H0_100_cgs * h * e_z;
    return 1.0 / H_z;
}

double
od_prefactor_hilya(const double h, const double e_z)
{
    return odf_const * odf_hilya * odf_cosmo(h, e_z);
}

void
optical_depth_midpoint(const double od_factor, const double m_x,
    const double v_domain,
    const int num_elements, const double * const n_array,
    const double * const vpara_array, const double * const t_array,
    const int num_pixels, double * const tau_array)
{
    // numerical floors
    const double t_min = 1.0;
    const double vth_min = 1e4;   // cm/s

    const double element_dv = v_domain / num_elements;
    const double pixel_dv = v_domain / num_pixels;

    // The thermal velocity prefactor.
    // v_th = sqrt( 2 k T / m ), so we compute sqrt( 2 k / m_x ) once.
    const double vth_factor = sqrt(2.0 * k_cgs / m_x);

    // Reset tau array to 0.0.
    for (int i = 0; i < num_pixels; ++i) { tau_array[i] = 0.0; }

    // Loop over elements, adding optical depth to nearby pixels
    for (int i = 0; i < num_elements; ++i) {
        // make sure this cell has positive density.
        if (n_array[i] > 0.0) {
            // thermal velocity.
            double v_doppler = vth_factor * sqrt(t_array[i] + t_min) + vth_min;
            // midpoint velocity
            double v_i = element_dv * (i + 0.5);
            // line-center velocity
            double v_lc = v_i + vpara_array[i];

            // The line center and thermal broadening in pixel index space.
            double pix_lc = v_lc / pixel_dv;
            double pix_doppler = v_doppler / pixel_dv;

            // Figure out which pixels we should restrict to.
            int num_integ_pixels = (int)(5 * pix_doppler + 0.5) + 1;
            int ipix_lo = pix_lc - num_integ_pixels - 1;
            int ipix_hi = pix_lc + num_integ_pixels + 1;

            // Note the <= to include the hi bin.
            for (int j = ipix_lo; j <= ipix_hi; ++j) {
                // Calculate the bin velocity with the original index to match v_lc.
                double v_pixel = pixel_dv * (j + 0.5);
                // Wrap j into bounds.
                int jw = index_wrap(j, num_pixels);

                // Add tau contribution to the pixel.
                double dx = (v_pixel - v_lc) / v_doppler;
                tau_array[jw] += n_array[i] / v_doppler * exp(-dx*dx);
            }
        }
    }

    // Don't forget prefactor.
    double f = od_factor * element_dv / sqrt(pi);
    for (int i = 0; i < num_pixels; ++i) {
        tau_array[i] *= f;
    }
}

void
optical_depth_midpoint_lss(const double od_factor, const double m_x,
    const double v_domain, const int num_samples,
    const int num_elements, const double * const n_array,
    const double * const v_array, const double * const t_array,
    const int num_pixels, double * const tau_array)
{
    // numerical floors
    const double t_min = 1.0;
    const double vth_min = 1e4;   // cm/s

    const double sample_dv = v_domain / num_samples;
    const double pixel_dv = v_domain / num_pixels;

    const double vth_factor = sqrt(2.0 * k_cgs / m_x);

    // Reset tau array to 0.0.
    for (int i = 0; i < num_pixels; ++i) { tau_array[i] = 0.0; }

    // Loop over elements, adding optical depth to nearby pixels
    for (int i = 0; i < num_samples; ++i) {
        // just do the linterp inline.
        double x = ((double)i + 0.5) / num_samples;
        // element index coords.
        double ex = x * num_elements - 0.5;
        int e0 = floor(ex);
        int e1 = e0 + 1;
        double w0 = e1 - ex;
        double w1 = 1.0 - w0;
        // fix for periodicity.
        if (e0 == -1) { e0 = num_elements - 1; }
        if (e1 == num_elements) { e1 = 0; }

        double n_i = w0 * n_array[e0] + w1 * n_array[e1];
        double v_i = w0 * v_array[e0] + w1 * v_array[e1];
        double t_i = w0 * t_array[e0] + w1 * t_array[e1];

        // make sure this cell has positive density.
        if (n_i > 0.0) {
            // thermal velocity.
            double v_doppler = vth_factor * sqrt(t_i + t_min) + vth_min;
            // line-center velocity
            double vlc_m = sample_dv * (i + 0.5) + v_i;

            // The line center and thermal broadening in pixel index space.
            double pix_lc = vlc_m / pixel_dv;
            double pix_doppler = v_doppler / pixel_dv;

            // Figure out which pixels we should restrict to.
            int num_integ_pixels = (int)(5 * pix_doppler + 0.5) + 1;
            int ipix_lo = pix_lc - num_integ_pixels - 1;
            int ipix_hi = pix_lc + num_integ_pixels + 1;

            // Note the <= to include the hi bin.
            for (int j = ipix_lo; j <= ipix_hi; ++j) {
                // Calculate the bin velocity with the original index to match v_lc.
                double v_pixel = pixel_dv * (j + 0.5);
                // Wrap j into bounds.
                int jw = index_wrap(j, num_pixels);

                // Add tau contribution to the pixel.
                double dx = (v_pixel - vlc_m) / v_doppler;
                tau_array[jw] += n_i / v_doppler * exp(-dx*dx);
            }
        }
    }

    // Don't forget prefactor.
    double f = od_factor * sample_dv / sqrt(pi);
    for (int i = 0; i < num_pixels; ++i) {
        tau_array[i] *= f;
    }
}



void
optical_depth_pwc_exact(const double od_factor, const double m_x,
    const double v_domain,
    const int num_elements, const double * const n_array,
    const double * const vpara_array, const double * const t_array,
    const int num_pixels, double * const tau_array)
{
    // numerical floors
    const double t_min = 1.0;
    const double vth_min = 1e4;   // cm/s

    const double element_dv = v_domain / num_elements;
    const double pixel_dv = v_domain / num_pixels;

    const double vth_factor = sqrt(2.0 * k_cgs / m_x);

    // Reset tau array to 0.0.
    for (int i = 0; i < num_pixels; ++i) { tau_array[i] = 0.0; }

    // Loop over elements, adding optical depth to nearby pixels
    for (int i = 0; i < num_elements; ++i) {
        // make sure this cell has positive density.
        if (n_array[i] > 0.0) {
            // thermal velocity.
            double v_doppler = vth_factor * sqrt(t_array[i] + t_min) + vth_min;
            // line-center velocity
            double vlc_l = element_dv * i + vpara_array[i];
            double vlc_m = element_dv * (i + 0.5) + vpara_array[i];
            double vlc_h = element_dv * (i + 1.0) + vpara_array[i];

            // The line center and thermal broadening in pixel index space.
            double pix_lc = vlc_m / pixel_dv;
            double pix_doppler = v_doppler / pixel_dv;

            // Figure out which pixels we should restrict to.
            int num_integ_pixels = (int)(5 * pix_doppler + 0.5) + 1;
            int ipix_lo = pix_lc - num_integ_pixels - 1;
            int ipix_hi = pix_lc + num_integ_pixels + 1;

            // Note the <= to include the hi bin.
            for (int j = ipix_lo; j <= ipix_hi; ++j) {
                // Calculate the bin velocity with the original index to match v_lc.
                double v_pixel = pixel_dv * (j + 0.5);
                // Wrap j into bounds.
                int jw = index_wrap(j, num_pixels);

                // Add tau contribution to the pixel.
                double dxl = (v_pixel - vlc_l) / v_doppler;
                double dxh = (v_pixel - vlc_h) / v_doppler;
                tau_array[jw] += n_array[i] * ( erf(dxl) - erf(dxh) );
            }
        }
    }

    // Don't forget prefactor.
    double f = 0.5 * od_factor;
    for (int i = 0; i < num_pixels; ++i) {
        tau_array[i] *= f;
    }
}
