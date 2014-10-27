/*
Big list of constants.
We use 'const double' here so the values will be inlined at compilation.
*/

#ifndef __cs_constants_h__
#define __cs_constants_h__

// Define pi so we don't have to worry about missing M_PI (C99).
// pi and multiples in double precision.
const double pi = 3.141592653589793;

// Speed of light, in cm/s.
const double c_cgs = 29979245800.0;
const double c_kms = 299792.458;

// Planck constant, in erg s.
const double h_cgs = 6.6260755e-27;
const double h_bar_cgs = 1.05457266e-27;

// Gravitational constant, in cm^3 g^-1 s^-2.
const double G_cgs = 6.67259e-8;

// Electron charge, in esu.
const double e_cgs = 4.8032068e-10;

// Electron mass, in g.
const double m_e_cgs = 9.1093897e-28;

// Proton mass, in g.
const double m_p_cgs = 1.672623e-24;

// Neutron mass, in g.
const double m_n_cgs = 1.6749286e-24;

// Hydrogen mass, in g.
const double m_H_cgs = 1.6733e-24;

// Atomic mass unit, in g.
const double amu_cgs = 1.6605402e-24;

// Boltzmann constant, in erg/K.
const double k_cgs = 1.38064e-16;

// Electron volt, in erg.
const double eV_cgs = 1.6021772e-12;

// Radiation density constant, in erg cm^-3 K^-4.
const double a_cgs = 7.5646e-15;

// Stefan-Boltzmann constant, in erg cm^-2 K^-4.
const double sigma_sb_cgs = 5.67051e-5;

// Fine structure constant.
const double alpha_cgs = 7.29735308e-3;

// Rydberg constant, in erg.
const double ryd_cgs = 2.1798741e-11;
const double rinf_cgs = 109737.31568538999;

//
// Astro units
//

// Astonomical unit, in cm.
const double AU_cgs = 1.496e13;

// Parsec, in cm.
const double pc_cgs  = 3.08568e18;
const double kpc_cgs = 3.08568e21;
const double Mpc_cgs = 3.08568e24;

// Solar mass, in g.
const double Msun_cgs = 1.99e33;

// Solar radius, in cm.
const double Rsun_cgs = 6.96e10;

// Solar luminosity, in erg/s.
const double Lsun_cgs = 3.9e33;

//
// Cosmological constants
//

// Hubble parameter with h = 1, 100 km/s/Mpc, in 1/s.
const double H0_100_cgs = 3.2407767493712897e-18;
const double H0_100_yr = 1.0220113556817299e-10;

// Critical density with h = 1, in g/cm^3.
// rho_crit = 3 H_0^2 / (8 pi G)
const double rho_crit_100_cgs = 1.8788200386793017e-29;
const double rho_crit_100_Msun_Mpc = 277386144370.4126;

#endif
