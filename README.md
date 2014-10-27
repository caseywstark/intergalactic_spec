# Intergalactic Spec

A small C++ library that calculates the optical depth to resonant line
scattering in a cosmological background. This was implemented for Ly$\alpha$
forest calculations, but can easily be extended to other atomic species.

## Building

The build system is just a `Makefile` that creates the static library
`libispec.a`. Please modify the make variables `CXX`, `CPPFLAGS`,
`CXXFLAGS`, and `LDFLAGS` as needed.

## Problem setup

We want to compute the optical depth to resonant line scattering $\tau$, given
the (proper) number density of species X $n_X$, the peculiar velocity
$v_\parallel$, and the temperature $T$ along some path. The optical depth is

\tau_\nu = \int n_X \sigma_\nu dr

where $\sigma_\nu$ is the cross section of the interaction, and $dr$ is the
(proper) path length.

The cross section for resonant line scattering is

\sigma_\nu = \frac{\pi e^2}{m_e c} \frac{f_{lu}}{\Delta \nu_D} \phi_\nu \\
    = \frac{\pi e^2}{m_e c} \lambda_{lu} f_{lu} \frac{\phi_\nu}{b}

where $\lambda_{lu}$ is the transition wavelength, $f_{lu}$ is the transition
oscillator strength, $\phi_\nu$ is the line profile, and $b$ is the Doppler
parameter, equal to the thermal velocity in our case (sometimes $b$ is larger
due to small-scale turbulence). The thermal velocity for species X is
$v_{\rm th} = \sqrt{2 k T / m_X}$.

We can evaluate the path integral in many coordinates systems, but for most
problems, it is simplest to use velocity coordinates. At fixed redshift, we
have $dr = c dt = a dx = dv / H$.

In general, the line profile is a Voigt profile $\phi_V = pi^{-1/2} H(a, x)$,
but we restrict to the Doppler profile $\phi_D = pi^{-1/2} \exp[-x^2]$, where
$x$ is the line shift in units of the Doppler broadening $x = (v - v_{lc}) / b$.
Note that the factors of $\pi$ are included so that $\int \phi dx = 1$.

Altogether the optical depth is

\tau_v = \frac{\pi e^2}{m_e c} \frac{\lambda_{lu} f_{lu}}{H(z)}
   \int \frac{n_X}{b} \pi^{-1/2}
   \exp[-\frac{(v - v' - v_\parallel)^2}{b^2}] dv'

The terms outside the integral are a constant for all optical depths $(\frac{\pi
e^2}{m_e c})$, the constant for the transition $(\lambda_{lu} f_{lu})$, and a
constant depending on cosmology $H(z)^{-1}$.

## Usage

The library provides three methods for evaluating the integral above. In all
cases, the user provides regularly spaced arrays of n, v, and T, and the
functions compute the optical depth with a regular spacing on the same domain.
All functions assume that the path data is cell-centered, $n_i = n(x = \Delta x
[i + 1/2])$.

The midpoint function is the simplest scheme, evaluating the integrand for each
input element and summing. This is fastest, but will underestimate the
contribution from thin lines for instance.

The midpoint linear subsample scheme is an improvement. This function assumes
that we can linearly interpolate between the element points and evaluate the
midpoint approximation to the integral with finer spacing. However, linear
interpolation might not be appropriate.

The piecewise constant scheme assumes the element data is constant between
points. In this case, the integral is analytic and the result does not depend
on sampling. This is the most robust method, although it has the disadvantage
of assuming PWC data.

## Tests

We include a test suite in the `test` directory. We use the `Catch` library for
testing. Each `cc` file contains test problems marked with the `TEST_CASE`
macro. The test suite includes a few unit tests and two toy problems.

### Uniform medium test

Integrating over a path where n, v, and T are constant. In this case,

\tau_v = A n \int b^{-1} \pi^{-1/2} \exp[-(v - v' - v_\parallel)^2 / b^2] dv'
       = A n / 2 [ erf( (v - v_a) / b ) - erf( (v - v_b) / b ) ]

where v_a and v_b are the path bounds.

### Gaussian bump test

Integrating over a path where v and T are constant and
n(v) = n_0 exp[-(v - v_0)^2 / a^2]. The solution is analytic, but ugly.

\tau_v = A n_0 \int b^{-1} \pi^{-1/2} \exp[-(v - v_0)^2 / a^2]
                    \exp[-(v - v' - v_\parallel)^2 / b^2] dv'
       = A n_0 a / (2 c) exp[-(v - v_0 - v_\parallel)^2 / c^2]
         [ erf( (a^2 (v_b + v_\parallel - v) + b^2 (v_b - v_0)) / (a b c) )
           - erf( (a^2 (v_a + v_\parallel - v) + b^2 (v_a - v_0)) / (a b c) ) ]

where $c = \sqrt{a^2 + b^2}$. Note that the code always assumes periodic
boundaries, but this test problem is isolated. We get around this by providing
separation between $v_0$ and the domain bounds and only testing near $v_0$.
