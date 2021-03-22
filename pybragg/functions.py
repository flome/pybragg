#########################################################
#
# modified based on this code: https://de.mathworks.com/matlabcentral/fileexchange/63405-bragg-peak-analysis
# original author: Jan Gajewski (https://de.mathworks.com/matlabcentral/profile/authors/4307259)
# email: jan.gajewski@ifj.edu.pl
#
# author of this file: Florian Mentzel
# email: florian.mentzel@tu-dortmund.de
# written: 20.03.2021
#
#########################################################

# needed for spline
from scipy.interpolate import interp1d
# parabolic cylinder function D
from scipy.special import pbdv
from scipy.optimize import curve_fit
import numpy as np

def fitBP(z, D, method='bortfeld', rel_resolution=0.01):
    # check for validity of relevant input arguments
    assert len(z) == len(D), f"z and D need to have same length but are len(z)={len(z)} and len(D)={len(D)}"
    assert method in ['spline', 'bortfeld'], f"method can only be 'spline' or 'bortfeld' but is {method}"

    # define some accuracy settings
    resolution = rel_resolution*np.min(np.diff(z))

    # fit spline with given precision to curve
    spline_func = interp1d(z, D, kind='cubic')
    z_spline    = np.linspace(min(z), max(z), round((max(z)-min(z)) / resolution ))
    quantities  = characterize_z_D_curve(z_spline, spline_func(z_spline))

    # return interpolation for reference as well
    quantities['spline_func'] = spline_func

    # if spline fit only, return already
    if method == 'spline':
        return quantities

    # use precomputed values for further computations using bortfeld fit
    if method == 'bortfeld':
        # create init parameters from
        p         = 1.77 # result from paper Fig. 1
        alpha     = 0.0022 # result from paper Fig. 1
        R0        = quantities['R80D']
        E0        = (R0/alpha)**(1/p) # paper: Eq. (4)
        sigmaMono = (0.012*R0**0.935)/10 # paper: Eq. (18)
        sigmaE0   = 0.01*E0 # assumtion that Delta E will be small
        sigma     = np.sqrt(sigmaMono**2+sigmaE0**2*alpha**2*p**2*E0**(2*p-2)) # paper: Eq. (19)
        eps       = 0.1 # assumption
        phi       = 0.1*quantities['D100'] # assumption

        # create two definition ranges according to paper Eq. (27)
        first_window  = z < R0 - 10*sigma
        second_window = (z >= R0 - 10*sigma) & (z <= R0 + 5*sigma)
        # fit only relevant part, rest will be zero anyway
        z_fit = z [ first_window | second_window ]
        D_fit = D [ first_window | second_window ]
        p, c = curve_fit(
            bortfeld, z_fit, D_fit,
            p0 = [R0, sigma, phi, eps],
            bounds=( # limits
                #     R0       |   sigma  |  phi  | eps
                (R0 - 1*sigma, 0.5*sigma,       0,  -10),
                (R0 + 1*sigma,   3*sigma,  np.inf,   10)
            ),

        )

        # return for easy access
        quantities['bortfeld_fit_results'] = {var: {'nominal': nom, 'std': std} for var, nom, std in zip(['R0', 'sigma', 'phi', 'epsilon'], p, np.diag(c))}

        # return parameter vector and cov matrix as well
        quantities['bortfeld_fit_p']       = p
        quantities['bortfeld_fit_cov']     = c

        # recalc some quantities if needed
        bortfeld_quantities = characterize_z_D_curve(z_spline, bortfeld(z_spline, *p))
        # overwrite results from spline fit
        quantities.update(bortfeld_quantities)

        return quantities

def bortfeld(z, R0, sigma, phi, eps):
    # create two definition ranges according to paper Eq. (27)
    first_window  = z < R0 - 10*sigma
    second_window = (z >= R0 - 10*sigma) & (z <= R0 + 5*sigma)

    D_Dhat = D_hat(z[first_window], R0, phi, eps)
    D_D    = D(z[second_window], R0, sigma, phi, eps)

    values = np.zeros_like(z)
    values[first_window]  = D_Dhat
    values[second_window] = D_D
    return values

def D_hat(z, R0, phi, eps): # paper: Eq. (28)
    first        = phi / ( 1 + 0.012 * R0)
    brack_first  = 17.93 * ( R0 - z )**(-0.435)
    brack_second = 0.444 + 31.7 * (eps/R0) * ( R0 - z )**(0.565)
    return first * (brack_first + brack_second)

def D(z, R0, sigma, phi, eps): # paper: Eq. (29)
    first        = phi * np.exp( - ( (R0 - z)/(2*sigma) )**2 ) * sigma**(0.565)
    second       = 1 + 0.012 * R0
    brack_first  = 11.26 / sigma * pbdv(-0.565, ( -(R0 - z)/sigma ) )[0]
    brack_second = (0.157 + 11.26 * eps/R0) * pbdv(-1.565, ( -(R0 - z)/sigma ) )[0]
    return first/second * ( brack_first + brack_second)

def characterize_z_D_curve(z, D):
    results = {}

    # compute quantities
    D100_index = np.argmax(D)
    D100       = D[D100_index]
    R100       = z[D100_index]
    results.update({
        'D100': D100,
        'R100': R100
    })

    # split at peak index into proximal and distal part
    z_proximal    = z[:D100_index]
    dose_proximal = D[:D100_index]
    z_distal      = z[D100_index:]
    dose_distal   = D[D100_index:]

    R90P = z_proximal[np.argmin(np.abs(dose_proximal - 0.9 * D100))]
    R90D = z_distal  [np.argmin(np.abs(dose_distal   - 0.9 * D100))]
    R80P = z_proximal[np.argmin(np.abs(dose_proximal - 0.8 * D100))]
    R80D = z_distal  [np.argmin(np.abs(dose_distal   - 0.8 * D100))]
    R50P = z_proximal[np.argmin(np.abs(dose_proximal - 0.5 * D100))]
    R50D = z_distal  [np.argmin(np.abs(dose_distal   - 0.5 * D100))]
    R20D = z_distal  [np.argmin(np.abs(dose_distal   - 0.2 * D100))]
    R10D = z_distal  [np.argmin(np.abs(dose_distal   - 0.1 * D100))]
    results.update({
        'R90P': R90P,
        'R90D': R90D,
        'R80P': R80P,
        'R80D': R80D,
        'R50P': R50P,
        'R50D': R50D,
        'R20D': R20D,
        'R10D': R10D
    })

    FWHM    = R50D  - R50P
    DFO2080 = R20D  - R80D
    DFO1090 = R10D  - R90D
    results.update({
        'FWHM': FWHM,
        'DFO2080': DFO2080,
        'DFO1090': DFO1090,
    })

    return results
