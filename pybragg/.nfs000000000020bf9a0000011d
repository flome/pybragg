#########################################################
#
# author of this file: Florian Mentzel
# email: florian.mentzel@tu-dortmund.de
# written: 20.03.2021
#
#########################################################

# needed for spline
from scipy.interpolate import interp1d
# parabolic cylinder function D
from scipy import special
from scipy.optimize import curve_fit
import numpy as np


def fitBP(z, D, method='bortfeld', rel_resolution=0.01):
    """ Automated fit and characterization of a pragg curve.
    
    Parameters
    -----------
    :param z: depth in phantom in cm
    :param D: dose at depth z
    :param method: "bortfeld" for full fit with Bortfeld approximation. "spline" for fast and simple spline fit (default "bortfeld")
    :param rel_resolution: fraction of z step width for the fit function characterization for range quantities like R80. (default 0.01)
  
    Returns
    --------
    :returns:  D(z) - depth dose in depth z
    """

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
        sigmaMono = (0.012*R0**0.935) # paper: Eq. (18)
        sigmaE0   = 0.01*E0 # assumtion that Delta E will be small
        sigma     = np.sqrt(sigmaMono**2+sigmaE0**2*alpha**2*p**2*E0**(2*p-2)) # paper: Eq. (19)

        # normalization constant for fit
        A = quantities['D100']
        # normalization constant for second part of equation, depends on the epsilon from the original publication
        k = 0.01394
        p, c = curve_fit(
            bortfeld, z, D,
            p0 = [A, R0, sigma, p, k],
            bounds=( # limits
            #               D100        |    R0     |  sigma   |   p  |  k  | 
                 ( .5*quantities['D100'], R0-3*sigma, 0.5*sigma,   0.5,    0),
                 (1.5*quantities['D100'], R0+3*sigma, 3  *sigma,   2.5,  0.1),
             )
           )
        
        # return for easy access
        quantities['bortfeld_fit_results'] = {var: {'nominal': nom, 'std': std} for var, nom, std in zip(['D100', 'R0', 'sigma', 'p', 'k'], p, np.diag(c))}

        # return parameter vector and cov matrix as well
        quantities['bortfeld_fit_p']       = p
        quantities['bortfeld_fit_cov']     = c

        # recalc some quantities if needed
        bortfeld_quantities = characterize_z_D_curve(z_spline, bortfeld(z_spline, *p))
        # overwrite results from spline fit
        quantities.update(bortfeld_quantities)

        return quantities


# adapted from https://gray.mgh.harvard.edu/attachments/article/293/BraggCurve.py
def cyl_gauss(a,x):
    """Calculate product of Gaussian and parabolic cylinder function"""
    y = np.zeros_like(x)
    branch = -12.0   #for large negative values of the argument we run into numerical problems, need to approximate result

    x1 = x[x<branch]
    y1 = np.sqrt(2*np.pi)/special.gamma(-a)*(-x1)**(-a-1)
    y[x<branch] = y1

    x2 = x[x>=branch]
    y2a = special.pbdv(a,x2)[0]     #special function yielding parabolic cylinder function, first array [0] is function itself
    y2b = np.exp(-x2*x2/4)
    y2 = y2a*y2b

    y[x>=branch] = y2

    return y


def bortfeld(z, D100, R0, sigma, p=1.77, k=0.01394):
    """ Bortfeld function to approximate a Bragg curve.
    
    Parameters
    -----------
    :param z: depth in phantom in cm
    :param D100: approximate maximum Dose (height of peak)
    :param R0: bragg curve range in cm (distance to 80% D100 on distal side)
    :param sigma: sigma term, measure for peak width
    :param p: exponent of "Bragg-Kleeman rule" (default 1.77)
    :param k: scaling factor w/ dependence on epsilon, (default 0.01394) 
  
    Returns
    --------
    :returns:  D(z) - depth dose in depth z
    """
    
    return 0.65 * D100 * ( cyl_gauss( -1/p, (z-R0)/sigma ) + sigma*k*cyl_gauss( -1/p-1, (z-R0)/sigma ) )


def characterize_z_D_curve(z, D):
    """ Method that computes dose and range quantities from a given bragg curve
    
    Parameters
    -----------
    :param z: depth in phantom in cm
    :param D: dose at depth z
    
    Returns
    --------
    :returns: 
      - results (dict): Ranges to certain fractions of Dmax on distal (D) and proximal (P) side of peak. Also: FWHM, DFO(1090)/(2080)
    """
    
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
