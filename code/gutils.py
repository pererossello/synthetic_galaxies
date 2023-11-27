import sympy as sp
import numpy as np
import astropy.units as u
import gala.potential as gp
from gala.units import galactic


# Define here your parameters for the potential. Check 
# https://gala.adrian.pw/en/latest/potential/index.html
# for info on the available potentials and parameters
# or use
# gp.__dict__

potentials = {
    'LogarithmicPotential': {'vals': [0.4, 1.2, 1, 1, 0.5, np.pi]},
    'LongMuraliBarPotential': {'vals': [1e12, 0.5, 1, 1, np.pi]},
    'NFWPotential': {'vals': [15e10, 1, 1, 1, 1]},
    'MiyamotoNagaiPotential': {'vals': [1e11, 1, 0.75]},
    'SatohPotential': {'vals': [2.5e11, 1, 0.75 ]}
}

for potential in potentials:
    potentials[potential]['pot'] = gp.__dict__[potential](*potentials[potential]['vals'], units=galactic)


def poisson(pot, symbols, simplify=True):
    G = sp.Symbol('G')
    lap = sum(sp.diff(sp.diff(pot, sym), sym) for sym in symbols) /(4*sp.pi*G) 
    return lap.simplify() if simplify else lap

def get_dens_func(potential):
    pot_sym = potential.to_sympy()
    cartesian = list(pot_sym[1].values())
    params = list(pot_sym[2].values())
    print(pot_sym[2])
    tot_params = cartesian + params
    dens = poisson(pot_sym[0], tot_params, simplify=False)
    dens_distr = sp.lambdify(tot_params, dens, 'numpy', )
    print(pot_sym[0])

    return dens_distr, pot_sym[2]

def get_dens_arr(dens_distr, params, lims=(1,1,1), S=50,
                 norm=True):
    x_ = np.linspace(-lims[0], lims[0], S)
    y_ = np.linspace(-lims[1], lims[1], S)
    z_ = np.linspace(-lims[2], lims[2], S)
    X, Y, Z = np.meshgrid(x_, y_, z_, indexing='ij')

    dens_arr = dens_distr(X, Y, Z, *params.values())
    # if for some reason there are negative densities, set to 0
    dens_arr[dens_arr < 0] = 0
    # if nans set to 0
    dens_arr[np.isnan(dens_arr)] = 0

    if norm:
        # normalize the density distribution
        dens_arr /= np.nansum(dens_arr)

    return dens_arr, X, Y, Z

def get_sample(X, Y, Z, dens_arr, N=10000, noise_amplitude=0.1):

    X_flat = X.ravel()
    Y_flat = Y.ravel()
    Z_flat = Z.ravel()

    dens_flat = dens_arr.ravel()
    sampled_indices = np.random.choice(len(dens_flat), size=N, p=dens_flat)

    sampled_X = X_flat[sampled_indices].value
    sampled_Y = Y_flat[sampled_indices].value
    sampled_Z = Z_flat[sampled_indices].value

    sampled_X += np.random.normal(0, noise_amplitude, N)
    sampled_Y += np.random.normal(0, noise_amplitude, N)
    sampled_Z += np.random.normal(0, noise_amplitude, N)

    sampled_points = np.column_stack((sampled_X, sampled_Y, sampled_Z))
    sampled_points = sampled_points*u.kpc
    return sampled_points.T


