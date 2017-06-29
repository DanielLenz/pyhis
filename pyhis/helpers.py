import numpy as np
from astropy import units as u
from scipy.constants import c
import pandas as pd

from pyhiframe.pyhiframe import HIConverter


def flux2HImass(S, D):
    """
    Converts observed fluxes into HI masses, depending on the distance

    Input
    -----
    S : float or ndarray
        Flux in Jy km/s
    D : float
        Distance in Mpc
    Returns
    -------
    m_hi : HI mass in solar masses

    """
    if hasattr(D, '__iter__'):
        D = np.array(D)
    if hasattr(S, '__iter__'):
        S = np.array(S)

    m_hi = 2.36e5 * D**2 * S

    return m_hi


def HImass2flux(m_hi, D):
    """
    Converts HI masses into observed fluxes, depending on the distance

    Input
    -----
    M : float or ndarray
        HI mass in solar masses
    D : float
        Distance in Mpc
    Returns
    -------
    S : Flux in Jy km/s

    """
    if hasattr(m_hi, '__iter__'):
        m_hi = np.array(m_hi)
    if hasattr(divmod, '__iter__'):
        D = np.array(D)


    # S = m_hi / 2.36e5 / D**2
    S = m_hi / flux2HImass(1., D)
    return S


def mhalo2mhi(M, z=0):
    """
    From https://arxiv.org/pdf/1608.00007.pdf, Eq. (2) and (14).
    Calculate the HI mass for a given halo mass and redshift.

    Inputs
    ------
    M : float or ndarray
        Mass of the dark matter halo in solar masses
    z : float
        Redshift of the dark matter halo

    Returns
    -------
    M_HI : float or ndarray, same shape as M
        HI mass for a given halo in solar masses
    """
    if hasattr(z, '__iter__'):
        z = np.array(z)
    if hasattr(M, '__iter__'):
        M = np.array(M)

    M10 = 4.58e11
    N10 = 9.89e-3
    b10 = 0.90
    y10 = 0.74

    log10M11 = 1.56
    N11 = 0.009
    b11 = -1.08
    y11 = 4.07
    cHI = 133.66

    log10M1 = np.log10(M10) + z / (z + 1.) * log10M11
    M1 = 10**log10M1
    N1 = N10 + z / (z + 1.) * N11
    b1 = b10 + z / (z + 1.) * b11
    y1 = y10 + z / (z + 1.) * y11

    M_HI = 2 * N1 * M / ((M / M1)**(-b1) + (M / M1)**y1)
    return M_HI


def extend_lightcone(lightcone):
    # HI mass
    mhi = pd.Series(
        mhalo2mhi(lightcone['mhalo'], z=lightcone['z']),
        name='mhi')
    lightcone = lightcone.join(mhi)

    # radial velocity
    converter = HIConverter()
    velos = pd.Series(
        converter.z2velo(lightcone['z']),
        name='velo')
    lightcone = lightcone.join(velos)

    # luminosity distances
    distances = pd.Series(
        converter.z2d(lightcone['z']),
        name='distance')
    lightcone = lightcone.join(distances)

    # HI fluxes
    fluxes = pd.Series(
        HImass2flux(lightcone['mhi'], distances),
        name='flux')
    lightcone = lightcone.join(fluxes)

    return lightcone
