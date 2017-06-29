import numpy as np
import tables

import cygrid as cg
from lineprofile.model import LineModel


def generate_cube(spectra, header, lightcone, outname=None):

    # set up the gridder and kernel
    gridder = cg.WcsGrid(header)

    # set kernel
    kernelsize_fwhm = header['BMAJ']
    kernelsize_sigma = kernelsize_fwhm / 2.355
    sphere_radius = 3. * kernelsize_sigma

    gridder.set_kernel(
        'gauss1d',
        (kernelsize_sigma,),
        sphere_radius,
        kernelsize_sigma / 2.)

    # grid and retrieve cube
    gridder.grid(
        lightcone['ra'].values,
        lightcone['dec'].values,
        spectra)
    cube = gridder.get_datacube()

    # write to disk
    if outname is not None:
        fits.writeto(outname, cube, header)

    return cube


def generate_spectra(lightcone, survey, galaxy, outname=None):
    """
    Generate HI emission spectra, based on the survey type, the light cone,
    and the HI galaxy shape parameters.

    Input
    -----
    lightcone : DataFrame, recArray or similar data structure. Must have
        keys ['ra', 'dec', 'z', 'mhi']
        Lightcone onto which the HI galaxies are placed.
    survey : Instance of a survey from pyhis.surveys
        Contains all relevant survey properties and the target header
    galaxy : Instance of a galaxy from pyhis.surveys
        Contains all galaxy shape parameters
    outname : str, optional
        Filename to save the spectra

    Returns
    -------
    specs : ndarray
        Array of HI spectra

    """
    n_halos = lightcone.shape[0]

    linemodel = LineModel(survey.v_grid, n_disks=1, n_baseline=0)

    specs = []
    for i, halo in enumerate(lightcone.iterrows()):
        spectrum = Spectrum(
            lon=halo[1]['ra'],
            lat=halo[1]['dec'],
            velo_center=halo[1]['velo'],
            flux=halo[1]['flux'],
            v_rot=galaxy.v_rot[i],
            turb_velo=galaxy.turb_velo[i],
            solid_rot=galaxy.solid_rot[i],
            skewness=galaxy.skewness[i],
            linemodel=linemodel)
        specs.append(spectrum.spectrum)

    specs = np.array(specs)

    return specs


class Spectrum(object):

    _model = None
    _spectrum = None
    _shape_parameters = None

    def __init__(
            self,
            lon, lat,
            velo_center, flux, v_rot, turb_velo, solid_rot, skewness,
            linemodel):

        self.velo_center = velo_center
        self.flux = flux
        self.logflux = np.log10(self.flux)
        self.v_rot = v_rot
        self.turb_velo = turb_velo
        self.solid_rot = solid_rot
        self.skewness = skewness
        self.model = linemodel

    @property
    def shape_parameters(self):
        """
        log10 of the flux in Jy.km/s
        apparent rotational velocity in km/s
        turbulent broadening in km/s
        fraction of the solid body rotation, between 0 and 1
        skewness, between -1 and +1
        """
        if self._shape_parameters is None:
            self._shape_parameters = np.array([
                self.logflux,
                self.velo_center,
                self.v_rot,
                self.turb_velo,
                self.solid_rot,
                self.skewness
            ])
        return self._shape_parameters

    @property
    def spectrum(self):
        if self._spectrum is None:
            self._spectrum = self.model.model(self.shape_parameters)
        return self._spectrum
