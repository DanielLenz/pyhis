import numpy as np
from astropy.io import fits

from pyhiframe.pyhiframe import HIConverter


class Survey(object):

    def __init__(self):

        self.z_grid = np.linspace(
            self.z_min,
            self.z_max,
            self.n_channels)

        # in frequency
        self.nu_max, self.nu_min = self.converter.z2nu([
            self.z_min, self.z_max])
        self.nu_grid = self.converter.z2nu(self.z_grid)

        # in v_radial, units of km/s
        self.v_min, self.v_max = self.converter.z2velo(
            [self.z_min, self.z_max])
        self.v_grid = self.converter.z2velo(self.z_grid)

    @property
    def header(self):
        return None

    def summary(self):
        print(vars(self))

    @property
    def converter(self):
        return HIConverter(mode='relativistic')


class GbtType(Survey):
    def __init__(self, *args, **kwargs):

        # survey properties
        # spatial
        self.boxsize = np.sqrt(2.)  # 10 * 10 deg field
        self.fwhm = 3. / 60.  # fwhm in degree
        self.n_pixels = 256

        # spectral
        self.n_channels = 512

        # redshift
        self.z_min, self.z_max = 0.5, 1.

        super(GbtType, self).__init__(*args, **kwargs)

    @property
    def header(self):

        # create header
        header = fits.Header()

        header['NAXIS'] = 3

        header['BMAJ'] = self.fwhm
        header['BMIN'] = self.fwhm

        header['NAXIS1'] = self.n_pixels
        header['NAXIS2'] = self.n_pixels
        header['NAXIS3'] = self.n_channels

        header['CDELT1'] = self.boxsize / self.n_pixels
        header['CDELT2'] = self.boxsize / self.n_pixels
        header['CDELT3'] = self.nu_grid[1] - self.nu_grid[0]

        header['CRPIX1'] = 0
        header['CRPIX2'] = 0
        header['CRPIX3'] = 0

        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CUNIT3'] = 'MHz'

        header['CRVAL1'] = 0.
        header['CRVAL2'] = 0.
        header['CRVAL3'] = self.nu_grid[0]
        header['LATPOLE'] = 90.

        header['CTYPE1'] = 'RA---SFL'
        header['CTYPE2'] = 'DEC--SFL'
        header['CTYPE3'] = 'FREQ'

        return header
