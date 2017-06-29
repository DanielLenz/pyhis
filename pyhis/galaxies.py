import numpy as np
from scipy import stats


class Galaxy(object):
    def __init__(self, n_samples=10):
        self.n_samples = n_samples

    @property
    def v_rot(self):
        return None

    @property
    def skewness(self):
        return None

    @property
    def flux(self):
        return None

    @property
    def turb_velo(self):
        return None

    @property
    def solid_rot(self):
        return None


class Stewart2014(Galaxy):
    """
    More realistic galaxy shape parameters, taken from Stewart et al. (2014)
    (REF?)
    """
    def __init__(self, *args, **kwargs):
        super(Simple, self).__init__(*args, **kwargs)

    @property
    def v_rot(self):
        return None

    @property
    def skewness(self):
        return None

    @property
    def flux(self):
        return None

    @property
    def turb_velo(self):
        return None

    @property
    def solid_rot(self):
        return None


class Simple(Galaxy):
    def __init__(self, *args, **kwargs):
        super(Simple, self).__init__(*args, **kwargs)

    @property
    def v_rot(self):
        inclination = np.random.uniform(0, 360., self.n_samples)
        # rotation velocities
        inclination = np.random.uniform(0, 360., self.n_samples)
        rot_velo = stats.norm.rvs(loc=200, scale=50, size=self.n_samples)
        v_rot = np.abs(rot_velo * np.sin(np.radians(inclination)))
        v_rot = np.clip(v_rot, a_min=50., a_max=None)

        return v_rot

    @property
    def skewness(self):
        # skewness, has to be within [-1, 1]
        skewness = stats.norm.rvs(loc=0, scale=0.3, size=self.n_samples)
        skewness = np.clip(skewness, a_min=-1., a_max=1.)

        return skewness

    @property
    def flux(self):
        # flux densities in Jy.km/s
        flux = stats.norm.rvs(loc=50, scale=10, size=self.n_samples)
        flux = np.clip(flux, a_min=0., a_max=None)

        return flux

    @property
    def turb_velo(self):
        # turbulent broadening in km/s
        turb_velo = stats.norm.rvs(loc=15., scale=5., size=self.n_samples)
        turb_velo = np.clip(turb_velo, a_min=0., a_max=None)

        return turb_velo

    @property
    def solid_rot(self):
        # fraction of solid body rotation, must be between 0 and 1
        solid_rot = stats.norm.rvs(loc=0.2, scale=0.05, size=self.n_samples)
        solid_rot = np.clip(solid_rot, a_min=0., a_max=None)

        return solid_rot
