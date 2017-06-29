"""
Galaxy shape parameters for use with the lightcones and the HI profile code
by Stewart+ (2014).
All units, unless specified otherwise, are in classical extragalactic HI
units i.e. km/s for velocities and Jy/(km/s) for fluxes
"""

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
    def v_turb(self):
        return None

    @property
    def solid_rot(self):
        return None


# Stewart+ (2014) priors
########################
class Stewart2014(Galaxy):
    """
    More realistic galaxy shape parameters, taken from Stewart et al. (2014)
    (REF?)
    """
    def __init__(self, *args, **kwargs):
        super(Stewart2014, self).__init__(*args, **kwargs)

    @property
    def v_rot(self):
        # rotation velocities
        inclination = np.random.uniform(0, 360., self.n_samples)
        log10_v_rot_true = stats.norm.rvs(
            2.2, 0.3, size=self.n_samples)
        v_rot_true = 10 ** log10_v_rot_true
        v_rot = np.abs(v_rot_true * np.sin(np.radians(inclination)))
        # v_rot = np.clip(v_rot, a_min=50., a_max=None)

        return v_rot

    @property
    def skewness(self):
        # skewness, has to be within [-1, 1]
        skewness = stats.norm.rvs(loc=0, scale=0.2, size=self.n_samples)
        skewness = np.clip(skewness, a_min=-1., a_max=1.)

        return skewness

    @property
    def v_turb(self):
        # turbulent broadening in km/s
        log10_v_turb = stats.norm.rvs(
            1.17, 0.14, size=self.n_samples)
        v_turb = 10 ** log10_v_turb
        v_turb = np.clip(v_turb, a_min=0., a_max=None)

        return v_turb

    @property
    def solid_rot(self):
        # fraction of solid body rotation, must be between 0 and 1
        solid_rot = stats.expon.rvs(scale=0.2, size=self.n_samples)
        solid_rot = np.clip(solid_rot, a_min=0., a_max=1.)
        return solid_rot


# Simple priors
########################
class Simple(Galaxy):
    def __init__(self, *args, **kwargs):
        super(Simple, self).__init__(*args, **kwargs)

    @property
    def v_rot(self):
        # rotation velocities
        inclination = np.random.uniform(0, 360., self.n_samples)
        v_rot = stats.norm.rvs(loc=200, scale=50, size=self.n_samples)
        v_rot = np.abs(v_rot_true * np.sin(np.radians(inclination)))
        v_rot = np.clip(v_rot, a_min=50., a_max=None)

        return v_rot

    @property
    def skewness(self):
        # skewness, has to be within [-1, 1]
        skewness = stats.norm.rvs(loc=0, scale=0.3, size=self.n_samples)
        skewness = np.clip(skewness, a_min=-1., a_max=1.)

        return skewness

    @property
    def v_turb(self):
        # turbulent broadening in km/s
        v_turb = stats.norm.rvs(loc=15., scale=5., size=self.n_samples)
        v_turb = np.clip(v_turb, a_min=0., a_max=None)

        return v_turb

    @property
    def solid_rot(self):
        # fraction of solid body rotation, must be between 0 and 1
        solid_rot = stats.norm.rvs(loc=0.2, scale=0.05, size=self.n_samples)
        solid_rot = np.clip(solid_rot, a_min=0., a_max=None)

        return solid_rot
