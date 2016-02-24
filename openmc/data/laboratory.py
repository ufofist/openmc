from collections import Iterable
from numbers import Real, Integral

import numpy as np

import openmc.checkvalue as cv
from openmc.stats import Tabular, Univariate, Discrete, Mixture
from .angle_energy import AngleEnergy
from .function import INTERPOLATION_SCHEME
from .endf import get_tab2_record, get_tab1_record


class LaboratoryAngleEnergy(AngleEnergy):
    def __init__(self, breakpoints, interpolation, energy, mu, energy_out):
        super(LaboratoryAngleEnergy).__init__()
        self.breakpoints = breakpoints
        self.interpolation = interpolation
        self.energy = energy
        self.mu = mu
        self.energy_out = energy_out

    @property
    def breakpoints(self):
        return self._breakpoints

    @property
    def interpolation(self):
        return self._interpolation
    @property
    def energy(self):
        return self._energy

    @property
    def mu(self):
        return self._mu

    @property
    def energy_out(self):
        return self._energy_out

    @breakpoints.setter
    def breakpoints(self, breakpoints):
        cv.check_type('laboratory angle-energy breakpoints', breakpoints,
                      Iterable, Integral)
        self._breakpoints = breakpoints

    @interpolation.setter
    def interpolation(self, interpolation):
        cv.check_type('laboratory angle-energy interpolation', interpolation,
                      Iterable, Integral)
        self._interpolation = interpolation

    @energy.setter
    def energy(self, energy):
        cv.check_type('laboratory angle-energy incoming energy', energy,
                      Iterable, Real)
        self._energy = energy

    @mu.setter
    def mu(self, mu):
        cv.check_type('laboratory angle-energy outgoing cosine', mu,
                      Iterable, Univariate)
        self._mu = mu

    @energy_out.setter
    def energy_out(self, energy_out):
        cv.check_iterable_type('laboratory angle-energy outgoing energy',
                               energy_out, Univariate, 2, 2)
        self._energy_out = energy_out

    @classmethod
    def from_endf(cls, file_obj):
        # Laboratory energy-angle distribution
        tab2 = get_tab2_record(file_obj)
        ne = tab2.params[5]
        energy = np.zeros(ne)
        mu = []
        energy_out = []
        for i in range(ne):
            tab2mu = get_tab2_record(file_obj)
            energy[i] = tab2mu.params[1]
            n_mu = tab2mu.params[5]
            mu_i = np.zeros(n_mu)
            p_mu_i = np.zeros(n_mu)
            energy_out_i = []
            for j in range(n_mu):
                params, f = get_tab1_record(file_obj)
                mu_i[j] = params[1]
                p_mu_i[j] = sum(f.y)
                energy_out_i.append(Tabular(f.x, f.y))
            mu.append(Tabular(mu_i, p_mu_i))
            energy_out.append(energy_out_i)

        return cls(tab2.NBT, tab2.INT, energy, mu, energy_out)
