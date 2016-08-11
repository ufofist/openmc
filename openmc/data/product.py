from collections import Iterable
from io import StringIO
from numbers import Real
import sys

import numpy as np

import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin
from openmc.stats import Tabular, Legendre
from .angle_energy import AngleEnergy
from .angle_distribution import AngleDistribution
from .correlated import CorrelatedAngleEnergy
from .data import ATOMIC_SYMBOL
from .endf import get_head_record, get_tab1_record, get_tab2_record, \
    get_list_record, get_cont_record
from .function import Tabulated1D, Polynomial, Function1D
from .kalbach_mann import KalbachMann
from .laboratory import LaboratoryAngleEnergy
from .nbody import NBodyPhaseSpace
from .uncorrelated import UncorrelatedAngleEnergy

if sys.version_info[0] >= 3:
    basestring = str


def get_products(ev, mt):
    file_obj = StringIO(ev.section[6, mt])

    # Read HEAD record
    items = get_head_record(file_obj)
    reference_frame = {1: 'laboratory', 2: 'center-of-mass',
                       3: 'light-heavy'}[items[3]]
    n_products = items[4]

    products = []
    for i in range(n_products):
        # Get yield for this product
        params, yield_ = get_tab1_record(file_obj)

        za = params[0]
        awr = params[1]
        lip = params[2]
        law = params[3]

        if za == 0:
            p = Product('photon')
        elif za == 1:
            p = Product('neutron')
        elif za == 1000:
            p = Product('electron')
        else:
            z = za // 1000
            a = za % 1000
            p = Product('{}{}'.format(ATOMIC_SYMBOL[z], a))

        p.yield_ = yield_

        """
        # Set reference frame
        if reference_frame == 'laboratory':
            p.center_of_mass = False
        elif reference_frame == 'center-of-mass':
            p.center_of_mass = True
        elif reference_frame == 'light-heavy':
            p.center_of_mass = (awr <= 4.0)
        """

        if law == 0:
            # No distribution given
            pass
        if law == 1:
            # Continuum energy-angle distribution
            tab2 = get_tab2_record(file_obj)
            lang = tab2.params[2]
            if lang == 1:
                p.distribution = [CorrelatedAngleEnergy.from_endf(
                    file_obj, tab2)]
            elif lang == 2:
                p.distribution = [KalbachMann.from_endf(file_obj, tab2)]

        elif law == 2:
            # Discrete two-body scattering
            tab2 = get_tab2_record(file_obj)
            ne = tab2.params[5]
            energy = np.zeros(ne)
            mu = []
            for i in range(ne):
                items, values = get_list_record(file_obj)
                energy[i] = items[1]
                lang = items[2]
                if lang == 0:
                    mu.append(Legendre(values))
                elif lang == 12:
                    mu.append(Tabular(values[::2], values[1::2]))
                elif lang == 14:
                    mu.append(Tabular(values[::2], values[1::2],
                                      'log-linear'))

            angle_dist = AngleDistribution(energy, mu)
            dist = UncorrelatedAngleEnergy(angle_dist)
            p.distribution = [dist]
            # TODO: Add level-inelastic info?

        elif law == 3:
            # Isotropic discrete emission
            p.distribution = [UncorrelatedAngleEnergy()]
            # TODO: Add level-inelastic info?

        elif law == 4:
            # Discrete two-body recoil
            pass

        elif law == 5:
            # Charged particle elastic scattering
            pass

        elif law == 6:
            # N-body phase-space distribution
            p.distribution = [NBodyPhaseSpace.from_endf(file_obj)]

        elif law == 7:
            # Laboratory energy-angle distribution
            p.distribution = [LaboratoryAngleEnergy.from_endf(file_obj)]

        products.append(p)

    return products


class Product(EqualityMixin):
    """Secondary particle emitted in a nuclear reaction

    Parameters
    ----------
    particle : str, optional
        What particle the reaction product is. Defaults to 'neutron'.

    Attributes
    ----------
    applicability : Iterable of openmc.data.Tabulated1D
        Probability of sampling a given distribution for this product.
    decay_rate : float
        Decay rate in inverse seconds
    distribution : Iterable of openmc.data.AngleEnergy
        Distributions of energy and angle of product.
    emission_mode : {'prompt', 'delayed', 'total'}
        Indicate whether the particle is emitted immediately or whether it
        results from the decay of reaction product (e.g., neutron emitted from a
        delayed neutron precursor). A special value of 'total' is used when the
        yield represents particles from prompt and delayed sources.
    particle : str
        What particle the reaction product is.
    yield_ : openmc.data.Function1D
        Yield of secondary particle in the reaction.

    """

    def __init__(self, particle='neutron'):
        self.particle = particle
        self.decay_rate = 0.0
        self.emission_mode = 'prompt'
        self.distribution = []
        self.applicability = []
        self.yield_ = Polynomial((1,))  # 0-order polynomial i.e. a constant

    def __repr__(self):
        if isinstance(self.yield_, Real):
            return "<Product: {}, emission={}, yield={}>".format(
                self.particle, self.emission_mode, self.yield_)
        elif isinstance(self.yield_, Tabulated1D):
            if np.all(self.yield_.y == self.yield_.y[0]):
                return "<Product: {}, emission={}, yield={}>".format(
                    self.particle, self.emission_mode, self.yield_.y[0])
            else:
                return "<Product: {}, emission={}, yield=tabulated>".format(
                    self.particle, self.emission_mode)
        else:
            return "<Product: {}, emission={}, yield=polynomial>".format(
                self.particle, self.emission_mode)

    @property
    def applicability(self):
        return self._applicability

    @property
    def decay_rate(self):
        return self._decay_rate

    @property
    def distribution(self):
        return self._distribution

    @property
    def emission_mode(self):
        return self._emission_mode

    @property
    def particle(self):
        return self._particle

    @property
    def yield_(self):
        return self._yield

    @applicability.setter
    def applicability(self, applicability):
        cv.check_type('product distribution applicability', applicability,
                      Iterable, Tabulated1D)
        self._applicability = applicability

    @decay_rate.setter
    def decay_rate(self, decay_rate):
        cv.check_type('product decay rate', decay_rate, Real)
        cv.check_greater_than('product decay rate', decay_rate, 0.0, True)
        self._decay_rate = decay_rate

    @distribution.setter
    def distribution(self, distribution):
        cv.check_type('product angle-energy distribution', distribution,
                      Iterable, AngleEnergy)
        self._distribution = distribution

    @emission_mode.setter
    def emission_mode(self, emission_mode):
        cv.check_value('product emission mode', emission_mode,
                       ('prompt', 'delayed', 'total'))
        self._emission_mode = emission_mode

    @particle.setter
    def particle(self, particle):
        cv.check_type('product particle type', particle, basestring)
        self._particle = particle

    @yield_.setter
    def yield_(self, yield_):
        cv.check_type('product yield', yield_, Function1D)
        self._yield = yield_

    def to_hdf5(self, group):
        """Write product to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """
        group.attrs['particle'] = np.string_(self.particle)
        group.attrs['emission_mode'] = np.string_(self.emission_mode)
        if self.decay_rate > 0.0:
            group.attrs['decay_rate'] = self.decay_rate

        # Write yield
        self.yield_.to_hdf5(group, 'yield')

        # Write applicability/distribution
        group.attrs['n_distribution'] = len(self.distribution)
        for i, d in enumerate(self.distribution):
            dgroup = group.create_group('distribution_{}'.format(i))
            if self.applicability:
                self.applicability[i].to_hdf5(dgroup, 'applicability')
            d.to_hdf5(dgroup)

    @classmethod
    def from_hdf5(cls, group):
        """Generate reaction product from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.Product
            Reaction product

        """
        particle = group.attrs['particle'].decode()
        p = cls(particle)

        p.emission_mode = group.attrs['emission_mode'].decode()
        if 'decay_rate' in group.attrs:
            p.decay_rate = group.attrs['decay_rate']

        # Read yield
        p.yield_ = Function1D.from_hdf5(group['yield'])

        # Read applicability/distribution
        n_distribution = group.attrs['n_distribution']
        distribution = []
        applicability = []
        for i in range(n_distribution):
            dgroup = group['distribution_{}'.format(i)]
            if 'applicability' in dgroup:
                applicability.append(Tabulated1D.from_hdf5(
                    dgroup['applicability']))
            distribution.append(AngleEnergy.from_hdf5(dgroup))

        p.distribution = distribution
        p.applicability = applicability

        return p
