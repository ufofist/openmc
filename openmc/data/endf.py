"""Module for parsing and manipulating data from ENDF evaluations.

All the classes and functions in this module are based on document
ENDF-102 titled "Data Formats and Procedures for the Evaluated Nuclear
Data File ENDF-6". The latest version from June 2009 can be found at
http://www-nds.iaea.org/ndspub/documents/endf/endf102/endf102.pdf

"""
from __future__ import print_function, division, unicode_literals

import io
import re
import os
from math import pi
from collections import OrderedDict, Iterable

import numpy as np
from numpy.polynomial.polynomial import Polynomial

from . import reaction_name
from .container import Tabulated1D, interpolation_scheme
from .energy_distribution import *
from .angle_distribution import AngleDistribution
from .thermal import CoherentElastic
from openmc.stats.univariate import Uniform, Tabular, Legendre


libraries = {0: 'ENDF/B', 1: 'ENDF/A', 2: 'JEFF', 3: 'EFF',
             4: 'ENDF/B High Energy', 5: 'CENDL', 6: 'JENDL',
             31: 'INDL/V', 32: 'INDL/A', 33: 'FENDL', 34: 'IRDF',
             35: 'BROND', 36: 'INGDB-90', 37: 'FENDL/A', 41: 'BROND'}


def radiation_type(value):
    p = {0: 'gamma', 1: 'beta-', 2: 'ec/beta+', 3: 'IT',
         4: 'alpha', 5: 'neutron', 6: 'sf', 7: 'proton',
         8: 'e-', 9: 'xray', 10: 'unknown'}
    if value % 1.0 == 0:
        return p[int(value)]
    else:
        return (p[int(value)], p[int(10*value % 10)])


def endftod(s):
    if 'e' in s or 'E' in s:
        return float(s)
    elif '+' in s[1:]:
        return float(s[0] + s[1:].replace('+', 'e+'))
    elif '-' in s[1:]:
        return float(s[0] + s[1:].replace('-', 'e-'))
    else:
        return float(s)


def at_end_of_tape(f):
    """Indicate whether file is positioned at the end of an ENDF tape.

    Parameters
    ----------
    f : file_like
        File to check

    Returns
    -------
    bool
        Whether the file is at the end of the ENDF tape

    """
    position = f.tell()
    line = f.readline()
    if line == '' or line[66:70] == '  -1':
        return True
    else:
        f.seek(position)
        return False


def seek_material_end(f):
    """Position the file at the end of the ENDF material (MAT) currently being read.

    Parameters
    ----------
    f : file_like
        File to position

    """
    while True:
        line = f.readline()
        if line[66:70] == '   0':
            break


def seek_file_end(f):
    """Position the file at the end of the ENDF file (MF) currently being read.

    Parameters
    ----------
    f : file_like
        File to position

    """
    while True:
        line = f.readline()
        if line[70:72] == ' 0':
            break


def seek_section_end(f):
    """Position the file at the end of the ENDF section (MT) currently being read.

    Parameters
    ----------
    f : file_like
        File to position

    """
    while True:
        line = f.readline()
        if line[72:75] == '  0':
            break


class Evaluation(object):
    """ENDF material evaluation with multiple files/sections

    The Evaluation class provides a means to parse data from an ENDF-6 format
    file and access it as stored internal Python objects. A summary of the
    parsing capabilities is as follows:

    == === =============================================== ========
    MF MT  Description                                     Complete
    == === =============================================== ========
    1  451 Descriptive data and directory                  Yes
    1  452 Number of neutrons per fission                  Yes
    1  455 Delayed neutron data                            Yes
    1  456 Number of prompt neutrons per fission           Yes
    1  458 Components of fission energy release            Yes
    1  460 Delayed photon data                             Yes
    2  151 Resonance parameters                            Yes
    3   -  Reaction cross sections                         Yes
    4   -  Angular distributions                           Yes
    5   -  Energy distributions                            Yes
    6   -  Product energy-angle distributions              Yes
    7  2   Thermal elastic scattering                      Yes
    7  4   Thermal inelastic scattering                    Yes
    8  454 Independent fission yields                      Yes
    8  457 Radioactive decay data                          Yes
    8  459 Cumulative fission yields                       Yes
    8   -  Radioactive nuclide production                  Yes
    9   -  Multiplicities of radioactive products          Yes
    10  -  Production cross sections for radionuclides     Yes
    12  -  Photon production yield data                    Yes
    13  -  Photon production cross sections                Yes
    14  -  Photon angular distributions                    Yes
    15  -  Continuous photon energy spectra                Yes
    23  -  Photon and electron cross sections              Yes
    26  -  Secondary distributions for electro-atomic data Yes
    27  -  Atomic form factors                             Yes
    28 533 Atomic relaxation data                          Yes
    30 1   Directory and correspondance table              No
    30 2   Covariance matrix                               No
    30  -  Sensitivities                                   No
    31  -  Covariances of fission                          No
    32  -  Covariances of resonance parameters             No
    33  -  Covariances of neutron cross sections           No
    34  -  Covariances for angular distributions           No
    35  -  Covariances for energy distributions            No
    40  -  Covariances for radionuclide production         No
    == === =============================================== ========

    Attributes
    ----------
    atomic_relaxation : dict
        Dictionary containing atomic relaxation data from MF=28, MT=533. If the
        evaluation is not an atomic relaxation sublibrary, the dictionary is
        empty.
    decay : dict
        Dictionary containing decay data from MF=8. If the evaluation is not
        from a decay sublibrary, the dictionary is empty.
    fission : dict
        Dictionary containing fission-related data, such as neutrons release
        from fission (MF=1, MT=452,455,456), components of energy release (MF=1,
        MT=458), delayed photons from fission (MF=1, MT=460), and
        cumulative/independent fission yields (MF=8, MT=454,459).
    info : dict
        Miscallaneous information about the evaluation.
    target : dict
        Information about the target material, such as its mass, isomeric state,
        whether it's stable, and whether it's fissionable.
    projectile : dict
        Information about the projectile such as its mass.
    reaction_list : list of 4-tuples
        List of sections in the evaluation. The entries of the tuples are the
        file (MF), section (MT), number of records (NC), and modification
        indicator (MOD).
    reactions : collections.OrderedDict
        Dictionary whose keys are MT numbers and values are Reaction instances.
    resonances : dict
        Resolved resonance data from MF=2, MT=151.
    thermal_elastic : dict
        Coherent and/or incoherent thermal elastic data from MF=7, MT=2.
    thermal_inelastic : dict
        Incoherent thermal inelastic data from MF=7, MT=4.

    """

    def __init__(self, filename_or_handle, verbose=True):
        if isinstance(filename_or_handle, io.IOBase):
            self._fh = filename_or_handle
        else:
            self._fh = open(filename_or_handle, 'rU')
        self._verbose = verbose
        self._veryverbose = False

        # Create public attributes
        self.atomic_relaxation = {}
        self.decay = {}
        self.fission = {'nu': {}, 'energy_release': {}, 'delayed_photon': {},
                        'yield_independent': {}, 'yield_cumulative': {}}
        self.info = {}
        self.target = {}
        self.projectile = {}
        self.reaction_list = []
        self.reactions = OrderedDict()
        self.resonances = {}
        self.thermal_elastic = {}
        self.thermal_inelastic = {}

        # Determine MAT number for this evaluation
        MF = 0
        while MF == 0:
            position = self._fh.tell()
            line = self._fh.readline()
            MF = int(line[70:72])
        self.material = int(line[66:70])

        # Save starting position for this evaluation
        self._fh.seek(position)

        # First we need to read MT=1, MT=451 which has a description of the ENDF
        # file and a list of what data exists in the file
        self._read_header()

        # Save starting position
        self._start_position = self._fh.tell()

    def read(self, reactions=None, skip_mf=[], skip_mt=[]):
        """Reads reactions from the ENDF file of the Evaluation object. If no
        arguments are provided, this method will read all the reactions in the
        file. A single reaction can be read if provided.

        Parameters
        ----------
        reactions : tuple or list of tuple, optional
            A single reaction in the following format: (MF, MT)
        skip_mf : list of int, optional
            Files (MF) which should not be read
        skip_mt : list of int, optional
            Reactions (MT) which should not be read

        """

        # Make sure file is positioned correctly
        self._fh.seek(self._start_position)

        if isinstance(reactions, tuple):
            reactions = [reactions]

        while True:
            # Find next section
            while True:
                position = self._fh.tell()
                line = self._fh.readline()
                MAT = int(line[66:70])
                MF = int(line[70:72])
                MT = int(line[72:75])
                if MT > 0 or MAT == 0:
                    self._fh.seek(position)
                    break

            # If end of material reached, exit loop
            if MAT == 0:
                break

            # If there are files/reactions requested to be skipped, check them
            if MF in skip_mf:
                seek_file_end(self._fh)
                continue
            if MT in skip_mt:
                seek_section_end(self._fh)
                continue

            # If reading is restricted to certain reactions, check here
            if reactions and (MF, MT) not in reactions:
                seek_section_end(self._fh)
                continue


            # File 1 data
            if MF == 1:
                if MT == 452:
                    # Number of total neutrons per fission
                    self._read_total_nu()
                elif MT == 455:
                    # Number of delayed neutrons per fission
                    self._read_delayed_nu()
                elif MT == 456:
                    # Number of prompt neutrons per fission
                    self._read_prompt_nu()
                elif MT == 458:
                    # Components of energy release due to fission
                    self._read_fission_energy()
                elif MT == 460:
                    self._read_delayed_photon()

            elif MF == 2:
                # Resonance parameters
                if MT == 151:
                    self._read_resonances()
                else:
                    seek_section_end(self._fh)

            elif MF == 3:
                # Reaction cross sections
                self._read_reaction_xs(MT)

            elif MF == 4:
                # Angular distributions
                self._read_angular_distribution(MT)

            elif MF == 5:
                # Energy distributions
                self._read_energy_distribution(MT)

            elif MF == 6:
                # Product energy-angle distributions
                self._read_product_energy_angle(MT)

            elif MF == 7:
                # Thermal scattering data
                if MT == 2:
                    self._read_thermal_elastic()
                if MT == 4:
                    self._read_thermal_inelastic()

            elif MF == 8:
                # decay and fission yield data
                if MT == 454:
                    self._read_independent_yield()
                elif MT == 459:
                    self._read_cumulative_yield()
                elif MT == 457:
                    self._read_decay()
                else:
                    self._read_radioactive_nuclide(MT)

            elif MF == 9:
                # multiplicities
                self._read_multiplicity(MT)

            elif MF == 10:
                # cross sections for production of radioactive nuclides
                self._read_production_xs(MT)

            elif MF == 12:
                # Photon production yield data
                self._read_photon_production_yield(MT)

            elif MF == 13:
                # Photon production cross sections
                self._read_photon_production_xs(MT)

            elif MF == 14:
                # Photon angular distributions
                self._read_photon_angular_distribution(MT)

            elif MF == 15:
                # Photon continuum energy distributions
                self._read_photon_energy_distribution(MT)

            elif MF == 23:
                # photon interaction data
                self._read_photon_interaction(MT)

            elif MF == 26:
                # secondary distributions for photon interactions
                self._read_electron_products(MT)

            elif MF == 27:
                # atomic form factors or scattering functions
                self._read_scattering_functions(MT)

            elif MF == 28:
                # atomic relaxation data
                self._read_atomic_relaxation()

            else:
                seek_file_end(self._fh)

    def _read_header(self):
        self._print_info(1, 451)

        # Information about target/projectile
        # First HEAD record
        items = self._get_head_record()
        self.target['ZA'] = items[0]
        self.target['mass'] = items[1]
        self._LRP = items[2]
        self.target['fissionable'] = (items[3] == 1)
        try:
            global libraries
            library = libraries[items[4]]
        except KeyError:
            library = 'Unknown'
        self.info['modification'] = items[5]

        # Control record 1
        items = self._get_cont_record()
        self.target['excitation_energy'] = items[0]
        self.target['stable'] = (int(items[1]) == 0)
        self.target['state'] = items[2]
        self.target['isomeric_state'] = items[3]
        self.info['format'] = items[5]
        assert self.info['format'] == 6

        # Control record 2
        items = self._get_cont_record()
        self.projectile['mass'] = items[0]
        self.info['energy_max'] = items[1]
        library_release = items[2]
        self.info['sublibrary'] = items[4]
        library_version = items[5]
        self.info['library'] = (library, library_version, library_release)

        # Control record 3
        items = self._get_cont_record()
        self.target['temperature'] = items[0]
        self.info['derived'] = (items[2] > 0)
        NWD = items[4]
        NXC = items[5]

        # Text records
        text = [self._get_text_record() for i in range(NWD)]
        if len(text) >= 5:
            self.target['zsymam'] = text[0][0:11]
            self.info['laboratory'] = text[0][11:22]
            self.info['date'] = text[0][22:32]
            self.info['author'] = text[0][32:66]
            self.info['reference'] = text[1][1:22]
            self.info['date_distribution'] = text[1][22:32]
            self.info['date_release'] = text[1][33:43]
            self.info['date_entry'] = text[1][55:63]
            self.info['identifier'] = text[2:5]
            self.info['description'] = text[5:]

        # File numbers, reaction designations, and number of records
        for i in range(NXC):
            items = self._get_cont_record(skipC=True)
            MF, MT, NC, MOD = items[2:6]
            self.reaction_list.append((MF, MT, NC, MOD))

    def _read_total_nu(self):
        self._print_info(1, 452)

        # Determine representation of total nu data
        items = self._get_head_record()
        LNU = items[3]

        # Polynomial representation
        if LNU == 1:
            coefficients = np.asarray(self._get_list_record(onlyList=True))
            self.fission['nu']['total'] = Polynomial(coefficients)

        # Tabulated representation
        elif LNU == 2:
            params, self.fission['nu']['total'] = self._get_tab1_record()

        # Skip SEND record
        self._fh.readline()

    def _read_delayed_nu(self):
        self._print_info(1, 455)

        # Create delayed nu reaction
        self.fission['nu']['delayed'] = {}

        # Determine representation of delayed nu data
        items = self._get_head_record()
        LDG = items[2]
        LNU = items[3]
        self.fission['nu']['delayed']['decay_energy_dependent'] = (LDG == 1)

        if LDG == 0:
            # Delayed-group constants energy independent
            self.fission['nu']['delayed']['decay_constants'] = np.asarray(
                self._get_list_record(onlyList=True))
        elif LDG == 1:
            # Delayed-group constants energy dependent
            raise NotImplementedError

        if LNU == 1:
            # Nu represented as polynomial
            coefficients = np.asarray(self._get_list_record(onlyList=True))
            self.fission['nu']['delayed']['values'] = Polynomial(coefficients)
        elif LNU == 2:
            # Nu represented by tabulation
            params, self.fission['nu']['delayed']['values'] = self._get_tab1_record()
        self._fh.readline()

    def _read_prompt_nu(self):
        self._print_info(1, 456)

        # Determine representation of delayed nu data
        items = self._get_head_record()
        LNU = items[3]

        if LNU == 1:
            # Polynomial representation (spontaneous fission)
            coefficients = np.asarray(self._get_list_record(onlyList=True))
            self.fission['nu']['prompt'] = Polynomial(coefficients)
        elif LNU == 2:
            # Tabulated values of nu
            params, self.fission['nu']['prompt'] = self._get_tab1_record()

        # Skip SEND record
        self._fh.readline()

    def _read_fission_energy(self):
        self._print_info(1, 458)
        er = self.fission['energy_release']

        # Skip HEAD record
        self._get_head_record()

        # Read LIST record containing components of fission energy release (or
        # coefficients)
        items, values = self._get_list_record()
        NPLY = items[3]
        er['order'] = NPLY

        values = np.asarray(values)
        values.shape = (NPLY + 1, 18)
        er['fission_products'] = np.vstack((values[:,0], values[:,1]))
        er['prompt_neutrons'] = np.vstack((values[:,2], values[:,3]))
        er['delayed_neutrons'] = np.vstack((values[:,4], values[:,5]))
        er['prompt_gammas'] = np.vstack((values[:,6], values[:,7]))
        er['delayed_gammas'] = np.vstack((values[:,8], values[:,9]))
        er['delayed_betas'] = np.vstack((values[:,10], values[:,11]))
        er['neutrinos'] = np.vstack((values[:,12], values[:,13]))
        er['total_less_neutrinos'] = np.vstack((values[:,14], values[:,15]))
        er['total'] = np.vstack((values[:,16], values[:,17]))

        # Skip SEND record
        self._fh.readline()

    def _read_reaction_xs(self, MT):
        self._print_info(3, MT)

        # Get Reaction instance
        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(3)

        # Read HEAD record with ZA and atomic mass ratio
        items = self._get_head_record()

        # Read TAB1 record with reaction cross section
        params, rx.xs = self._get_tab1_record()
        rx.Q_mass_difference = params[0]
        rx.Q_reaction = params[1]
        rx.complex_breakup_flag = params[3]

        # Skip SEND record
        self._fh.readline()

    def _read_angular_distribution(self, MT):
        # Find energy distribution
        self._print_info(4, MT)

        # Get Reaction instance
        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(4)

        # Read HEAD record
        items = self._get_head_record()
        ltt = items[3]

        # Read CONT record
        items =self._get_cont_record()
        li = items[2]
        rx.center_of_mass = (items[3] == 2)

        if ltt == 0 and li == 1:
            # Purely isotropic
            energy = np.array([0, self.info['energy_max']])
            mu = [Uniform(-1., 1.), Uniform(-1., 1.)]

        elif ltt == 1 and li == 0:
            # Legendre polynomial coefficients
            tab2 = self._get_tab2_record()
            n_energy = tab2.params[5]

            energy = np.zeros(n_energy)
            mu = []
            for i in range(n_energy):
                items, al = self._get_list_record()
                temperature = items[0]
                energy[i] = items[1]
                coefficients = np.asarray([1.0] + al)
                mu.append(Legendre(coefficients))

        elif ltt == 2 and li == 0:
            # Tabulated probability distribution
            tab2 = self._get_tab2_record()
            n_energy = tab2.params[5]

            energy = np.zeros(n_energy)
            mu = []
            for i in range(n_energy):
                params, f = self._get_tab1_record()
                temperature = params[0]
                energy[i] = params[1]
                if f.n_regions > 1:
                    raise NotImplementedError('Angular distribution with multiple '
                                              'interpolation regions not supported.')
                mu.append(Tabular(f.x, f.y, interpolation_scheme[f.interpolation[0]]))

        elif ltt == 3 and li == 0:
            # Legendre for low energies / tabulated for high energies
            tab2_legendre = self._get_tab2_record()
            n_energy_legendre = tab2_legendre.params[5]

            energy_legendre = np.zeros(n_energy_legendre)
            mu = []
            for i in range(n_energy_legendre):
                items, al = self._get_list_record()
                temperature = items[0]
                energy_legendre[i] = items[1]
                coefficients = np.asarray([1.0] + al)
                mu.append(Legendre(coefficients))

            tab2_tabulated = self._get_tab2_record()
            n_energy_tabulated = tab2_tabulated.params[5]

            energy_tabulated = np.zeros(n_energy_tabulated)
            for i in range(n_energy_tabulated):
                params, f = self._get_tab1_record()
                temperature = params[0]
                energy_tabulated[i] = params[1]
                if f.n_regions > 1:
                    raise NotImplementedError('Angular distribution with multiple '
                                              'interpolation regions not supported.')
                mu.append(Tabular(f.x, f.y, interpolation_scheme[f.interpolation[0]]))

            energy = np.concatenate((energy_legendre, energy_tabulated))

        rx.angular_distribution = AngleDistribution(energy, mu)

    def _read_energy_distribution(self, MT):
        # Find energy distribution
        self._print_info(5, MT)

        # Get Reaction instance
        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(5)

        # Read HEAD record
        items = self._get_head_record()
        nk = items[4]

        for i in range(nk):
            # Read TAB1 record for p(E)
            params, applicability = self._get_tab1_record()
            lf = params[3]
            if lf == 1:
                tab2 = self._get_tab2_record()
                n_energies = tab2.params[5]

                energy = np.zeros(n_energies)
                pdf = []
                for j in range(n_energies):
                    params, func = self._get_tab1_record()
                    energy[j] = params[1]
                    pdf.append(func)
                edist = ArbitraryTabulated(energy, pdf)
            elif lf == 5:
                # General evaporation spectrum
                u = params[0]
                params, theta = self._get_tab1_record()
                params, g = self._get_tab1_record()
                edist = GeneralEvaporation(theta, g, u)
            elif lf == 7:
                # Simple Maxwellian fission spectrum
                u = params[0]
                params, theta = self._get_tab1_record()
                edist = MaxwellEnergy(theta, u)
            elif lf == 9:
                # Evaporation spectrum
                u = params[0]
                params, theta = self._get_tab1_record()
                edist = Evaporation(theta, u)
            elif lf == 11:
                # Energy-dependent Watt spectrum
                u = params[0]
                params, a = self._get_tab1_record()
                params, b = self._get_tab1_record()
                edist = WattEnergy(a, b, u)
            elif lf == 12:
                # Energy-dependent fission neutron spectrum (Madland-Nix)
                params, tm = self._get_tab1_record()
                efl, efh = params[0:2]
                edist = MadlandNix(efl, efh, tm)

            edist.applicability = applicability
            rx.energy_distribution.append(edist)

    def _read_product_energy_angle(self, MT):
        # Find distribution
        self._print_info(6, MT)

        # Get Reaction instance
        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(6)

        # Read HEAD record
        items = self._get_head_record()
        rx.reference_frame = {1: 'laboratory', 2: 'center-of-mass',
                               3: 'light-heavy'}[items[3]]
        n_products = items[4]

        for i in range(n_products):
            product = {}
            rx.product_distribution.append(product)

            # Read TAB1 record for product yield
            params, product['yield'] = self._get_tab1_record()
            product['za'] = params[0]
            product['mass'] = params[1]
            product['lip'] = params[2]
            product['law'] = params[3]

            if product['law'] == 1:
                # Continuum energy-angle distribution
                tab2 = self._get_tab2_record()
                product['lang'] = tab2.params[2]
                product['lep'] = tab2.params[3]
                ne = tab2.params[5]
                product['energy'] = np.zeros(ne)
                product['n_discrete_energies'] = np.zeros(ne)
                product['energy_out'] = []
                product['b'] = []
                for i in range(ne):
                    items, values = self._get_list_record()
                    product['energy'][i] = items[1]
                    product['n_discrete_energies'][i] = items[2]
                    n_angle = items[3]
                    n_energy_out = items[5]
                    values = np.array(values)
                    values.shape = (n_energy_out, n_angle + 2)
                    product['energy_out'].append(values[:,0])
                    product['b'].append(values[:,1:])

            elif product['law'] == 2:
                # Discrete two-body scattering
                tab2 = self._get_tab2_record()
                ne = tab2.params[5]
                product['energy'] = np.zeros(ne)
                product['lang'] = np.zeros(ne, dtype=int)
                product['Al'] = []
                for i in range(ne):
                    items, values = self._get_list_record()
                    product['energy'][i] = items[1]
                    product['lang'][i] = items[2]
                    product['Al'].append(np.asarray(values))

            elif product['law'] == 5:
                # Charged particle elastic scattering
                tab2 = self._get_tab2_record()
                product['spin'] = tab2.params[0]
                product['identical'] = (tab2.params[2] == 1)
                ne = tab2.params[5]
                product['energies'] = np.zeros(ne)
                product['ltp'] = np.zeros(ne, dtype=int)
                product['coefficients'] = []
                for i in range(ne):
                    items, values = self._get_list_record()
                    product['energies'][i] = items[1]
                    product['lpt'][i] = items[2]
                    product['coefficients'].append(values)

            elif product['law'] == 6:
                # N-body phase-space distribution
                items = self._get_cont_record()
                product['total_mass'] = items[0]
                product['n_particles'] = items[5]

            elif product['law'] == 7:
                # Laboratory energy-angle distribution
                tab2 = self._get_tab2_record()
                ne = tab2.params[5]
                product['energies'] = np.zeros(ne)
                product['mu'] = []
                product['distribution'] = []
                for i in range(ne):
                    tab2mu = self._get_tab2_record()
                    product['energies'][i] = tab2mu.params[1]
                    nmu = tab2mu.params[5]
                    mu = np.zeros(nmu)
                    dists = []
                    for j in range(nmu):
                        params, f = self._get_tab1_record()
                        mu[j] = params[1]
                        dists.append(f)
                    product['mu'].append(mu)
                    product['distribution'].append(dists)

    def _read_delayed_photon(self):
        self._print_info(1, 460)
        dp = self.fission['delayed_photon']

        # Determine whether discrete or continuous representation
        items = self._get_head_record()
        LO = items[2]
        NG = items[4]

        # Discrete representation
        if LO == 1:
            dp['form'] = 'discrete'

            # Initialize lists for energies of photons and time dependence of
            # photon multiplicity
            dp['energy'] = np.zeros(NG)
            dp['multiplicity'] = []
            for i in range(NG):
                # Read TAB1 record with multiplicity as function of time
                params, mult = self._get_tab1_record()
                dp['multiplicity'].append(mult)

                # Determine energy
                dp['energy'][i] = params[0]

        # Continuous representation
        elif LO == 2:
            # Determine decay constant and number of precursor families
            dp['form'] = 'continuous'
            dp['decay_constant'] = self._get_list_record(onlyList=True)

    def _read_resonances(self):
        self._print_info(2, 151)
        res = self.resonances

        # Determine whether discrete or continuous representation
        items = self._get_head_record()
        NIS = items[4] # Number of isotopes
        res['isotopes'] = []

        for iso in range(NIS):
            # Create dictionary for this isotope
            isotope = {}
            res['isotopes'].append(isotope)

            items = self._get_cont_record()
            isotope['abundance'] = items[1]
            LFW = items[3] # average fission width flag
            NER = items[4] # number of resonance energy ranges

            isotope['ranges'] = []

            for j in range(NER):
                items = self._get_cont_record()
                emin, emax = items[0:2]  # min/max energies of range
                resonance_flag = items[2]  # flag for resolved (1)/unresolved (2)
                resonance_formalism = items[3]  # resonance formalism
                nro = items[4]  # flag for energy dependence of scattering radius
                naps = items[5]  # flag controlling use of channel/scattering radius

                if resonance_flag == 0 and nro == 0:
                    # Only scattering radius specified
                    erange = ScatteringRadius(emin, emax, nro, naps)
                    items = self._get_cont_record()
                    erange.spin = items[0]
                    erange.scattering_radius = items[1]

                elif resonance_flag == 1:
                    # resolved resonance region
                    erange = _formalisms[resonance_formalism](
                        emin, emax, nro, naps)
                    erange.read(self)
                    if NIS == 1:
                        res['resolved'] = erange

                elif resonance_flag == 2:
                    # unresolved resonance region
                    erange = Unresolved(emin, emax, nro, naps)
                    erange.fission_widths = (LFW == 1)
                    erange.LRF = resonance_formalism
                    erange.read(self)
                    if NIS == 1:
                        res['unresolved'] = erange

                erange.material = self
                isotope['ranges'].append(erange)

    def _read_thermal_elastic(self):
        self._print_info(7, 2)
        elast = self.thermal_elastic

        # Get head record
        items = self._get_head_record()
        LTHR = items[2]  # coherent/incoherent flag
        elast['S'] = {}

        if LTHR == 1:
            elast['type'] = 'coherent'
            params, sdata = self._get_tab1_record()
            temperature = params[0]
            bragg_edges = sdata.x
            LT = params[2]
            elast['S'][temperature] = CoherentElastic(bragg_edges, sdata.y)

            for t in range(LT):
                params, sdata = self._get_list_record()
                temperature = params[0]
                LT = params[2]
                elast['S'][temperature] = CoherentElastic(bragg_edges, sdata)

        elif LTHR == 2:
            elast['type'] = 'incoherent'
            params, wt = self._get_tab1_record()
            elast['bound_xs'] = params[0]
            elast['debye_waller'] = wt

    def _read_thermal_inelastic(self):
        self._print_info(7, 4)
        inel = self.thermal_inelastic

        # Get head record
        items = self._get_head_record()
        inel['temperature_used'] = 'actual' if items[3] == 0 else '0.0253 eV'  # Temperature flag
        inel['symmetric'] = (items[4] == 0)  # Symmetry flag
        header, B = self._get_list_record()
        inel['ln(S)'] = (header[2] == 1)
        inel['num_non_principal'] = header[5]
        inel['B'] = B
        if B[0] != 0.0:
            tab2 = self._get_tab2_record()
            n_beta = tab2.NBT[0]
            for i_beta in range(n_beta):
                #Read record for first temperature (always present)
                params, sab0 = self._get_tab1_record()
                n_temps = params[2] + 1

                # Create arrays on first pass through -- note that alphas and
                # temperatures only need to be stored on first beta
                if i_beta == 0:
                    alpha_values = sab0.x
                    n_alpha = alpha_values.shape[0]
                    beta_values = np.zeros(n_beta)
                    temp_values = np.zeros(n_temps)
                    temp_values[0] = params[0]
                    sab_values = np.zeros((n_alpha, n_beta, n_temps), order='F')

                # Store beta and S(a,b,0) for first beta
                beta_values[i_beta] = params[1]
                sab_values[:, i_beta, 0] = sab0.y

                for i_temp in range(1, n_temps):
                    # Read records for all the other temperatures
                    params, sabt = self._get_list_record()
                    if i_beta == 0:
                        temp_values[i_temp] = params[0]
                    sab_values[:, i_beta, i_temp] = sabt

            # Store arrays in dictionary
            inel['scattering_law'] = sab_values
            inel['alpha'] = alpha_values
            inel['beta'] = beta_values
            inel['temperature'] = temp_values

        params, teff = self._get_tab1_record()
        inel['teff'] = teff

    def _read_independent_yield(self):
        self._print_info(8, 454)
        iyield = self.fission['yield_independent']

        # Initialize energies and yield dictionary
        iyield['energies'] = []
        iyield['data'] = {}
        iyield['interp'] = []

        items = self._get_head_record()
        LE = items[2] - 1  # Determine energy-dependence

        for i in range(LE + 1):
            items, itemList = self._get_list_record()
            E = items[0]  # Incident particle energy
            iyield['energies'].append(E)
            NFP = items[5]  # Number of fission product nuclide states
            if i > 0:
                iyield['interp'].append(items[2]) # Interpolation scheme

            # Get data for each yield
            iyield['data'][E] = {}
            iyield['data'][E]['zafp'] = [int(i) for i in itemList[0::4]] # ZA for fission products
            iyield['data'][E]['fps'] = itemList[1::4] # State designator
            iyield['data'][E]['yi'] = zip(itemList[2::4],itemList[3::4]) # Independent yield

        # Skip SEND record
        self._fh.readline()

    def _read_cumulative_yield(self):
        self._print_info(8, 459)
        cyield = self.fission['yield_cumulative']

        # Initialize energies and yield dictionary
        cyield['energies'] = []
        cyield['data'] = {}
        cyield['interp'] = []

        items = self._get_head_record()
        LE = items[2] - 1  # Determine energy-dependence

        for i in range(LE + 1):
            items, itemList = self._get_list_record()
            E = items[0]  # Incident particle energy
            cyield['energies'].append(E)
            NFP = items[5]  # Number of fission product nuclide states
            if i > 0:
                cyield['interp'].append(items[2]) # Interpolation scheme

            # Get data for each yield
            cyield['data'][E] = {}
            cyield['data'][E]['zafp'] = [int(i) for i in itemList[0::4]] # ZA for fission products
            cyield['data'][E]['fps'] = itemList[1::4] # State designator
            cyield['data'][E]['yc'] = zip(itemList[2::4],itemList[3::4]) # Cumulative yield

        # Skip SEND record
        self._fh.readline()

    def _read_decay(self):
        self._print_info(8, 457)
        decay = self.decay

        # Get head record
        items = self._get_head_record()
        decay['ZA'] = items[0]  # ZA identifier
        decay['awr'] = items[1]  # AWR
        decay['state']= items[2]  # State of the original nuclide
        decay['isomeric_state'] = items[3]  # Isomeric state for the original nuclide
        decay['stable'] = (items[4] == 1)  # Nucleus stability flag

        # Determine if radioactive/stable
        if not decay['stable']:
            NSP = items[5]  # Number of radiation types

            # Half-life and decay energies
            items, itemList = self._get_list_record()
            decay['half_life'] = (items[0], items[1])
            decay['NC'] = items[4]//2
            decay['energies'] = zip(itemList[0::2], itemList[1::2])

            # Decay mode information
            items, itemList = self._get_list_record()
            decay['spin'] = items[0]  # Spin of the nuclide
            decay['parity'] = items[1]  # Parity of the nuclide
            NDK = items[5]  # Number of decay modes

            # Decay type (beta, gamma, etc.)
            decay['modes'] = []
            for i in range(NDK):
                mode = {}
                mode['type'] = radiation_type(itemList[6*i])
                mode['isomeric_state'] = itemList[6*i + 1]
                mode['energy'] = tuple(itemList[6*i + 2:6*i + 4])
                mode['branching_ratio'] = tuple(itemList[6*i + 4:6*(i + 1)])
                decay['modes'].append(mode)

            discrete_type = {0.0: None, 1.0: 'allowed', 2.0: 'first-forbidden',
                             3.0: 'second-forbidden'}

            # Read spectra
            decay['spectra'] = {}
            for i in range(NSP):
                spectrum = {}

                items, itemList = self._get_list_record()
                # Decay radiation type
                spectrum['type'] = radiation_type(items[1])
                # Continuous spectrum flag
                spectrum['continuous_flag'] = {0: 'discrete', 1: 'continuous',
                                               2: 'both'}[items[2]]
                spectrum['discrete_normalization'] = tuple(itemList[0:2])
                spectrum['energy_average'] = tuple(itemList[2:4])
                spectrum['continuous_normalization'] = tuple(itemList[4:6])

                NER = items[5]  # Number of tabulated discrete energies

                if not spectrum['continuous_flag'] == 'continuous':
                    # Information about discrete spectrum
                    spectrum['discrete'] = []
                    for j in range(NER):
                        items, itemList = self._get_list_record()
                        di = {}
                        di['energy'] = tuple(items[0:2])
                        di['from_mode'] = radiation_type(itemList[0])
                        di['type'] = discrete_type[itemList[1]]
                        di['intensity'] = tuple(itemList[2:4])
                        if spectrum['type'] == 'ec/beta+':
                            di['positron_intensity'] = tuple(itemList[4:6])
                        elif spectrum['type'] == 'gamma':
                            di['internal_pair'] = tuple(itemList[4:6])
                            di['total_internal_conversion'] = tuple(itemList[6:8])
                            di['k_shell_conversion'] = tuple(itemList[8:10])
                            di['l_shell_conversion'] = tuple(itemList[10:12])
                        spectrum['discrete'].append(di)

                if not spectrum['continuous_flag'] == 'discrete':
                    # Read continuous spectrum
                    ci = {}
                    params, ci['probability'] = self._get_tab1_record()
                    ci['type'] = radiation_type(params[0])

                    # Read covariance (Ek, Fk) table
                    LCOV = params[3]
                    if LCOV != 0:
                        items, itemList = self._get_list_record()
                        ci['covariance_lb'] = items[3]
                        ci['covariance'] = zip(itemList[0::2], itemList[1::2])

                    spectrum['continuous'] = ci

                # Add spectrum to dictionary
                decay['spectra'][spectrum['type']] = spectrum

        else:
            items, itemList = self._get_list_record()
            items, itemList = self._get_list_record()
            decay['spin'] = items[0]
            decay['parity'] = items[1]

        # Skip SEND record
        self._fh.readline()

    def _read_radioactive_nuclide(self, MT):
        self._print_info(8, MT)

        # Get Reaction instance
        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(8)

        # Get head record
        items = self._get_head_record()
        NS = items[4]
        complete_chain = (items[5] == 0)

        if complete_chain:
            # If complete chain is specified here rather than in MF=8, get all
            # levels of radionuclide and corresponding decay modes
            rx.radionuclide_production = []
            for i in range(NS):
                items, values = self._get_list_record()
                radionuclide = {}
                radionuclide['za'] = items[0]
                radionuclide['excitation_energy'] = items[1]
                radionuclide['mf_multiplicity'] = items[2]
                radionuclide['level'] = items[3]
                radionuclide['modes'] = []
                for j in range(items[4]//6):
                    mode = {}
                    mode['half_life'] = values[6*j]
                    mode['type'] = radiation_type(values[6*j + 1])
                    mode['za'] = values[6*j + 2]
                    mode['branching_ratio'] = values[6*j + 3]
                    mode['endpoint_energy'] = values[6*j + 4]
                    mode['chain_terminator'] = values[6*j + 5]
                    radionuclide['modes'].append(mode)
                rx.radionuclide_production.append(radionuclide)
        else:
            items = self._get_cont_record()
            rx.radionuclide_production = radionuclide = {}
            radionuclide['za'] = items[0]
            radionuclide['excitation_energy'] = items[1]
            radionuclide['mf_multiplicity'] = items[2]
            radionuclide['level'] = items[3]

    def _read_multiplicity(self, MT):
        self._print_info(9, MT)

        # Get Reaction instance
        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(9)

        # Get head record
        items = self._get_head_record()
        NS = items[4]  # Number of final states

        for i in range(NS):
            params, state = self._get_tab1_record()
            QM = params[0] # Mass difference Q value (eV)
            QI = params[1] # Reaction Q value (eV)
            IZAP = params[2] # 1000Z + A
            LFS = params[3] # Level number of the nuclide
            rx.multiplicities[LFS] = {'QM': QM, 'QI': QI, 'ZA': IZAP,
                                       'values': state}

    def _read_production_xs(self, MT):
        self._print_info(10, MT)

        # Get Reaction instance
        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(10)

        # Get head record
        items = self._get_head_record()
        NS = items[4]  # Number of final states

        for i in range(NS):
            params, state = self._get_tab1_record()
            QM = params[0]  # Mass difference Q value (eV)
            QI = params[1]  # Reaction Q value (eV)
            IZAP = params[2]  # 1000Z + A
            LFS = params[3]  # Level number of the nuclide
            rx.production[LFS] = {'QM': QM, 'QI': QI, 'ZA': IZAP,
                                   'values': state}

    def _read_photon_production_yield(self, MT):
        self._print_info(12, MT)

        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(12)
        rx.photon_production['yield'] = ppyield = {}

        # Determine option
        items = self._get_head_record()
        option = items[2]

        if option == 1:
            # Multiplicities given
            ppyield['type'] = 'multiplicity'
            n_discrete_photon = items[4]
            if n_discrete_photon > 1:
                items, ppyield['total'] = self._get_tab1_record()
            ppyield['discrete'] = []
            for k in range(n_discrete_photon):
                y = {}
                items, y['yield'] = self._get_tab1_record()
                y['energy_photon'] = items[0]
                y['energy_level'] = items[1]
                y['lp'] = items[2]
                y['law'] = items[3]
                ppyield['discrete'].append(y)

        elif option == 2:
            # Transition probability arrays given
            ppyield['type'] = 'transition'
            ppyield['transition'] = transition = {}

            # Determine whether simple (LG=1) or complex (LG=2) transitions
            lg = items[3]

            # Get transition data
            items, values = self._get_list_record()
            transition['energy_start'] = items[0]
            transition['energies'] = np.array(values[::lg + 1])
            transition['direct_probability'] = np.array(values[1::lg + 1])
            if lg == 2:
                # Complex case
                transition['conditional_probability'] = np.array(
                    values[2::lg + 1])

    def _read_photon_production_xs(self, MT):
        self._print_info(13, MT)

        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(13)
        rx.photon_production['xs'] = ppxs = {}

        # Determine option
        items = self._get_head_record()
        n_discrete_photon = items[4]
        if n_discrete_photon > 1:
            items, ppxs['total'] = self._get_tab1_record()
        ppxs['discrete'] = []
        for k in range(n_discrete_photon):
            xs = {}
            items, xs['xs'] = self._get_tab1_record()
            xs['energy_photon'] = items[0]
            xs['energy_level'] = items[1]
            xs['lp'] = items[2]
            xs['law'] = items[3]
            ppxs['discrete'].append(xs)

    def _read_photon_angular_distribution(self, MT):
        self._print_info(14, MT)

        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(14)
        rx.photon_production['angular_distribution'] = ppad = {}

        # Determine format for angular distributions
        items = self._get_head_record()
        ppad['isotropic'] = (items[2] == 1)

        if not ppad['isotropic']:
            ltt = items[3]
            n_discrete_photon = items[4]
            n_isotropic = items[5]
            ppad['discrete'] = []
            for i in range(n_isotropic):
                adist = AngularDistribution()
                adist.type = 'isotropic'

                items = self._get_cont_record()
                adist.energy_photon = items[0]
                adist.energy_level = items[1]
                ppad['discrete'].append(adist)

            if ltt == 1:
                # Legendre polynomial coefficients
                for i in range(n_isotropic, n_discrete_photon):
                    adist = AngularDistribution()
                    adist.type = 'legendre'

                    adist.tab2 = self._get_tab2_record()
                    adist.energy_photon = adist.tab2.params[0]
                    adist.energy_level = adist.tab2.params[1]
                    n_energy = adist.tab2.params[5]

                    adist.energy = np.zeros(n_energy)
                    adist.probability = []
                    for i in range(n_energy):
                        items, al = self._get_list_record()
                        adist.energy[i] = items[1]
                        coefficients = np.asarray([1.0] + al)
                        for i in range(len(coefficients)):
                            coefficients[i] *= (2.*i + 1.)/2.
                        adist.probability.append(Legendre(coefficients))
                    ppad['discrete'].append(adist)

            elif ltt == 2:
                # Tabulated probability distribution
                for i in range(n_isotropic, n_discrete_photon):
                    adist.type = 'tabulated'

                    adist.tab2 = self._get_tab2_record()
                    adist.energy_photon = adist.tab2.params[0]
                    adist.energy_level = adist.tab2.params[1]
                    n_energy = adist.tab2.params[5]

                    adist.energy = np.zeros(n_energy)
                    adist.probability = []
                    for i in range(n_energy):
                        params, f = self._get_tab1_record()
                        adist.energy[i] = params[1]
                        adist.probability.append(f)
                    ppad['discrete'].append(adist)

    def _read_photon_energy_distribution(self, MT):
        self._print_info(15, MT)

        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(15)
        rx.photon_production['energy_distribution'] = pped = []

        # Read HEAD record
        items = self._get_head_record()
        nc = items[4]

        for i in range(nc):
            # Read TAB1 record for p(E)
            params, applicability = self._get_tab1_record()
            lf = params[3]
            if lf == 1:
                # Arbitrary tabulated function -- only format currently
                # available in ENDF-102 for photon continuum
                tab2 = self._get_tab2_record()
                n_energies = tab2.params[5]

                energy = np.zeros(n_energies)
                pdf = []
                for j in range(n_energies):
                    params, func = self._get_tab1_record()
                    energy[j] = params[1]
                    pdf.append(func)
                edist = ArbitraryTabulated(energy, pdf)

            edist.applicability = applicability
            pped.append(edist)

    def _read_photon_interaction(self, MT):
        self._print_info(23, MT)

        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(23)

        # Skip HEAD record
        self._get_head_record()

        # Read cross section
        params, rx.xs = self._get_tab1_record()
        if MT >= 534 and MT <= 599:
            rx.subshell_binding_energy = params[0]
        if MT >= 534 and MT <= 572:
            rx.fluorescence_yield = params[1]

        # Skip SEND record
        self._fh.readline()

    def _read_electron_products(self, MT):
        self._print_info(26, MT)

        if MT not in self.reactions:
            self.reactions[MT] = Reaction(MT)
        rx = self.reactions[MT]
        rx.files.append(26)

        # Read HEAD record
        items = self._get_head_record()
        n_products = items[4]

        for i in range(n_products):
            product = {}
            rx.products.append(product)

            # Read TAB1 record for product yield
            params, product['yield'] = self._get_tab1_record()
            product['za'] = params[0]
            product['law'] = params[3]

            if product['law'] == 1:
                # Continuum energy-angle distribution
                tab2 = self._get_tab2_record()
                product['lang'] = tab2.params[2]
                product['lep'] = tab2.params[3]
                ne = tab2.params[5]
                product['energy'] = np.zeros(ne)
                product['n_discrete_energies'] = np.zeros(ne)
                product['energy_out'] = []
                product['b'] = []
                for i in range(ne):
                    items, values = self._get_list_record()
                    product['energy'][i] = items[1]
                    product['n_discrete_energies'][i] = items[2]
                    n_angle = items[3]
                    n_energy_out = items[5]
                    values = np.array(values)
                    values.shape = (n_energy_out, n_angle + 2)
                    product['energy_out'].append(values[:,0])
                    product['b'].append(values[:,1:])

            elif product['law'] == 2:
                # Discrete two-body scattering
                tab2 = self._get_tab2_record()
                ne = tab2.params[5]
                product['energy'] = np.zeros(ne)
                product['lang'] = np.zeros(ne, dtype=int)
                product['Al'] = []
                for i in range(ne):
                    items, values = self._get_list_record()
                    product['energy'][i] = items[1]
                    product['lang'][i] = items[2]
                    product['Al'].append(np.asarray(values))

            elif product['law'] == 8:
                # Energy transfer for excitation
                params, product['energy_transfer'] = self._get_tab1_record()

    def _read_scattering_functions(self, MT):
        self._print_info(27, MT)

        # Skip HEAD record
        self._get_head_record()

        # Get scattering function
        params, func = self._get_tab1_record()

        # Store in appropriate place
        if MT in (502, 504):
            rx = self.reactions[MT]
            rx.scattering_factor = func
        elif MT == 505:
            rx = self.reactions[502]
            rx.anomalous_scattering_imaginary = func
        elif MT == 506:
            rx = self.reactions[502]
            rx.anomalous_scattering_real = func

        # Skip SEND record
        self._fh.readline()

    def _read_atomic_relaxation(self):
        self._print_info(28, 533)

        # Read HEAD record
        params = self._get_head_record()
        n_subshells = params[4]

        # Read list of data
        subshells = {1: 'K', 2: 'L1', 3: 'L2', 4: 'L3', 5: 'M1',
                     6: 'M2', 7: 'M3', 8: 'M4', 9: 'M5', 10: 'N1',
                     11: 'N2', 12: 'N3', 13: 'N4', 14: 'N5', 15: 'N6',
                     16: 'N7', 17: 'O1', 18: 'O2', 19: 'O3', 20: 'O4',
                     21: 'O5', 22: 'O6', 23: 'O7', 24: 'O8', 25: 'O9',
                     26: 'P1', 27: 'P2', 28: 'P3', 29: 'P4', 30: 'P5',
                     31: 'P6', 32: 'P7', 33: 'P8', 34: 'P9', 35: 'P10',
                     36: 'P11', 37: 'Q1', 38: 'Q2', 39: 'Q3', 0: None}
        for i in range(n_subshells):
             params, list_items = self._get_list_record()
             subi = subshells[int(params[0])]
             n_transitions = int(params[5])
             ebi = list_items[0]
             eln = list_items[1]
             data = {'binding_energy': ebi, 'number_electrons': eln, 'transitions': []}
             for j in range(n_transitions):
                 subj = subshells[int(list_items[6*(j+1)])]
                 subk = subshells[int(list_items[6*(j+1) + 1])]
                 etr = list_items[6*(j+1) + 2]
                 ftr = list_items[6*(j+1) + 3]
                 data['transitions'].append((subj, subk, etr, ftr))
             self.atomic_relaxation[subi] = data

        # Skip SEND record
        self._fh.readline()

    def _get_text_record(self, line=None):
        if not line:
            line = self._fh.readline()
        if self._veryverbose:
            print('Get TEXT record')
        return line[0:66]

    def _get_cont_record(self, line=None, skipC=False):
        if self._veryverbose:
            print('Get CONT record')
        if not line:
            line = self._fh.readline()
        if skipC:
            C1 = None
            C2 = None
        else:
            C1 = endftod(line[:11])
            C2 = endftod(line[11:22])
        L1 = int(line[22:33])
        L2 = int(line[33:44])
        N1 = int(line[44:55])
        N2 = int(line[55:66])
        return [C1, C2, L1, L2, N1, N2]

    def _get_head_record(self, line=None):
        if not line:
            line = self._fh.readline()
        if self._veryverbose:
            print('Get HEAD record')
        ZA = int(endftod(line[:11]))
        AWR = endftod(line[11:22])
        L1 = int(line[22:33])
        L2 = int(line[33:44])
        N1 = int(line[44:55])
        N2 = int(line[55:66])
        return [ZA, AWR, L1, L2, N1, N2]

    def _get_list_record(self, onlyList=False):
        # determine how many items are in list
        if self._veryverbose:
            print('Get LIST record')
        items = self._get_cont_record()
        NPL = items[4]

        # read items
        itemsList = []
        m = 0
        for i in range((NPL-1)//6 + 1):
            line = self._fh.readline()
            toRead = min(6, NPL-m)
            for j in range(toRead):
                val = endftod(line[0:11])
                itemsList.append(val)
                line = line[11:]
            m = m + toRead
        if onlyList:
            return itemsList
        else:
            return (items, itemsList)

    def _get_tab1_record(self):
        if self._veryverbose:
            print('Get TAB1 record')

        # Determine how many interpolation regions and total points there are
        line = self._fh.readline()
        C1 = endftod(line[:11])
        C2 = endftod(line[11:22])
        L1 = int(line[22:33])
        L2 = int(line[33:44])
        n_regions = int(line[44:55])
        n_pairs = int(line[55:66])
        params = [C1, C2, L1, L2]

        # Read the interpolation region data, namely NBT and INT
        nbt = np.zeros(n_regions)
        interp = np.zeros(n_regions)
        m = 0
        for i in range((n_regions - 1)//3 + 1):
            line = self._fh.readline()
            toRead = min(3, n_regions - m)
            for j in range(toRead):
                nbt[m] = int(line[0:11])
                interp[m] = int(line[11:22])
                line = line[22:]
                m += 1

        # Read tabulated pairs x(n) and y(n)
        x = np.zeros(n_pairs)
        y = np.zeros(n_pairs)
        m = 0
        for i in range((n_pairs - 1)//3 + 1):
            line = self._fh.readline()
            toRead = min(3, n_pairs - m)
            for j in range(toRead):
                x[m] = endftod(line[:11])
                y[m] = endftod(line[11:22])
                line = line[22:]
                m += 1

        return params, Tabulated1D(x, y, nbt, interp)

    def _get_tab2_record(self):
        if self._veryverbose:
            print('Get TAB2 record')
        r = ENDFTab2Record()
        r.read(self._fh)
        return r


    def _print_info(self, MF, MT):
        if self._verbose:
            if MT in reaction_name:
                print('Reading MF={0}, MT={1} {2}'.format(MF, MT, reaction_name[MT]))
            else:
                print('Reading MF={0}, MT={1}'.format(MF, MT))

    def __repr__(self):
        try:
            name = self.target['zsymam'].replace(' ', '')
            library = self.info['library']
        except:
            name = 'Undetermined'
            library = 'None'
        return '<Evaluation: {0}, {1}>'.format(name, library)




class ENDFTab2Record(object):
    def __init__(self):
        self.NBT = []
        self.INT = []

    def read(self, fh):
        # Determine how many interpolation regions and total points there are
        line = fh.readline()
        C1 = endftod(line[:11])
        C2 = endftod(line[11:22])
        L1 = int(line[22:33])
        L2 = int(line[33:44])
        NR = int(line[44:55])
        NZ = int(line[55:66])
        self.params = [C1, C2, L1, L2, NR, NZ]

        # Read the interpolation region data, namely NBT and INT
        m = 0
        for i in range((NR-1)//3 + 1):
            line = fh.readline()
            toRead = min(3,NR-m)
            for j in range(toRead):
                NBT = int(line[0:11])
                INT = int(line[11:22])
                self.NBT.append(NBT)
                self.INT.append(INT)
                line = line[22:]
            m = m + toRead


class AngularDistribution(object):
    def __init__(self):
        pass


class Reaction(object):
    """Data for a single reaction including its cross section and secondary
    angle/energy distribution.

    Parameters
    ----------
    MT : int
        The MT number from the ENDF file.

    Attributes
    ----------
    angular_distribution : AngularDistribution
        Angular distribution represented as a tabulated function or as moments
        of a Legendre polynomial expansion from MF=4.
    complex_breakup_flag : int
        Complex breakup flag.
    xs : Tabulated1D
        Tabulated cross section as a function of incident energy from MF=3.
    energy_distribution : list of EnergyDistribution
        List of partial energy distributions for the reaction from MF=5.
    files : list of int
        List of files (MF) that have been read for this reaction
    MT : int
        The MT number from the ENDF file.
    multiplicities : dict
        Multiplicities for production of radioactive nuclides as given in MF=9.
    product_distribution : list of dict
        Secondary energy or correlated energy-angle distribution from MF=6.
    production : dict
        Cross sections for production of radioactive nuclides as given in MF=10.
    Q_mass_difference : float
        Mass difference Q value in eV
    Q_reaction : float
        Reaction Q value in eV
    radionuclide_production : list of dict
        List of radioactive nuclides produced as given in MF=8.
    reference_frame : {'laboratory', 'center-of-mass', 'light-heavy'}
        Indicates what reference frame is used for outgoing energies and
        angles. Only relevant for product energy-angle distributions read in
        MF=6.
    subshell_binding_energy : float
        Subshell binding energy in eV as given in MF=23. Only relevant for
        photoelectric subshell ionization reactions in a photo-atomic
        sublibrary.
    fluorescence_yield : float
        Fluorescence yield (eV/photoionization) as given in MF=23. Only relevant
        for photoelectric subshell ionization reactions in a photo-atomic
        sublibrary.
    products : list of dict
        Secondary photon and electron distributions for electro-atomic
        reactions as given in MF=26.
    scattering_factor : Tabulated1D
        Coherent or incoherent form factor as given in MF=27.
    anomalous_scattering_imaginary : Tabulated1D
        Imaginary component of the anomalous scattering factor as given in
        MF=27.
    anomalous_scattering_real : Tabulated1D
        Real component of the anomalous scattering factor as given in MF=27.

    """

    def __init__(self, MT):
        self.MT = MT
        self.files = []
        self.photon_production = {}
        self.xs = None
        self.Q_mass_difference = None
        self.Q_reaction = None
        self.complex_breakup_flag = None
        self.angular_distribution = None
        self.energy_distribution = []
        self.product_distribution = []
        self.reference_frame = None
        self.radionuclide_production = None
        self.multiplicities = {}
        self.production = {}
        self.subshell_binding_energy = None
        self.fluorescence_yield = None
        self.products = []
        self.scattering_factor = None
        self.anomalous_scattering_imaginary = None
        self.anomalous_scattering_real = None

    def __repr__(self):
        if self.MT in reaction_name:
            return '<ENDF Reaction: MT={0}, {1}>'.format(self.MT, reaction_name[self.MT])
        else:
            return '<ENDF Reaction: MT={0}>'.format(self.MT)


def _wave_number(A, E):
    return 2.196807122623e-3*A/(A + 1)*sqrt(abs(E))

def get_rhos(erange, A, E, l=None):
    k = _wave_number(A, E)

    # Calculate scattering radius
    a = 0.123*(1.008664904*A)**(1.0/3.0) + 0.08

    # Calculate energy-independent channel radius, possibly l-dependent
    if l and erange.apl[l] > 0.:
        AP = erange.apl[l]
    else:
        AP = erange.scattering_radius_ind

    if erange.nro == 0:
        if erange.naps == 0:
            rho = k*a
            rhohat = k*AP
        elif erange.naps == 1:
            rho = rhohat = k*AP
    elif erange.nro == 1:
        APE = erange.scattering_radius(E)
        if erange.naps == 0:
            rho = k*a
            rhohat = k*APE
        elif erange.naps == 1:
            rho = rhohat = k*APE
        elif erange.naps == 2:
            rho = k*AP
            rhohat = k*APE

    return rho, rhohat

def phaseshift(l, rho):
    """Calculate hardsphere phase shift as given in ENDF-102, Equation D.13

    Parameters
    ----------
    l : int
        Angular momentum quantum number
    rho : float
        Product of the wave number and the channel radius

    """

    if l == 0:
        return rho
    elif l == 1:
        return rho - atan(rho)
    elif l == 2:
        return rho - atan(3*rho/(3 - rho**2))
    elif l == 3:
        return rho - atan((15*rho - rho**3)/(15 - 6*rho**2))
    elif l == 4:
        return rho - atan((105*rho - 10*rho**3)/(105 - 45*rho**2 + rho**4))

def penetration_shift(l, rho):
    """Calculate shift and penetration factors as given in ENDF-102, Equations D.11
    and D.12.

    Parameters
    ----------
    l : int
        Angular momentum quantum number
    rho : float
        Product of the wave number and the channel radius

    """

    if l == 0:
        return rho, 0.
    elif l == 1:
        den = 1 + rho**2
        return rho**3/den, -1/den
    elif l == 2:
        den = 9 + 3*rho**2 + rho**4
        return rho**5/den, -(18 + 3*rho**2)/den
    elif l == 3:
        den = 225 + 45*rho**2 + 6*rho**4 + rho**6
        return rho**7/den, -(675 + 90*rho**2 + 6*rho**4)/den
    elif l == 4:
        den = 11025 + 1575*rho**2 + 135*rho**4 + 10*rho**6 + rho**8
        return rho**9/den, -(44100 + 4725*rho**2 + 270*rho**4 + 10*rho**6)/den


class ResonanceRange(object):
    """Resolved resonance formalism data as given in MF=2, MT=151.

    Parameters
    ----------
    emin : float
        Minimum energy of the resolved resonance range in eV
    emax : float
        Maximum energy of the resolved resonance range in eV
    nro : int
        Flag designating energy-dependence of scattering radius (NRO). A value
        of 0 indicates it is energy-independent, a value of 1 indicates it is
        energy-dependent.
    naps : int
        Flag controlling use of the channel and scattering radius (NAPS)

    Attributes
    ----------
    energy_min : float
        Minimum energy of the resolved resonance range in eV
    energy_max : float
        Maximum energy of the resolved resonance range in eV
    nro : int
        Flag designating energy-dependence of scattering radius (NRO). A value
        of 0 indicates it is energy-independent, a value of 1 indicates it is
        energy-dependent.
    naps : int
        Flag controlling use of the channel and scattering radius (NAPS)

    """

    def __init__(self, emin, emax, nro, naps):
        self.energy_min = emin
        self.energy_max = emax
        self.nro = nro
        self.naps = naps
        self._prepared = False

    def read(self, ev):
        raise NotImplementedError

    def reconstruct(self):
        raise NotImplementedError


class MultiLevelBreitWigner(ResonanceRange):
    """Multi-level Breit-Wigner resolved resonance formalism data. This is
    identified by LRF=2 in the ENDF-6 format.

    Parameters
    ----------
    emin : float
        Minimum energy of the resolved resonance range in eV
    emax : float
        Maximum energy of the resolved resonance range in eV
    nro : int
        Flag designating energy-dependence of scattering radius (NRO). A value
        of 0 indicates it is energy-independent, a value of 1 indicates it is
        energy-dependent.
    naps : int
        Flag controlling use of the channel and scattering radius (NAPS)

    Attributes
    ----------
    spin : float
        Spin of the target nucleus
    scattering_radius_ind : float
        Scattering radius in units of 10^-12 cm
    resonances : list of Resonance
        List of resolved resonances
    l_values : ndarray of int
        Neutron orbital angular momentum values
    q : ndarray of float
        Q-value to be added to incident particle's center-of-mass energy to
        determine the channel energy for use in the penetrability factor. Given
        as a function of the l-value.
    competitive : ndarray of bool
        Whether a competitive width is given for each l-value.

    """

    def __init__(self, emin, emax, nro, naps):
        super(MultiLevelBreitWigner, self).__init__(emin, emax, nro, naps)
        self.spin = None
        self.scattering_radius_ind = None
        self.resonances = []
        self.l_values = None
        self.q = None
        self.competitive = None

    def read(self, ev):
        """Read multi-level Breit-Wigner resonance data from an ENDF evaluation.

        Parameters
        ----------
        ev : Evaluation
            ENDF evaluation positioned at the start of a resonance range
            subsection

        """

        # Read energy-dependent scattering radius if present
        if self.nro != 0:
            params, self.scattering_radius = ev._get_tab1_record()

        # Other scatter radius parameters
        items = ev._get_cont_record()
        self.spin = items[0]
        self.scattering_radius_ind = items[1]
        NLS = items[4]  # Number of l-values

        self.l_values = np.zeros(NLS, int)
        self.q = np.zeros(NLS)
        self.competitive = np.zeros(NLS, bool)

        # Read resonance widths, J values, etc
        for l in range(NLS):
            items, values = ev._get_list_record()
            self.q[l] = items[1]
            self.l_values[l] = items[2]
            self.competitive[l] = items[3]
            energy = values[0::6]
            spin = values[1::6]
            GT = values[2::6]
            GN = values[3::6]
            GG = values[4::6]
            GF = values[5::6]

            resonances = []
            for i, E in enumerate(energy):
                resonance = Resonance()
                resonance.energy = energy[i]
                resonance.spin = spin[i]
                resonance.width_total = GT[i]
                resonance.width_neutron = GN[i]
                resonance.width_gamma = GG[i]
                resonance.width_fissionA = GF[i]
                resonances.append(resonance)
            self.resonances.append(resonances)


class SingleLevelBreitWigner(MultiLevelBreitWigner):
    """Single-level Breit-Wigner resolved resonance formalism data. This is
    identified by LRF=1 in the ENDF-6 format.

    Parameters
    ----------
    emin : float
        Minimum energy of the resolved resonance range in eV
    emax : float
        Maximum energy of the resolved resonance range in eV
    nro : int
        Flag designating energy-dependence of scattering radius (NRO). A value
        of 0 indicates it is energy-independent, a value of 1 indicates it is
        energy-dependent.
    naps : int
        Flag controlling use of the channel and scattering radius (NAPS)

    """

    def __init__(self, emin, emax, nro, naps):
        super(SingleLevelBreitWigner, self).__init__(emin, emax, nro, naps)


class ReichMoore(ResonanceRange):
    """Reich-Moore resolved resonance formalism data. This is identified by LRF=3 in
    the ENDF-6 format.

    Parameters
    ----------
    emin : float
        Minimum energy of the resolved resonance range in eV
    emax : float
        Maximum energy of the resolved resonance range in eV
    nro : int
        Flag designating energy-dependence of scattering radius (NRO). A value
        of 0 indicates it is energy-independent, a value of 1 indicates it is
        energy-dependent.
    naps : int
        Flag controlling use of the channel and scattering radius (NAPS)

    Attributes
    ----------
    resonances : dict
        Dictionary whose keys are (l, J) pairs and values are instances are
        lists of Resonance instances.
    spin : float
        Spin of the target nucleus
    scattering_radius_ind : float
        Energy-independent scattering radius in units of 10^-12 cm
    scattering_radius : Tabulated1D
        Scattering radius in units of 10^-12 cm as a function of energy
    LAD : int
        Indicate whether parameters can be used to compute angular distributions
    NLSC : int
        Number of l-values which must be used to converge the calculation
    apl : ndarray of float
        l-dependent scattering radius
    l_values : ndarray of int
        Neutron orbital angular momentum values

    """

    def __init__(self, emin, emax, nro, naps):
        super(ReichMoore, self).__init__(emin, emax, nro, naps)
        self.resonances = {}
        self.spin = None
        self.scattering_radius_ind = None
        self.scattering_radius = None
        self.LAD = None
        self.NLSC = None
        self.apl = None
        self.l_values = None

    def read(self, ev):
        """Read Reich-Moore resonance data from an ENDF evaluation.

        Parameters
        ----------
        ev : Evaluation
            ENDF evaluation positioned at the start of a resonance range
            subsection

        """
        # Read energy-dependent scattering radius if present
        if self.nro != 0:
            params, self.scattering_radius = ev._get_tab1_record()

        # Other scatter radius parameters
        items = ev._get_cont_record()
        self.spin = items[0]
        self.scattering_radius_ind = items[1]
        self.LAD = items[3]  # Flag for angular distribution
        NLS = items[4]  # Number of l-values
        self.NLSC = items[5]  # Number of l-values for convergence

        self.apl = np.zeros(NLS)
        self.l_values = np.zeros(NLS, int)

        # Read resonance widths, J values, etc
        for i in range(NLS):
            items, values = ev._get_list_record()
            self.apl[i] = items[1]
            self.l_values[i] = l = items[2]
            energy = values[0::6]
            spin = values[1::6]
            GN = values[2::6]
            GG = values[3::6]
            GFA = values[4::6]
            GFB = values[5::6]
            for i, E in enumerate(energy):
                resonance = Resonance()
                resonance.energy = energy[i]
                resonance.spin = J = spin[i]
                resonance.width_neutron = GN[i]
                resonance.width_gamma = GG[i]
                resonance.width_fissionA = GFA[i]
                resonance.width_fissionB = GFB[i]

                if not (l, abs(J)) in self.resonances:
                    self.resonances[l, abs(J)] = []
                self.resonances[l, abs(J)].append(resonance)


class AdlerAdler(ResonanceRange):
    """Adler-Adler resolved resonance formalism data. This is identified by LRF=4 in
    the ENDF-6 format.

    Parameters
    ----------
    emin : float
        Minimum energy of the resolved resonance range in eV
    emax : float
        Maximum energy of the resolved resonance range in eV
    nro : int
        Flag designating energy-dependence of scattering radius (NRO). A value
        of 0 indicates it is energy-independent, a value of 1 indicates it is
        energy-dependent.
    naps : int
        Flag controlling use of the channel and scattering radius (NAPS)

    Attributes
    ----------
    spin : float
        Spin of the target nucleus
    scattering_radius_ind : float
        Scattering radius in units of 10^-12 cm
    LI : float
        Flag indicating the kind of parameters
    AT, BT : ndarray
        Background constants for the total cross section
    AC, BC : ndarray
        Background constants for the radiative cross section
    AF, BF : ndarray
        Background constants for the fission cross section
    resonances : list of AdlerResonance
        List of resonances defined by Adler-Adler parameters

    """

    def __init__(self, emin, emax, nro, naps):
        super(AdlerAdler, self).__init__(emin, emax, nro, naps)

    def read(self, ev):
        """Read Adler-Adler resonance data from an ENDF evaluation.

        Parameters
        ----------
        ev : Evaluation
            ENDF evaluation positioned at the start of a resonance range
            subsection

        """
        # Read energy-dependent scattering radius if present
        if self.nro != 0:
            params, self.scattering_radius = ev._get_tab1_record()

        # Other scatter radius parameters
        items = ev._get_cont_record()
        self.spin = items[0]
        self.scattering_radius_ind = items[1]
        NLS = items[4]  # Number of l-values

        # Get AT, BT, AF, BF, AC, BC constants
        items, values = ev._get_list_record()
        self.LI = items[2]
        NX = items[5]
        self.AT = np.asarray(values[:4])
        self.BT = np.asarray(values[4:6])
        if NX == 2:
            self.AC = np.asarray(values[6:10])
            self.BC = np.asarray(values[10:12])
        elif NX == 3:
            self.AF = np.asarray(values[6:10])
            self.BF = np.asarray(values[10:12])
            self.AC = np.asarray(values[12:16])
            self.BC = np.asarray(values[16:18])

        self.resonances = []

        for ls in range(NLS):
            items = ev._get_cont_record()
            l_value = items[2]
            NJS = items[4]
            for j in range(NJS):
                items, values = ev._get_list_record()
                AJ = items[0]
                NLJ = items[5]
                for res in range(NLJ):
                    resonance = AdlerResonance()
                    resonance.L, resonance.J = l_value, AJ
                    resonance.DET = values[12*res]
                    resonance.DWT = values[12*res + 1]
                    resonance.DRT = values[12*res + 2]
                    resonance.DIT = values[12*res + 3]
                    resonance.DEF_ = values[12*res + 4]
                    resonance.DWF = values[12*res + 5]
                    resonance.GRF = values[12*res + 6]
                    resonance.GIF = values[12*res + 7]
                    resonance.DEC = values[12*res + 8]
                    resonance.DWC = values[12*res + 9]
                    resonance.GRC = values[12*res + 10]
                    resonance.GIC = values[12*res + 11]
                    self.resonances.append(resonance)


class RMatrixLimited(ResonanceRange):
    """R-Matrix Limited resolved resonance formalism data. This is identified by
    LRF=7 in the ENDF-6 format.

    Parameters
    ----------
    emin : float
        Minimum energy of the resolved resonance range in eV
    emax : float
        Maximum energy of the resolved resonance range in eV
    nro : int
        Flag designating energy-dependence of scattering radius (NRO). A value
        of 0 indicates it is energy-independent, a value of 1 indicates it is
        energy-dependent.
    naps : int
        Flag controlling use of the channel and scattering radius (NAPS)

    Attributes
    ----------
    IFG : int
        Flag indicating whether channel widths in eV or reduced-width amplitudes
        in eV^1/2 are given
    KRM : int
        Flag to specify which formulae for the R-matrix are to be used
    particle_pairs : list of dict
        List of particle pairs. Each particle pair is represented by a
        dictionary that contains the mass, atomic number, spin, and parity of
        each particle as well as other characteristics.
    spin_groups : list of dict
        List of spin groups. Each spin group is characterized by channels,
        resonance energies, and resonance widths.

    """

    def __init__(self, emin, emax, nro, naps):
        super(RMatrixLimited, self).__init__(emin, emax, nro, naps)
        self.IFG = None
        self.KRM = None
        self.particle_pairs = []
        self.spin_groups = []

    def read(self, ev):
        """Read R-Matrix limited resonance data from an ENDF evaluation.

        Parameters
        ----------
        ev : Evaluation
            ENDF evaluation positioned at the start of a resonance range
            subsection

        """
        items = ev._get_cont_record()
        self.IFG = items[2]  # reduced width amplitude?
        self.KRM = items[3]  # Specify which formulae are used
        n_spin_groups = items[4]  # Number of Jpi values (NJS)
        KRL = items[5]  # Flag for non-relativistic kinematics

        items, values = ev._get_list_record()
        n_pairs = items[5]//2  # Number of particle pairs (NPP)
        for i in range(n_pairs):
            pp = {'mass_a': values[12*i],
                  'mass_b': values[12*i + 1],
                  'z_a': values[12*i + 2],
                  'z_b': values[12*i + 3],
                  'spin_a': values[12*i + 4],
                  'spin_b': values[12*i + 5],
                  'q': values[12*i + 6],
                  'pnt': values[12*i + 7],
                  'shift': values[12*i + 8],
                  'MT': values[12*i + 9],
                  'parity_a': values[12*i + 10],
                  'parity_b': values[12*i + 11]}
            self.particle_pairs.append(pp)

        # loop over spin groups
        for i in range(n_spin_groups):
            sg = {'channels': []}
            self.spin_groups.append(sg)

            items, values = ev._get_list_record()
            J = items[0]
            parity = items[1]
            kbk = items[2]
            kps = items[3]
            n_channels = items[5]
            for j in range(n_channels):
                channel = {}
                channel['particle_pair'] = values[6*j]
                channel['l'] = values[6*j + 1]
                channel['spin'] = values[6*j + 2]
                channel['boundary'] = values[6*j + 3]
                channel['effective_radius'] = values[6*j + 4]
                channel['true_radius'] = values[6*j + 5]
                sg['channels'].append(channel)

            items, values = ev._get_list_record()
            n_resonances = items[3]
            sg['resonance_energies'] = np.zeros(n_resonances)
            sg['resonance_widths'] = np.zeros((n_channels, n_resonances))
            m = n_channels//6 + 1
            for j in range(n_resonances):
                sg['resonance_energies'][j] = values[m*j]
                for k in range(n_channels):
                    sg['resonance_widths'][k,j] = values[m*j + k + 1]

            # Optional extension (Background R-Matrix)
            if kbk > 0:
                items, values = ev._get_list_record()
                lbk = items[4]
                if lbk == 1:
                    params, rbr = ev._get_tab1_record()
                    params, rbi = ev._get_tab1_record()

            # Optional extension (Tabulated phase shifts)
            if kps > 0:
                items, values = ev._get_list_record()
                lps = items[4]
                if lps == 1:
                    params, psr = ev._get_tab1_record()
                    params, psi = ev._get_tab1_record()


class Unresolved(ResonanceRange):
    """Unresolved resonance parameters as identified by LRU=2 in MF=2.

    Parameters
    ----------
    emin : float
        Minimum energy of the unresolved resonance range in eV
    emax : float
        Maximum energy of the unresolved resonance range in eV
    nro : int
        Flag designating energy-dependence of scattering radius (NRO). A value
        of 0 indicates it is energy-independent, a value of 1 indicates it is
        energy-dependent.
    naps : int
        Flag controlling use of the channel and scattering radius (NAPS)

    Attributes
    ----------
    scattering_radius : float
        Scattering radius in units of 10^-12 cm
    LSSF : int
        Flag governing interpretation of file 3 cross sections
    spin : float
        Spin of the target nucleus
    l_values : ndarray of int
        Neutron orbital angular momentum values
    parameters : dict of dict
        Dictionary whose keys are l-values and whose values are dictionaries
        containing unresolved resonance parameters

    """

    def __init__(self, emin, emax, nro, naps):
        super(Unresolved, self).__init__(emin, emax, nro, naps)

    def read(self, ev):
        """Read unresolved resonance data from an ENDF evaluation.

        Parameters
        ----------
        ev : Evaluation
            ENDF evaluation positioned at the start of a resonance range
            subsection

        """
        # Read energy-dependent scattering radius if present
        if self.nro != 0:
            params, self.scattering_radius = ev._get_tab1_record()

        # Get SPI, AP, and LSSF
        if not (self.fission_widths and self.LRF == 1):
            items = ev._get_cont_record()
            self.spin = items[0]
            if self.nro == 0:
                self.scattering_radius = items[1]
            self.LSSF = items[2]

        if not self.fission_widths and self.LRF == 1:
            # Case A -- fission widths not given, all parameters are
            # energy-independent
            NLS = items[4]
            self.l_values = np.zeros(NLS)
            self.parameters = {}
            for ls in range(NLS):
                items, values = ev._get_list_record()
                l = items[2]
                NJS = items[5]
                self.l_values[ls] = l
                params = {}
                self.parameters[l] = params
                params['d'] = np.asarray(values[0::6])
                params['j'] = np.asarray(values[1::6])
                params['amun'] = np.asarray(values[2::6])
                params['gn0'] = np.asarray(values[3::6])
                params['gg'] = np.asarray(values[4::6])
                # params['gf'] = np.zeros(NJS)

        elif self.fission_widths and self.LRF == 1:
            # Case B -- fission widths given, only fission widths are
            # energy-dependent
            items, self.energies = ev._get_list_record()
            self.spin = items[0]
            if self.nro == 0:
                self.scatter_radius = items[1]
            self.LSSF = items[2]
            NE, NLS = items[4:6]
            self.l_values = np.zeros(NLS, int)
            self.parameters = {}
            for ls in range(NLS):
                items = ev._get_cont_record()
                l = items[2]
                NJS = items[4]
                self.l_values[ls] = l
                params = {}
                self.parameters[l] = params
                params['d'] = np.zeros(NJS)
                params['j'] = np.zeros(NJS)
                params['amun'] = np.zeros(NJS)
                params['gn0'] = np.zeros(NJS)
                params['gg'] = np.zeros(NJS)
                params['gf'] = []
                for j in range(NJS):
                    items, values = ev._get_list_record()
                    muf = items[3]
                    params['d'][j] = values[0]
                    params['j'][j] = values[1]
                    params['amun'][j] = values[2]
                    params['gn0'][j] = values[3]
                    params['gg'][j] = values[4]
                    params['gf'].append(np.asarray(values[6:]))

        elif self.LRF == 2:
            # Case C -- all parameters are energy-dependent
            NLS = items[4]
            self.l_values = np.zeros(NLS)
            self.parameters = {}
            for ls in range(NLS):
                items = ev._get_cont_record()
                l = items[2]
                NJS = items[4]
                self.l_values[ls] = l
                params = {}
                self.parameters[l] = params
                params['j'] = np.zeros(NJS)
                params['amux'] = np.zeros(NJS)
                params['amun'] = np.zeros(NJS)
                params['amug'] = np.zeros(NJS)
                params['amuf'] = np.zeros(NJS)
                params['energies'] = []
                params['d'] = []
                params['gx'] = []
                params['gn0'] = []
                params['gg'] = []
                params['gf'] = []
                for j in range(NJS):
                    items, values = ev._get_list_record()
                    ne = items[5]
                    params['j'][j] = items[0]
                    params['amux'][j] = values[2]
                    params['amun'][j] = values[3]
                    params['amug'][j] = values[4]
                    params['amuf'][j] = values[5]
                    params['energies'].append(np.asarray(values[6::6]))
                    params['d'].append(np.asarray(values[7::6]))
                    params['gx'].append(np.asarray(values[8::6]))
                    params['gn0'].append(np.asarray(values[9::6]))
                    params['gg'].append(np.asarray(values[10::6]))
                    params['gf'].append(np.asarray(values[11::6]))


class ScatteringRadius(ResonanceRange):
    """Energy range with no resonances and only a scattering radius.

    Parameters
    ----------
    emin : float
        Minimum energy of the resolved resonance range in eV
    emax : float
        Maximum energy of the resolved resonance range in eV
    nro : int
        Flag designating energy-dependence of scattering radius (NRO). A value
        of 0 indicates it is energy-independent, a value of 1 indicates it is
        energy-dependent.
    naps : int
        Flag controlling use of the channel and scattering radius (NAPS)

    """
    def __init__(self, emin, emax, nro, naps):
        super(ScatteringRadius, self).__init__(emin, emax, nro, naps)


class AdlerResonance(object):
    pass


class Resonance(object):
    pass


_formalisms = {1: SingleLevelBreitWigner, 2: MultiLevelBreitWigner,
               3: ReichMoore, 4: AdlerAdler, 7: RMatrixLimited}
