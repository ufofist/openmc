import warnings

from six import string_types

import openmc.checkvalue as cv
from openmc.checkvalue import check_iterable_type


class Nuclide(object):
    """A nuclide that can be used in a material.

    Parameters
    ----------
    name : str
        Name of the nuclide, e.g. 'U235'

    Attributes
    ----------
    name : str
        Name of the nuclide, e.g. 'U235'
    scattering : 'data' or 'iso-in-lab' or None
        The type of angular scattering distribution to use

    """

    def __init__(self, name=''):
        # Initialize class attributes
        self._name = ''
        self._scattering = None

        self._poly_coeffs = None
        self._poly_type = None

        # Set the Material class attributes
        self.name = name

    def __eq__(self, other):
        if isinstance(other, Nuclide):
            if self.name != other.name:
                return False
            else:
                return True
        elif isinstance(other, string_types) and other == self.name:
            return True
        else:
            return False

    def __ne__(self, other):
        return not self == other

    def __gt__(self, other):
        return repr(self) > repr(other)

    def __lt__(self, other):
        return not self > other

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        string = 'Nuclide    -    {0}\n'.format(self._name)
        if self.scattering is not None:
            string += '{0: <16}{1}{2}\n'.format('\tscattering', '=\t',
                                                self.scattering)
        return string

    @property
    def name(self):
        return self._name

    @property
    def scattering(self):
        return self._scattering

    @property
    def poly_coeffs(self):
        return self._poly_coeffs

    @property
    def poly_type(self):
        return self._poly_type

    @name.setter
    def name(self, name):
        cv.check_type('name', name, string_types)
        self._name = name

        if '-' in name:
            self._name = name.replace('-', '')
            self._name = self._name.replace('Nat', '0')
            if self._name.endswith('m'):
                self._name = self._name[:-1] + '_m1'

            msg = 'OpenMC nuclides follow the GND naming convention. Nuclide ' \
                  '"{}" is being renamed as "{}".'.format(name, self._name)
            warnings.warn(msg)

    @scattering.setter
    def scattering(self, scattering):
        if not scattering in ['data', 'iso-in-lab', None]:
            msg = 'Unable to set scattering for Nuclide to {0} which ' \
                  'is not "data", "iso-in-lab", or None'.format(scattering)
            raise ValueError(msg)

        self._scattering = scattering

    @poly_coeffs.setter
    def poly_coeffs(self, poly_coeffs):
        check_iterable_type('poly_coeffs)', poly_coeffs, float)
        self._poly_coeffs = poly_coeffs

    @poly_type.setter
    def poly_type(self, poly_type):
        check_type('polynomial type', poly_type, string_types)
        self._poly_type = poly_type
