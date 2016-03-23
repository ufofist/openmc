#!/usr/bin/env python

from functools import reduce
from operator import mul
import os
import sys

import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

from openmc import *
from openmc.source import Source
from openmc.stats import Box
from openmc.mgxs_library import XSdata, MGXSLibraryFile
from openmc.mgxs import EnergyGroups


# Create 1g cross sections
onegroup = EnergyGroups([1.0e-5, 20.])
toy_xs = XSdata('toy.71c', onegroup)
toy_xs.total = np.array([0.3])
toy_xs.kT = 2.53e-8
toy_xs.alias = 'toy.71c'
toy_xs.order = 0
toy_xs.fissionable = True
toy_xs.absorption = np.array([0.03])
toy_xs.nu_fission = np.array([0.03])
toy_xs.chi = np.array([1.0])
toy_xs.fission = np.array([0.012])
toy_xs.scatter = np.array([[[0.27]]])
toy_xs.multiplicity = np.array([[1.0]])

library = MGXSLibraryFile(onegroup)
library.add_xsdata(toy_xs)
library.export_to_xml('mgxs.xml')

# Define material
toy = Material(material_id=1, name='toy')
toy.set_density('macro', 1.0)
toy.add_macroscopic(Macroscopic('toy', '71c'))
materials_file = MaterialsFile()
materials_file.add_material(toy)
materials_file.export_to_xml()

# Define surfaces
left = XPlane(surface_id=1, x0=-200., boundary_type='reflective')
right = XPlane(surface_id=2, x0=200., boundary_type='reflective')
back = YPlane(surface_id=3, y0=-200., boundary_type='reflective')
front = YPlane(surface_id=4, y0=200., boundary_type='reflective')
bottom = ZPlane(surface_id=5, z0=-200., boundary_type='reflective')
top = ZPlane(surface_id=6, z0=200., boundary_type='reflective')

# Define box
box = Cell(cell_id=1)
box.fill = toy
box.region = +left & -right & +back & -front & +bottom & -top

# Create root universe, assign to geometry
root = Universe(universe_id=0)
root.add_cell(box)
geometry = Geometry()
geometry.root_universe = root
geometry_file = GeometryFile()
geometry_file.geometry = geometry
geometry_file.export_to_xml()

# Settings
settings = SettingsFile()
settings.particles = 100000
settings.batches = 10
settings.inactive = 0
settings.source = Source(space=Box([-200., -200., -200.], [200., 200., 200.]))
settings.energy_mode = 'multi-group'
settings.cross_sections = 'mgxs.xml'

# Tallies
mesh = Mesh(mesh_id=1)
mesh.lower_left = (-200., -200., -200.)
mesh.upper_right = (200., 200., 200.)
mesh.dimension = (4, 4, 4)
tally = Tally(tally_id=1)
mesh_filter = Filter(type='mesh')
mesh_filter.mesh = mesh
tally.filters = [mesh_filter]
tally.estimator = 'collision'
tally.scores = ['fission']
tallies_file = TalliesFile()
tallies_file.add_mesh(mesh)
tallies_file.add_tally(tally)
tallies_file.export_to_xml()

exe = Executor()

n = reduce(mul, mesh.dimension)
s = np.zeros(n)
s2 = np.zeros(n)

analytical_answer = 0.012/0.03/n

n_realization = 10
for i in range(n_realization):
    settings.seed = i + 1
    settings.export_to_xml()

    exe.run_simulation()
    os.rename('statepoint.{}.h5'.format(settings.batches),
              'realization_{}.h5'.format(i + 1))

    # Get fission scores
    sp = StatePoint('realization_{}.h5'.format(i + 1))
    df = sp.tallies[1].get_pandas_dataframe()
    fission = df[df['score'] == 'fission']

    # Accumulate
    s += fission['mean']
    s2 += fission['mean']**2

mean = s/n_realization
std_dev = np.sqrt((s2/n_realization - mean*mean)/(n_realization - 1))

alpha = 0.05
t_value = scipy.stats.t.ppf(1 - alpha/2, n_realization - 1)

plt.errorbar(range(n), mean, yerr=t_value*std_dev, fmt='bo')
plt.plot((0, n), (analytical_answer, analytical_answer), 'k--')
plt.show()
