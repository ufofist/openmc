#!/usr/bin/env python

import os
import sys
import glob
import hashlib
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
from input_set import InputSet
import openmc

class TestPyAPIPlot(PyAPITestHarness):
    def _build_inputs(self):
        self._input_set.build_default_materials_and_geometry()

        root = self._input_set.geometry.root_universe
        root.plot(width=(50, 50), filename='plot.png', seed=1)

    def _test_output_created(self):
        """Make sure the plot has been created."""
        assert os.path.isfile(os.path.join(os.getcwd(), 'plot.png')), \
            'Python API-generated plot was not found.'

    def execute_test(self):
        try:
            self._build_inputs()
            self._test_output_created()
            results = self._get_results()
            self._write_results(results)
            self._compare_results()
        finally:
            self._cleanup()

    def update_results(self):
        """Update the results_true using the current version of OpenMC."""
        try:
            self._build_inputs()
            self._test_output_created()
            results = self._get_results()
            self._write_results(results)
            self._overwrite_results()
        finally:
            self._cleanup()

    def _get_results(self):
        data = open('plot.png', 'rb').read()
        sha512 = hashlib.sha512()
        sha512.update(data)
        return sha512.hexdigest()

    def _cleanup(self):
        """Delete XMLs, statepoints, tally, and test files."""
        super(PyAPITestHarness, self)._cleanup()
        output = [os.path.join(os.getcwd(), 'plot.png')]
        for f in output:
            if os.path.exists(f):
                os.remove(f)


if __name__ == '__main__':
    test = TestPyAPIPlot('')
    test.main()
