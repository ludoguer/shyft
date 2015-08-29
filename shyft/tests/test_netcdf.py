﻿"""
Tests for the netcdf datasets.
"""

from __future__ import print_function
from __future__ import absolute_import

import os

import unittest
import numpy as np

from shyft import api
from shyft.orchestration2 import config_constructor, cell_extractor, CalibrationConfig
from shyft.orchestration2.shyft_runner import Simulator, Calibrator


class Simulation():
    def setUp(self):
        # Get the configuration section
        config_file, section = self.config_file, self.section
        config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "netcdf/%s" % config_file)
        config = config_constructor(config_file, section)

        # Build the simulator
        self.simulator = simulator = Simulator(config)
        time_axis = api.Timeaxis(config.start_time, config.run_time_step, config.number_of_steps)
        simulator.build_model(time_axis.start(), time_axis.delta(), time_axis.size())
        simulator.run_model()

    def test_simulation(self):
        cells = self.simulator.model.get_cells()
        assert len(cells) == 4
        expected_results = {
            "total_discharge": [0.00844, 0.00586, 0.00836, 0.00824],
            "discharge": [0.0431, 0.0444, 0.0429, 0.0423],
            "snow_storage": [41.8, 49.7, 41.5, 41.2],
            "temperature": [9.14, 8.41, 9.27, 9.62],
            "precipitation": [0.173, 0.178, 0.172, 0.170],
        }

        # for the fun of it, demonstrate how to use cell_statistics
        cids = api.IntVector()
        temperature = self.simulator.model.statistics.temperature(cids)
        precipitation = self.simulator.model.statistics.precipitation(cids)
        discharge = self.simulator.model.statistics.discharge(cids)
        assert discharge.size() > 0
        assert temperature.size() > 0
        assert precipitation.size() > 0
        assert discharge.size() > 0

        for i, cell in enumerate(cells):
            for param in expected_results:
                value = cell_extractor[param](cell)
                if type(value) is np.ndarray:
                    assert len(np.where(value != value)[0]) == 0
                    # Take the mean value as a gross estimator
                    value = value.sum() / len(value)
                #print("param, cell, value:", param, i, value)
                np.testing.assert_allclose(np.float64(value), np.float64(expected_results[param][i]), rtol=1e-2)


# Some examples of simulation.  Feel free to add more.
class Simulation1(Simulation, unittest.TestCase):
    config_file = "configuration.yaml"
    section = "Himalayas"


class Calibration():
    def setUp(self):
        # Get the configuration section
        config_file, section = self.config_file, self.section
        config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "netcdf/%s" % config_file)
        config = CalibrationConfig(config_file, section)

        # Build the calibrator
        time_axis = api.Timeaxis(config.model_config.start_time, config.model_config.run_time_step,
                                 config.model_config.number_of_steps)
        self.calibrator = Calibrator(config)
        self.calibrator.init(time_axis)

    def test_calibration(self):
        calibr_results = self.calibrator.calibrate(tol=1.0e-5)
        # print("calibrated results:", calibr_results)
        expected_results = {
            'wind_const': 1.0, 'max_albedo': 0.9, 'p_corr_scale_factor': 1.0, 'fast_albedo_decay_rate': 10.0,
            'TX': -0.5, 'glacier_albedo': 0.4, 'surface_magnitude': 30.0, 'snowfall_reset_depth': 5.0,
            'wind_scale': 3.5, 'slow_albedo_decay_rate': 30.0, 'ae_scale_factor': 1.5, 'c3': -0.1,
            'c2': 1.0, 'c1': -0.5, 'snow_cv': 0.4,
            'min_albedo': 0.6, 'max_water': 0.1}
        # Check that the results are as expected
        for key, val in calibr_results.iteritems():
            assert val == expected_results[key]
            np.testing.assert_allclose(np.float64(val), np.float64(expected_results[key]), rtol=1e-2)


# Some examples of calibration.  Feel free to add more.
class Calibration1(Calibration, unittest.TestCase):
    config_file = "calibration.yaml"
    section = "Himalayas"