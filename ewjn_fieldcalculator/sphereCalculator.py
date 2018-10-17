import os
import subprocess
import time
import logging

import ewjn_fieldcalculator.fieldCalculator as fieldCalculator

logger = logging.getLogger(__name__)

class CalculationGroup(object):
	def __init__(self, directory_base, sphere_radius):
		self.sphere_radius = sphere_radius
		timestamp = time.strftime("%Y%m%d-%H%M%S")
		self.directory_name = "{}-{}".format(directory_base, timestamp)
		logger.info("Will output to directory {}".format(self.directory_name))
		os.makedirs(self.directory_name, exist_ok=True)
	def calculateJ(self, dx):
		j_filename = "sphereJ--r{}--dx{}.csv".format(str(self.sphere_radius), str(dx))
		return SphereJCalculation(self, j_filename, dx)


class SphereJCalculation(object):
	def __init__(self, calculation_group, current_filename, dx):
		self.calculation_group = calculation_group
		self.sphere_radius = calculation_group.sphere_radius
		self.current_filename = current_filename
		self.dx = dx
		
		logger.debug("##################### MATHEMATICA OUTPUT ###############################")
		subprocess.run(["wolframscript", "-f", "wl/sphereCurrent.wl", os.path.abspath(self.calculation_group.directory_name), self.current_filename, str(self.sphere_radius), str(self.dx)])
		logger.debug("##################### END MATHEMATICA OUTPUT ###########################")
	def calculateBAlongX(self):
		return fieldCalculator.BCalculationAlongXAxis(self)