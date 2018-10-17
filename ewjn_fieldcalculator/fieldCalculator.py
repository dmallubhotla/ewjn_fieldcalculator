import logging
import numpy
import os

logger = logging.getLogger(__name__)

ZERO_THRES = .0000005
class BCalculationAlongXAxis(object):
	def __init__(self, j_calculation, maxX=10):
		self.j_calculation = j_calculation
		self.calculation_group = j_calculation.calculation_group
		self.dV = self.j_calculation.dx ** 3
		self.maxX = maxX
		j_filename = os.path.join(self.calculation_group.directory_name, self.j_calculation.current_filename)
		logger.info("Will read from {}, scaling with dV = {}".format(j_filename, str(self.dV)))
		self.j_data = numpy.loadtxt(j_filename, delimiter=",")
	def get_integrand_coord(self, val, coord, field_point):
		xp, yp, zp = val[0:3]
		jx, jy, jz = val[3:6]
		x, y, z = field_point
		numeratorz = numpy.cross((jx, jy, jz), (x - xp, y - yp, z - zp))[coord]
		denom = (numpy.linalg.norm([x - xp, y - yp, z]))**3
		return self.dV * numeratorz / denom
	def get_integral_at_point(self, field_point):
		ixs = []
		iys = []
		izs = []
		for point in self.j_data:
			ixs.append(self.get_integrand_coord(point, 0, field_point))
			iys.append(self.get_integrand_coord(point, 1, field_point))
			izs.append(self.get_integrand_coord(point, 2, field_point))
		return sum(ixs), sum(iys), sum(izs)
	def calculate(self):
		logger.info("Calculating B...")
		r = numpy.arange(0, self.maxX + 1, 0.5)
		points = [(x, 0, 0) for x in r]
		fieldpts = [pt + self.get_integral_at_point(pt) for pt in points]
		b_filename = "sphereB--r{}--dx{}--range{}.csv".format(str(self.j_calculation.sphere_radius), str(self.j_calculation.dx), str(self.maxX))
		numpy.savetxt(os.path.join(self.calculation_group.directory_name, b_filename), fieldpts, delimiter=",")
		self.estimate_sphere_volume()
	def estimate_sphere_volume(self):
		vol = 0
		for point in self.j_data:
			jx, jy, jz = point[3:6]
			if (abs(jx) > ZERO_THRES) or (abs(jy) > ZERO_THRES) or (abs(jz) > ZERO_THRES):
				vol += 1
		logger.info("Got a volume of {}".format(str(self.dV * vol)))