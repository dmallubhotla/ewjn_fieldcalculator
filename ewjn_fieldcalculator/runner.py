import argparse
import logging
import ewjn_fieldcalculator.sphereCalculator as sphereCalculator

def get_parser():
	parser = argparse.ArgumentParser(description="Outputs a correction field for sphere as a file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-r", "--radius", type=int, default=3, help="The radius of the sphere or cylinder")
	parser.add_argument("--working-folder", type=str, default="data", help="The prefix to use for the data folder for current and field values")
	parser.add_argument("output_filename", type=str, help="The filename of the output field")
	return parser

def setUpLogging():
	formatter = logging.Formatter("%(asctime)s | %(levelname)-7s | %(name)s | %(message)s")
	console = logging.StreamHandler()
	console.setFormatter(formatter)
	#Getting root logger cause we want all the logs.
	rootLogger = logging.getLogger()
	rootLogger.setLevel(logging.DEBUG)
	rootLogger.addHandler(console)
	

def main():
	setUpLogging()
	LOGGER = logging.getLogger(__name__)
	LOGGER.info("Running runner")
	
	parser = get_parser()
	args = parser.parse_args()
	LOGGER.debug(args)
	calculator = sphereCalculator.CalculationGroup(args.working_folder, args.radius)
	calculator.calculateJ(1).calculateBAlongX().calculate()
	calculator.calculateJ(0.5).calculateBAlongX().calculate()
	calculator.calculateJ(0.25).calculateBAlongX().calculate()
	calculator.calculateJ(0.1).calculateBAlongX().calculate()
	calculator.calculateJ(0.05).calculateBAlongX().calculate()

	
if __name__ == "__main__":
	main()