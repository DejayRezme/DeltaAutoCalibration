deltaCalibrateLBFGS: deltaCalibrateLBFGS.o lbfgs.o
	g++ -o deltaCalibrateLBFGS deltaCalibrateLBFGS.o lbfgs.o

deltaCalibrateLBFGS.o: src/deltaCalibrateLBFGS.c src/lbfgs.h
	g++ -c src/deltaCalibrateLBFGS.c

lbfgs.o: src/lbfgs.c src/lbfgs.h
	g++ -c src/lbfgs.c