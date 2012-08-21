all: testfileIO testimageIO testMatrix project1

INCLUDE = -I../include
LIB = -L../lib

testfileIO: testfileIO.o
	g++ -o testfileIO testfileIO.o $(LIB) -lmatrix

testfileIO.o: testfileIO.cpp
	g++ -c testfileIO.cpp $(INCLUDE)

testimageIO: testimageIO.o
	g++ -o testimageIO testimageIO.o $(LIB) -lmatrix

testimageIO.o: testimageIO.cpp
	g++ -c testimageIO.cpp $(INCLUDE)

testMatrix: testMatrix.o
	g++ -o testMatrix testMatrix.o $(LIB) -lmatrix

testMatrix.o: testMatrix.cpp
	g++ -c testMatrix.cpp $(INCLUDE)

project1.o: project1.cpp
	g++ -c project1.cpp $(INCLUDE)

project1: project1.o
	g++ -o project1 project1.o $(LIB) -lmatrix




