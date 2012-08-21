OBJ = Matrix.o fileIO.o matrixProcessing.o eigen.o

AR = ar
INCLUDE = -I../include
all:
	${MAKE} libmatrix.a

libmatrix.a: $(OBJ)
	$(AR) rvu $@ $(OBJ)
	ranlib $@

eigen.o: eigen.cpp
	g++ -c eigen.cpp $(INCLUDE)

matrixProcessing.o: matrixProcessing.cpp
	g++ -c matrixProcessing.cpp $(INCLUDE)

fileIO.o: fileIO.cpp
	g++ -c fileIO.cpp $(INCLUDE)

Matrix.o: Matrix.cpp
	g++ -c Matrix.cpp $(INCLUDE)

clean:
	-rm *.o *~ 	