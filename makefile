f = calQ.f90

FC      = ifort
FCFLAGS = -shared -fpic -O2 -xHost -openmp 
LDFLAG = -liomp5 

symmetry_function.so: ${f}
	${FC} ${FCFLAGS} ${LDFLAG} -o calQ.so ${f}

.PHONY: clean

clean:
	-rm -rf *.o *.so *.mod 

