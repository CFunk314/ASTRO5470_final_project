FC=gfortran
FLAGS=-c -w

OBJECTS = defs.o grids.o integrate.o solve.o dump.o main.o

a.out: $(OBJECTS)
	$(FC) $(OBJECTS) -o a.out

defs.o: defs.f90
	$(FC) $(FLAGS) defs.f90

grids.o: grids.f90
	$(FC) $(FLAGS) grids.f90

integrate.o: integrate.f90
	$(FC) $(FLAGS) integrate.f90

solve.o: solve.f90
	$(FC) $(FLAGS) solve.f90

dump.o: dump.f90
	$(FC) $(FLAGS) dump.f90

main.o: main.f90
	$(FC) $(FLAGS) main.f90

clean:
	rm *.o
	rm *.mod
	rm a.out
