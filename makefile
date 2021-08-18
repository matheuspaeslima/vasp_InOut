FC=ifort

.SUFFIXES: .f .f90 .o

EXE_NAME = main.x
FLAGS =

default: executable

OBJS = m_lu.o m_contcar.o uni64.o main.o 

executable: $(OBJS)
	$(FC) -o $(EXE_NAME) $(OBJS) $(FLAGS)

.f.o:
	$(FC) -c $(FLAGS) $<
.f90.o:
	$(FC) -c $(FLAGS) $<

clean:
	rm -rf *.o $(EXE_NAME)                                                                                                                                                                               
