PS_PATH = -I ./src/ -I ./my-src -I ./

CC = CC
CFLAGS = -O3 -std=c++17
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard ./my-src/*.cpp))
CPPOBJM = $(patsubst %.cpp, %.o, $(wildcard ./*.cpp))
CPPHDRS = $(wildcard *.h) $(wildcard ./my-src/*.h)
PROGRAM = sph.out

.PHONY:	clean all

all:	$(CPPOBJS) $(CPPOBJM) $(CPPHDRS)
	@echo "Linking object files..."
	@$(CC) $(CFLAGS) $(WARNINGS) $(CPPOBJM) $(CPPOBJS) -o $(PROGRAM) $(LIBS) $(PS_PATH)
	@echo "Link Success! [$(PROGRAM)]"

%.o:	%.cpp $(CPPHDRS) Makefile
	@echo "Building $< ..."	
	@$(CC) -c $< -o $(patsubst %.cpp, %.o, $<) $(CFLAGS) $(WARNINGS) $(PS_PATH)
	@echo "[$< OK]"

clean:
	-rm *.out *.o ./my-src/*.o

makeinitial:
	@$(CC) make-sphere-from-relaxed-condition-for-ANEOS.c -o make-initial.out $(PS_PATH)
	@echo "make make-initial!"

makeimpactcondition:
	@$(CC) make-impact-condition-from-two-spheres.c -o make-impact-condition-from-two-spheres.out $(PS_PATH)
	@echo "make make-impact-condition-from-two-spheres!"

convertcgstomks:
	@$(CC) convert-from-cgs-to-mks.c -o convert-from-cgs-to-mks.out $(PS_PATH)
	@echo "make convert-from-cgs-to-mks!"


