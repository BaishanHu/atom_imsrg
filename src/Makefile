#CC = mpiCC
CC = g++
#CC = icc

# Compiler flags:
#  -DNO_ODE=1 compiles without boost/ode package used for flow equation solver
#  -DOPENBLAS_NOUSEOMP=1 removes parallel blocks which take threads away from OPENBLAS
#                        to be used if OpenBlas was compiled without the USE_OMP flat

INCLUDE   = -I./armadillo
FLAGS     = -O3 -march=native -std=c++11 -fopenmp -fPIC #-flto
SOFLAGS   = -O3 -march=native -std=c++11 -fopenmp -fPIC #-flto
WARNFLAGS = -Wall -Wno-comment -Wno-deprecated-declarations
FLAGS    += $(WARNFLAGS)
SOFLAGS  += $(WARNFLAGS)
INSTDIR = $(HOME)

ifeq ($(DEBUG),on)
 FLAGS     = -march=native -std=c++11 -fopenmp -fPIC  
 SOFLAGS   = -march=native -std=c++11 -fopenmp -fPIC  
 FLAGS += -g
 SOFLAGS += -g
endif

PYTHONFLAGS = -I/usr/include/python2.7 -lboost_python #-lpython 

THEHOST = $(shell if [ `hostname|grep jrl` ]; then echo jureca; elif [ `hostname|grep cougar` ]; then echo cougar; elif [ `hostname|grep cronos` ]; then echo cronos; elif [ `hostname|grep oak` ]; then echo oak; elif [ `hostname|grep cedar` ]; then echo cedar; else echo other; fi)

ALL =  libIMSRG.so imsrg++
LIBS = -lopenblas -lgsl -lm -lz -lhdf5_cpp -lhdf5 -fext-numeric-literals -I -json -fPIC

ifeq ($(HDF5),on)
 LIBS += -lhdf5_cpp -lhdf5
else
 FLAGS += -DNO_HDF5   # By default, don't bother building with HDF5, since it's not used very often and can make building a pain
endif

ifeq ($(THEHOST),jureca) # specific options for jureca cluster
# FLAGS += -DOPENBLAS_NOUSEOMP=1
 SOFLAGS += -fuse-ld=bfd
endif
ifeq ($(THEHOST),cougar) # specific options for cougar cluster
 LIBS += -llapack
 PYTHONFLAGS = -I/opt/anaconda/include/python2.7 -L/opt/boost/1.58.0/lib -lboost_python #-lpython 
endif
ifeq ($(THEHOST),cronos)
 FLAGS += -DOLD_BOOST=1
 SOFLAGS += -DOLD_BOOST=1
 ifneq ($(PYTHON),off)
  ALL += pyIMSRG.so
  endif
endif

ifeq ($(THEHOST),oak)
#  FLAGS += -DOLD_BOOST
#  FLAGS += -DOPENBLAS_NOUSEOMP=1
  NEWLIBS := $(filter-out -lopenblas,$(LIBS))
  LIBS = $(NEWLIBS) -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
 ifneq ($(PYTHON),off)
  ALL += pyIMSRG.so
  endif
endif

ifeq ($(THEHOST),cedar)
 LIBS += -lmkl -lmkl_gnu_thread -lmkl_core -lpthread -lgomp -liomp5 
 ifneq ($(PYTHON),off)
  ALL += pyIMSRG.so
 endif
endif

ifeq ($(THEHOST),other)  # default options. assumes boost and python are set up nicely.
 LIBS += -llapack
 ifneq ($(PYTHON),off)
  ALL += pyIMSRG.so
 endif
endif

ifeq ($(HDF5),on)
# The new version of Ubuntu puts the hdf5 libraries in a different location
 NEWHDF5LOC = $(shell if [ -d /usr/include/hdf5/serial ]; then echo yes; else echo no; fi)
 ifeq ($(NEWHDF5LOC),yes)
   INCLUDE += -I/usr/include/hdf5/serial
   NEWLIBS := $(filter-out -lhdf5,$(LIBS))
   LIBS = $(NEWLIBS) -lhdf5_serial
 endif
endif

all: $(ALL)

install:
	@if [ ! -d $(INSTDIR)/lib ] ; then \
	  mkdir $(INSTDIR)/lib; \
	fi
	@if [ ! -d $(INSTDIR)/include ] ; then \
	  mkdir $(INSTDIR)/include; \
	fi 
	ln -sf $(PWD)/libIMSRG.so $(INSTDIR)/lib/libIMSRG.so
	ln -sfT $(PWD)/armadillo $(INSTDIR)/include/armadillo
	@if [ -f pyIMSRG.so ] ; then\
	  ln -sf $(PWD)/pyIMSRG.so $(INSTDIR)/lib/pyIMSRG.so;\
	fi
	@for x in *.hh; do \
	 echo linking $(PWD)/$$x  '=>'  $(INSTDIR)/include/$$x;\
	 ln -sf $(PWD)/$$x $(INSTDIR)/include/$$x; \
	done
	@if [ -d $(INSTDIR)/bin ] ; then \
	  echo linking $(PWD)/imsrg++  '=>' $(INSTRDIR)/bin/imsrg++;\
	  ln -sf $(PWD)/imsrg++ $(INSTDIR)/bin/imsrg++;\
	fi
	@printf "\n\nDone installing.\n\n"
	@echo '*********************************************************************'
	@echo '* Make sure libIMSRG.so is in your LIBRARY_PATH and LD_LIBRARY_PATH *'
	@echo '*********************************************************************'
	 
OBJ = ModelSpace.o TwoBodyME.o ThreeBodyME.o Operator.o  ReadWrite.o\
      HartreeFock.o imsrg_util.o Generator.o IMSRGSolver.o AngMom.o\
      boost_src/gzip.o boost_src/zlib.o  IMSRGProfiler.o

mysrg: main.cc $(OBJ)
	$(CC) $^ -o $@ $(INCLUDE) $(LIBS) $(FLAGS) 

boost_src/gzip.o: boost_src/gzip.cpp
	$(CC) -c $^ -o $@ $(INCLUDE) $(LIBS) $(FLAGS) 

boost_src/zlib.o: boost_src/zlib.cpp
	$(CC) -c $^ -o $@ $(INCLUDE) $(LIBS) $(FLAGS) 


%.o: %.cc %.hh
	$(CC) -c $*.cc -o $@ $(INCLUDE) $(LIBS) $(FLAGS) 



libIMSRG.so: $(OBJ)
	$(CC) $^ -shared -o $@ $(SOFLAGS)

python: pyIMSRG.so

pyIMSRG.so: $(OBJ)  pyIMSRG.o
	$(CC) $^ -o $@ $(FLAGS) $(INCLUDE) $(LIBS) -shared $(PYTHONFLAGS)

pyIMSRG.o: pyIMSRG.cc
	$(CC) -c $^ -o $@ $(FLAGS) $(INCLUDE) $(LIBS) $(PYTHONFLAGS)

imsrg++: imsrg++.cc libIMSRG.so Parameters.hh
	$(CC) $(INCLUDE) $< -o $@ $(FLAGS) -L$(PWD) -lIMSRG $(LIBS)

clean:
	rm -f *.o *.so boost_src/*.o
