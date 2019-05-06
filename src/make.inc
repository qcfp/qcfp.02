# UNIFIED MAKEFILE FOR QCFP
# NOT USED IN contructor-excitons


SYSTEM=GNULINUX
#SYSTEM=WINDOWS
MPI_COMPILE=FALSE


# set the compilers
ifeq ($(MPI_COMPILE),TRUE)
CCMPI=/Users/dariusa/usr/bin/mpicxx
CCMPIOPT=-O2 -DMPIPROCESSING
#CCMPIOPT=--debug
endif


# FOR MACOSX AND GNULINUX:
ifeq ($(SYSTEM),GNULINUX)
COPY=cp
DEL=rm -f
#CC=llvm-g++
#CC=g++-5
CC=g++
CCOPT=-O2 -c  -I ~/usr/include
#CCOPT=--debug -c  -I ~/usr/include
CLOPT= -lfftw3 -llapack -lblas  -L ~/usr/lib
LIBO=-o
BINO=-o
LIBS=../../lib/
BINS=../../bin/
endif

# FOR WINDOWS 10 MS Visual Studio:
ifeq ($(SYSTEM),WINDOWS)
COPY=copy
DEL=del /Q
CC=cl
#CCOPT=/W4 /EHsc /c /O2  /I ..\..\..\fftw-3.3.5-dll32
CCOPT=/W4 /EHsc /c  /I ..\..\..\fftw-3.3.5-dll32
CLOPT=/link \
	/LIBPATH ..\..\..\fftw-3.3.5-dll32\libfftw3-3.lib \
	/LIBPATH ..\..\..\lapack\liblapack.lib \
	/LIBPATH ..\..\..\lapack\libblas.lib
LIBO=/Fo
BINO=/Fe
LIBS=..\..\lib\\
BINS=..\..\bin\\
endif


# general rules

ifeq ($(BINARY),TRUE)
default: depends $(BINS)z.$(PROJECT) binaries
else
default: depends $(BINS)z.$(PROJECT)
endif
	echo "Done"

$(BINS)z.$(PROJECT): \
 $(LIBS)$(PROJECT).obj \
 $(LIBS)$(PROJECT).main.obj \
 $(foreach val,$(DEPENDENCIES), $(LIBS)$(val).obj )
	$(CC) $(BINO)$(BINS)z.$(PROJECT) \
	$(LIBS)$(PROJECT).obj $(LIBS)$(PROJECT).main.obj  \
	$(foreach val,$(DEPENDENCIES), $(LIBS)$(val).obj) $(CLOPT)

$(LIBS)$(PROJECT).obj: $(PROJECT).cpp $(PROJECT).hpp
	$(CC) $(PROJECT).cpp $(LIBO)$(LIBS)$(PROJECT).obj $(CCOPT)

$(LIBS)$(PROJECT).main.obj: main.cpp $(PROJECT).hpp
	$(CC) main.cpp  $(LIBO)$(LIBS)$(PROJECT).main.obj $(CCOPT)

$(foreach val,$(DEPENDENCIES), $(LIBS)$(val).obj ):

depends:
	$(foreach val,$(DEPENDENCIES), $(MAKE) $(LIBS)$(val).obj -C ../$(val) && )  echo "Done"

clean:
	$(DEL) $(LIBS)$(PROJECT)* $(BINS)*$(PROJECT)*

binaries: $(BINS)z.$(PROJECT)
	$(COPY) $(BINS)z.$(PROJECT) $(BINS)qcfp.$(PROJECT).exe
