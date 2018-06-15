# UNIFIED MAKEFILE FOR QCFP
# NOT USED IN contructor-excitons


#SYSTEM=GNULINUX
SYSTEM=WINDOWS
MPI_COMPILE=FALSE


# set the compilers
ifeq ($(MPI_COMPILE),TRUE)

CCMPI=/Users/dariusa/usr/bin/mpicxx
CCMPIOPT=-O2 -DMPIPROCESSING
#CCMPIOPT=--debug

endif


# FOR MACOSX AND GNULINUX
ifeq ($(SYSTEM),GNULINUX)
CC=llvm-g++
#CC=g++-5
#CC=g++
#CCOPT=-O2 -c
CCOPT=--debug -c
CLOPT= -lfftw3 -llapack -lblas 
INCPATH= -I ~/usr/include 
LIBPATH= -L ~/usr/lib 
LIBS=-o ../aaa-lib/
BINS=-o ../aaa-bin/
endif

# FOR WINDOWS 10 MS Visual Studio
ifeq ($(SYSTEM),WINDOWS)
CC=cl
#CCOPT=/W4 /EHsc /O2 
CCOPT=/W4 /EHsc /c
CLOPT=
LIBPATH=
LIBO=/Fo
BINO=/Fe
LIBS=..\aaa-lib\\
BINS=..\aaa-bin\\

endif


# general rules

qcfp.$(PROJECT): $(PROJECT).obj main.obj depends
	$(CC) $(CLOPT) $(BINO)$(BINS)z.$(PROJECT) $(LIBS)$(PROJECT).obj $(LIBS)$(PROJECT).main.obj  \
	$(foreach val,$(DEPENDENCIES), $(LIBS)$(val).obj) $(LIBPATH) 

$(PROJECT).obj: $(PROJECT).cpp $(PROJECT).hpp
	$(CC) $(CCOPT) $(PROJECT).cpp $(LIBO)$(LIBS)$(PROJECT).obj $(INCPATH) 

main.obj: main.cpp $(PROJECT).hpp
	$(CC) $(CCOPT) main.cpp  $(LIBO)$(LIBS)$(PROJECT).main.obj $(INCPATH) 

depends:
	$(foreach val,$(DEPENDENCIES), $(MAKE) $(val).obj -C ../$(val) &&) echo



clean:
	rm -f $(LIBS)$(PROJECT)* $(BINS)z.$(PROJECT)*