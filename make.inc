# UNIFIED MAKEFILE FOR QCFP
# NOT USED IN contructor-excitons

MPI_COMPILE=FALSE


# set the compilers
ifeq ($(MPI_COMPILE),TRUE)

CCMPI=/Users/dariusa/usr/bin/mpicxx
CCMPIOPT=-O2 -DMPIPROCESSING
#CCMPIOPT=--debug

endif

CC=llvm-g++
#CC=g++-5
#CCOPT=-O2
CCOPT=--debug


# FOR MACOSX AND GNU
CLOPT= -lfftw3 -llapack -lblas 


INCPATH= -I ~/usr/include 
LIBPATH= -L ~/usr/lib 



# general rules

qcfp.$(PROJECT): $(PROJECT).o main.o depends
	$(CC) $(LIBPATH) -o ../aaa-bin/z.$(PROJECT) ../aaa-lib/$(PROJECT).o ../aaa-lib/$(PROJECT).main.o  \
	$(foreach val,$(DEPENDENCIES), ../aaa-lib/$(val).o)  $(CLOPT)

$(PROJECT).o: $(PROJECT).cpp $(PROJECT).hpp
	$(CC) $(CCOPT)  $(INCPATH) -c $(PROJECT).cpp  -o ../aaa-lib/$(PROJECT).o

main.o: main.cpp $(PROJECT).hpp
	$(CC) $(CCOPT)  $(INCPATH) -c main.cpp  -o ../aaa-lib/$(PROJECT).main.o

depends:
	$(foreach val,$(DEPENDENCIES), $(MAKE) $(val).o -C ../$(val) &&) echo



clean:
	rm -f ../aaa-lib/$(PROJECT)* ../aaa-bin/z.$(PROJECT)*