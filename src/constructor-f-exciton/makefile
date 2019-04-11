PROJECT=constructor-f-exciton

DEPENDENCIES= toolsIO storage dvector3d eigen constants complexv dtensor3x3 dvectorNd

include ../make.inc

main-ham-1: main-ham-1.cpp $(PROJECT).hpp
	$(CC) $(CCOPT) $(INCPATH) -c main-ham-1.cpp -o ../aaa-lib/$(PROJECT).main-ham-1.o
	$(CC) $(LIBPATH) -o ../aaa-bin/z.$(PROJECT)-ham-1 ../aaa-lib/$(PROJECT).main-ham-1.o  ../aaa-lib/$(PROJECT).o  $(CLOPT)

main-ham-dips-1: main-ham-dips-1.cpp $(PROJECT).hpp
	$(CC) $(CCOPT) $(INCPATH) -c main-ham-dips-1.cpp  -o ../aaa-lib/$(PROJECT).main-ham-dips-1.o
	$(CC) $(LIBPATH) -o ../aaa-bin/z.$(PROJECT)-ham-dips-1 ../aaa-lib/$(PROJECT).main-ham-dips-1.o  ../aaa-lib/$(PROJECT).o  $(CLOPT)

main-ham-2: main-ham-2.cpp $(PROJECT).hpp
	$(CC) $(CCOPT) $(INCPATH) -c main-ham-2.cpp  -o ../aaa-lib/$(PROJECT).main-ham-2.o
	$(CC) $(LIBPATH) -o ../aaa-bin/z.$(PROJECT)-ham-2 ../aaa-lib/$(PROJECT).main-ham-2.o  ../aaa-lib/$(PROJECT).o  $(CLOPT)

main-ham-dips-2: main-ham-dips-2.cpp $(PROJECT).hpp
	$(CC) $(CCOPT) $(INCPATH) -c main-ham-dips-2.cpp   -o ../aaa-lib/$(PROJECT).main-ham-dips-2.o
	$(CC) $(LIBPATH) -o ../aaa-bin/z.$(PROJECT)-ham-dips-2 ../aaa-lib/$(PROJECT).main-ham-dips-2.o  ../aaa-lib/$(PROJECT).o  $(CLOPT)

main-ham-3: main-ham-3.cpp $(PROJECT).hpp
	$(CC) $(CCOPT) $(INCPATH) -c main-ham-3.cpp  -o ../aaa-lib/$(PROJECT).main-ham-3.o
	$(CC) $(LIBPATH) -o ../aaa-bin/z.$(PROJECT)-ham-3 ../aaa-lib/$(PROJECT).main-ham-3.o  ../aaa-lib/$(PROJECT).o  $(CLOPT)

main-ham-dips-3: main-ham-dips-3.cpp $(PROJECT).hpp
	$(CC) $(CCOPT) $(INCPATH) -c main-ham-dips-3.cpp   -o ../aaa-lib/$(PROJECT).main-ham-dips-3.o
	$(CC) $(LIBPATH) -o ../aaa-bin/z.$(PROJECT)-ham-dips-3 ../aaa-lib/$(PROJECT).main-ham-dips-3.o  ../aaa-lib/$(PROJECT).o $(CLOPT)



