

SUNDIALS_DIR=/usr/local
SUNDIALS_LIB_DIR=$(SUNDIALS_DIR)/lib
SUNDIALS_INC_DIR=$(SUNDIALS_DIR)/include
SUNDIALS_LIBS=-lsundials_cvode -lsundials_nvecserial
SUNDIALS_INCS=-I$(SUNDIALS_INC_DIR)
LIBS=-lm

all: demain

demain: demain.o de.o
	$(CC) $(LDFLAGS) -o demain demain.o de.o -L$(SUNDIALS_LIB_DIR) $(SUNDIALS_LIBS) $(LIBS)

demain.o: demain.c de.h
	$(CC) $(CPPFLAGS) $(SUNDIALS_INCS) -c demain.c

de.o: de.c de.h
	$(CC) $(CPPFLAGS) $(SUNDIALS_INCS) -c de.c

clean:
	rm -f demain demain.o de.o

