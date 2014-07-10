#======================= CoMD related ===================================
#SHELL = /bin/bash
# double precision (ON/OFF)
DOUBLE_PRECISION = ON
# MPI for parallel (ON/OFF)
DO_MPI = OFF

CC = gcc
CFLAGS = -std=c99
OPTFLAGS = -g -O5
INCLUDES = 
C_LIB = -lm

LDFLAGS += ${C_LIB} ${OTHER_LIB}
CFLAGS  += ${OPTFLAGS} ${INCLUDES} ${OTHER_INCLUDE}

# clear all suffixes
.SUFFIXES:
# list only those that we use 
.SUFFIXES: .c .o

.PHONY: DEFAULT clean distclean depend

# Check for double precision
ifeq ($(DOUBLE_PRECISION), ON)
CFLAGS += -DDOUBLE
else
CFLAGS += -DSINGLE
endif

CoMD_VARIANT = CnC
#======================= CoMD related ===================================

#=======================CnC-OCR ========================================= 	
TARGET=comd.exe
CFLAGS +=-g -I$(CNCOCR_INSTALL)/include -I$(OCR_INSTALL)/include -D__OCR__  -I./

include steplist.mk
SRCS= Common.c Context.c Dispatch.c $(STEP_SRCS) cmdLineParser.c mycommand.c random.c parallel.c
SRCS += haloExchange.c initAtoms.c linkCells.c yamlOutput.c timestep.c ljForce.c decomposition.c Main.c
#SRCS=$(wildcard *.c)
OBJS=$(patsubst %.c,%.o,$(SRCS))


CFLAGS += -include CoMDTypes.h

compile: $(TARGET)

# building source files
%.o: %.c %.h
	gcc $(CFLAGS) -c $<

# linking - creating the executable
$(TARGET): $(OBJS)
	gcc $(CFLAGS) -L"$(OCR_INSTALL)/lib" \
		-L"$(CNCOCR_INSTALL)/lib" \
		$(OBJS) \
		-locr -lcncocr ${LDFLAGS} -o$@


# delete binaries
clean:
	rm -f $(OBJS) $(TARGET)

# delete binaries and scaffolding files
squeaky: clean
	rm {Context,Dispatch,Common}.[ch] steplist.mk 
	
#=======================CnC-OCR =========================================
