source = $(wildcard *.c)
objects = $(source:.c=.o)

#===============================================================================
# Compiler
#===============================================================================

# Set default compiler as gnu C
CC = gcc

ifeq ($(COMPILER),gnu)
  CC = gcc
endif
ifeq ($(COMPILER),intel)
  CC = icc
endif
ifeq ($(COMPILER),pgi)
  CC = pgcc
endif
ifeq ($(COMPILER),ibm)
  CC = xlc
endif
ifeq ($(COMPILER),cray)
  CC = cc
endif

#===============================================================================
# Targets
#===============================================================================

all: $(objects)
clean:
	@rm -f *.o

#===============================================================================
# Rules
#===============================================================================

.SUFFIXES: .c .o
.PHONY: all clean
.PRECIOUS: %.c

%.o: %.c
	$(CC) $(CCFLAGS) -c $<
