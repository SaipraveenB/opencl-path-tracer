# MAKEFILE

UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
# do something Linux-y
	LDFLAGS = -L/opt/AMDAPP/lib/x86_64 -lOpenCL
	CCFLAGS =
endif

ifeq ($(UNAME), Darwin)
# do something Solaris-y
	LDFLAGS = -framework OpenCL
	CCFLAGS = -framework OpenCL
endif

all: test

test: main.o
	g++ -o $@ $^ $(LDFLAGS)

main.o: main.cpp
	g++ -c $< $(CCFLAGS)

.PHONY: all clean

clean:
	rm *.o
	rm test
