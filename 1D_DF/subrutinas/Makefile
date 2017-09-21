FC=gfortran
CFLAG=-c

SOURCE=$(shell ls | grep .*90)
OBJS=$(patsubst %.f90, %.o , $(wildcard *.f90))

all:$(OBJS)

$(OBJS):$(SOURCE)
	$(FC) $(CFLAG) $(SOURCE)

clean:
	rm *.o
