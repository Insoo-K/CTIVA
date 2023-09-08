# Makefile
CC=g++
CFLAGS=-m64 -c -Wall -O2
LDFLAGS=-m64

MAINOBJS=main.o
OBJS=CenData.o Density.o WKM.o PHI.o OPT.o MyMath.o ExpData.o student_t_distribution.o beta_distribution.o beta_function.o gamma_function.o ln_beta_function.o ln_gamma_function.o GAIT.o
TARGET=GAIT

all: $(OBJS) $(MAINOBJS)
	$(CC) -o $(TARGET) $(LDFLAGS) $(MAINOBJS) $(OBJS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o *.so tags $(TARGET) $(TARGET_TOOL) $(TARGET_TEST) $(OBJS)


