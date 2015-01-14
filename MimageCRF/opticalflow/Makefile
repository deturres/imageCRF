SRC =  GaussianPyramid.cpp OpticalFlow.cpp Stochastic.cpp

CC = g++

WARN   =
# for linux
OPT    =  -O3 -funroll-loops --fast-math -mfpmath=sse -msse2
# for mac
# OPT = -O3 -funroll-loops
# for linux
CPPFLAGS = $(OPT) $(WARN) -I/usr/include/opencv
# for mac
# CPPFLAGS = $(OPT) $(WARN) -I/Users/chetan/libraries/OpenCV-2.4.3/include/opencv

#-L/lib -lpng -lpng12 
# for linux
LDLIBS = -L/usr/lib/x86_64-linux-gnu -lopencv_core -lopencv_highgui -lopencv_ml
# for mac
# LDLIBS = -L/Users/chetan/libraries/OpenCV-2.4.3/bin/lib/ -lopencv_core -lopencv_highgui -lopencv_ml

OBJ = $(SRC:.cpp=.o)

all: libopflow.a

libopflow.a: $(OBJ)
	rm -f libopflow.a
	ar ruc libopflow.a $(OBJ)
	ranlib libopflow.a

clean: 
	rm -f $(OBJ)
