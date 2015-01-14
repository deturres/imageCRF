# OS_TYPE = linux32
OS_TYPE = linux64
# OS_TYPE = macosx

ifeq ($(OS_TYPE),macosx)
	OPENCV_INCLUDE = /Users/chetan/libraries/OpenCV-2.4.3/include/opencv
	OPENCV_LIB = /Users/chetan/libraries/OpenCV-2.4.3/bin/lib/
	MATH_OPT =
else ifeq ($(OS_TYPE),linux32)
	OPENCV_INCLUDE = /usr/include/opencv
	OPENCV_LIB = /usr/lib
	MATH_OPT = -funroll-loops --fast-math -mfpmath=sse
else ifeq ($(OS_TYPE),linux64)
	OPENCV_INCLUDE = /usr/include/opencv
	OPENCV_LIB = /usr/lib/x86_64-linux-gnu
	MATH_OPT = -funroll-loops --fast-math -mfpmath=sse
endif

SRC = evaluation.cpp Global.cpp crf.cpp bfgs.cpp baseModel.cpp readParamFiles.cpp helperLib.cpp features.cpp featuresM.cpp intensityNB.cpp intensityClass.cpp appearanceMV.cpp locationMV.cpp appearance.cpp crfmodel.cpp crosscheck.cpp evaldisps.cpp logreg.cpp logregpair.cpp logregpl.cpp logregpd.cpp logregpd4.cpp crf.cpp IO.cpp modelmain.cpp expectations.cpp Parameters.cpp

BIN = modelmain

MRF = MRF2.0float
OF = opticalflow

CC = g++

WARN   = -w -W -Wall

CPPFLAGS = -O3 -DNDEBUG $(MATH_OPT) $(WARN) -I$(MRF) -IimageLib -Iboost_1_41_0 -Ilibconfig-1.4.6/lib -I$(OPENCV_INCLUDE)

LDLIBS = -L$(OPENCV_LIB) -lopencv_core -lopencv_highgui -lopencv_ml -L$(MRF) -lMRF -LimageLib -lImg -L$(OPENCV_LIB) -lpng -Lzlib-1.2.3 -lz -Llibconfig-1.4.6/lib/.libs -lconfig++ -L$(OF) -lopflow

OBJ = $(SRC:.cpp=.o)
EXE = $(SRC:.cpp=.exe)

all: $(BIN)

modelmain: modelmain.o Global.o evaluation.o crf.o bfgs.o baseModel.o readIVMparam.o logreg.o logregpl.o logregpair.o logregpd.o logregpd4.o crf.o featuresM.o Parameters.o readParamFiles.o helperLib.o crfmodel.o expectations.o evaldisps.o appearance.o locationMV.o appearanceMV.o intensityNB.o intensityClass.o IO.o features.o

baseModel.o: baseModel.h
crf.o: crf.h
bfgs.o: bfgs.h
logreg.o: logreg.h
logregpl.o: logregpl.h
featuresM.o:featuresM.h
crfmodel.o: crfmodel.h
Global.o: Global.h
modelmain.o: modelmain.h
Parameters.o: Parameters.h


clean: 
	rm -f $(OBJ) core *.stackdump

allclean: clean
	rm -f $(BIN) $(EXE)
