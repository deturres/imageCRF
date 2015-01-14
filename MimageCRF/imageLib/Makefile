WARN   = -W -Wall
CPPFLAGS = -O2 $(WARN) -I/usr/include

SRC = Convert.cpp Convolve.cpp Image.cpp ImageIO.cpp ImageIOpng.cpp RefCntMem.cpp
OBJ    = $(SRC:.cpp=.o)

libImg.a: $(OBJ)
	rm -f libImg.a
	ar ruc libImg.a $(OBJ)
	ranlib libImg.a

clean: 
	rm -f $(OBJ)

allclean: clean
	rm -f libImg.a

zip:
	cd .. ;	rm -f imageLib.zip ; zip imageLib.zip \
	imageLib/README.txt imageLib/Makefile \
	imageLib/*cpp imageLib/*h \
	imageLib/imageLib.sln imageLib/imageLib.vcproj \
	imageLib/lpng/*
