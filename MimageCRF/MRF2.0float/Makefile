# OS_TYPE = linux32
# OS_TYPE = linux64
OS_TYPE = macosx

ifeq ($(OS_TYPE),macosx)
	MATH_OPT = -funroll-loops
else ifeq ($(OS_TYPE),linux32)
	MATH_OPT = -funroll-loops --fast-math -mfpmath=sse -msse2
else ifeq ($(OS_TYPE),linux64)
	MATH_OPT = -funroll-loops --fast-math -mfpmath=sse -msse2
endif

VERSION = MRF1.9

QF = QPBO-v1.3.src

SRC =  mrf.cpp ICM.cpp GCoptimization.cpp graph.cpp maxflow.cpp \
       MaxProdBP.cpp LinkedBlockList.cpp regions-maxprod.cpp \
       TRW-S.cpp BP-S.cpp MaxProdBPTree.cpp DualDecomposition.cpp \
       MRFWrapper.cpp Flipper.cpp FlipperForDD.cpp DualDecompositionFlipper.cpp \
       FlipperGeneral.cpp \
       MaxProdBPTreeDD.cpp DualDecompositionSubGraph.cpp DualDecompositionTree.cpp \
       DualDecompositionTreeES.cpp DualDecompositionFlipperES.cpp \
       $(QF)/QPBO.cpp $(QF)/QPBO_maxflow.cpp $(QF)/QPBO_extra.cpp $(QF)/QPBO_postprocessing.cpp       

JSRC = LabelImageWriter.cpp MRFEnergy.cpp ArrayMath.cpp NodeMessages.cpp \
	   SyncSumProd.cpp AsyncSumProd.cpp SparseSyncSumProd.cpp SparseAsyncSumProd.cpp \
	   MeanField.cpp MeanField.cpp SparseMeanField.cpp SparseMeanField.cpp

CC = g++

WARN   = -W -Wall
OPT    = -std=c++11 -O3 $(MATH_OPT)
CPPFLAGS = $(OPT) $(WARN) -I../imageLib -DUSE_64_BIT_PTR_CAST
LDLIBS = -L../imageLib -lImg -lpng

OBJ = $(SRC:.cpp=.o)
JOBJ = $(JSRC:.cpp=.o)
JHEAD = $(JSRC:.cpp=.h)

all: libMRF.a example_dd_flipper_es example_dd_tree_es example_dd_segmentation example_dd_realworld2 example_dd_flipper example_dd_realworld example example_flipper example_flipper2 example_tree example_gc example_dd_mix example_dd_mix2 example_dd_nonsub example_dd_tree example_dd_subgraph example_compare
jall: $(JOBJ)

libMRF.a: $(OBJ) $(JOBJ)
	rm -f libMRF.a
	ar ruc libMRF.a $(OBJ) $(JOBJ)
	ranlib libMRF.a

example_dd_nonsub: libMRF.a example_dd_nonsub.o
	$(CC) $(OPT) -o example_dd_nonsub example_dd_nonsub.o -L. -lMRF

example_dd_realworld: libMRF.a example_dd_realworld.o
	$(CC) $(OPT) -o example_dd_realworld example_dd_realworld.o $(LDLIBS) -L. -lMRF

example_dd_segmentation: libMRF.a example_dd_segmentation.o
	$(CC) $(OPT) -o example_dd_segmentation example_dd_segmentation.o $(LDLIBS) -L. -lMRF

example_dd_realworld2: libMRF.a example_dd_realworld2.o
	$(CC) $(OPT) -o example_dd_realworld2 example_dd_realworld2.o $(LDLIBS) -L. -lMRF

example_flipper: libMRF.a example_flipper.o
	$(CC) $(OPT) -o example_flipper example_flipper.o $(LDLIBS) -L. -lMRF

example_flipper_general: libMRF.a example_flipper_general.o
	$(CC) $(OPT) -o example_flipper_general example_flipper_general.o $(LDLIBS) -L. -lMRF

enumerate_flippers: libMRF.a enumerate_flippers.o
	$(CC) $(OPT) -o enumerate_flippers enumerate_flippers.o $(LDLIBS) -L. -lMRF

example_dd_flipper: libMRF.a example_dd_flipper.o
	$(CC) $(OPT) -o example_dd_flipper example_dd_flipper.o $(LDLIBS) -L. -lMRF

example_flipper2: libMRF.a example_flipper2.o
	$(CC) $(OPT) -o example_flipper2 example_flipper2.o $(LDLIBS) -L. -lMRF

example_gc: libMRF.a example_gc.o
	$(CC) $(OPT) -o example_gc example_gc.o -L. -lMRF

example_dd_subgraph: libMRF.a example_dd_subgraph.o
	$(CC) $(OPT) -o example_dd_subgraph example_dd_subgraph.o -L. -lMRF

example_dd_mix: libMRF.a example_dd_mix.o
	$(CC) $(OPT) -o example_dd_mix example_dd_mix.o -L. -lMRF

example_dd_mix2: libMRF.a example_dd_mix2.o
	$(CC) $(OPT) -o example_dd_mix2 example_dd_mix2.o -L. -lMRF

example_dd_mix_es: libMRF.a example_dd_mix_es.o
	$(CC) $(OPT) -o example_dd_mix_es example_dd_mix_es.o -L. -lMRF

example_dd_tree: libMRF.a example_dd_tree.o
	$(CC) $(OPT) -o example_dd_tree example_dd_tree.o -L. -lMRF

example_dd_tree_es: libMRF.a example_dd_tree_es.o
	$(CC) $(OPT) -o example_dd_tree_es example_dd_tree_es.o -L. -lMRF

example_dd_flipper_es: libMRF.a example_dd_flipper_es.o
	$(CC) $(OPT) -o example_dd_flipper_es example_dd_flipper_es.o -L. -lMRF

example_tree: libMRF.a example_tree.o
	$(CC) $(OPT) -o example_tree example_tree.o -L. -lMRF

example: libMRF.a example.o
	$(CC) -o example example.o -L. -lMRF

example_compare: libMRF.a example_compare.o
	$(CC) -o example_compare example_compare.o -L. -lMRF

test: libMRF.a testSumProd.o
	$(CC) -o testSumProd testSumProd.o -L. -lMRF $(LDLIBS)

clean:
	rm -f $(OBJ) example.o example_dd_mix_es.o example_dd_flipper_es.o example_dd_tree_es.o example_dd_segmentation.o example_dd_realworld2.o example_dd_realworld.o example_dd_flipper.o example_flipper2.o example_flipper.o example_gc.o example_tree.o example_tree2.o example_dd_mix.o example_dd_mix2.o example_dd_nonsub.o example_dd.o example_dd_tree.o example_dd_subgraph.o example_compare.o

jclean:
	rm -f $(JOBJ)


allclean: clean jclean
	rm -f libMRF.a example example_dd_mix_es example_dd_flipper_es example_dd_tree_es example_dd_segmentation example_dd_realworld2 example_dd_flipper example_flipper example_flipper2 example_dd_realworld example_gc example_tree example_tree2 example_dd_mix example_dd_mix2 example_dd_nonsub example_dd example_dd_tree example_dd_subgraph example_compare

zip:
	cd .. ;	rm -f $(VERSION).zip ; zip $(VERSION).zip \
	$(VERSION)/*txt $(VERSION)/Makefile $(VERSION)/SConstruct \
	$(VERSION)/*cpp $(VERSION)/*h \
	$(VERSION)/MRFvcc6.dsw $(VERSION)/MRFvcc6.dsp \
	$(VERSION)/MRFvcc7/MRFvcc7.sln $(VERSION)/MRFvcc7/MRFvcc7.vcproj
