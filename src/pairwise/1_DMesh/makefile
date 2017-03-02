include $(TA_HOME)/src/pairwise/aln_makefile_std_defs

CFLAGS += -DCREATEMESH_WRITE_DEBUG_IMAGES

files =\
 ApproximateMatch_NoCR.cpp\
 ApproximateMatch.cpp\
 CGBL_dmesh.cpp\
 CreateMesh.cpp\
 CThmUtil.cpp\
 dmesh.cpp\
 dmeshdriver.cpp\
 ImproveMesh.cpp\
 InSectionOverlap.cpp\
 RegionToRegionMap.cpp

objs = ${files:.cpp=.o}

all : align2     

clean :
	rm -f *.o

align2: .CHECK_GENLIB ${objs}
	@mkdir -p ../../../bin
	$(CC) $(LFLAGS) ${objs} $(LINKS_STD) $(OUTPUT)

