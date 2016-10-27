
include $(TA_HOME)/src/pairwise/aln_makefile_std_defs

libname = xmllib

files =\
 tinyxml.cpp\
 tinyxmlerror.cpp\
 tinyxmlparser.cpp

objs = ${files:.cpp=.o}

all : $(libname).a

clean :
	rm -f *.{a,o}

$(libname).a : ${objs}
	ar rvs $(libname).a ${objs}

