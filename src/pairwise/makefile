
include ./aln_makefile_std_defs

0dirs =\
 0_GEN\
 0_UTL\
 0_XML
1dirs =\
 1_DMesh\

all : xdirs = ${0dirs} ${1dirs} 
all : build_dirs

clean : xdirs = ${0dirs} ${1dirs} 
clean : clean_dirs clean_libs

build_dirs : .CHECK_ALL_LIBS
	@ \
	for i in $(xdirs); \
	do \
		cd $$i && $(MAKE) && cd ..; \
	done

clean_dirs :
	@ \
	for i in $(xdirs); \
	do \
		cd $$i && $(MAKE) clean && cd ..; \
	done

clean_libs:
	@ \
	for i in 0_GEN; \
	do \
		cd $$i && $(MAKE) clean && cd ..; \
	done

.PHONY : build_dirs clean_dirs clean_libs

