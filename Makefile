# Top level Makefile for TileAlignment source code
# # Copyright (c) 2016 Howard Hughes Medical Institute

SUBDIRS = src/pairwise src/multi-tile_pastix 

all:
	for dir in $(SUBDIRS); \
            do (cd $$dir; $(MAKE) $@); \
        done
	 
clean:
	cd src/pairwise; $(MAKE) clean
	cd src/multi-tile_pastix; $(MAKE) clean
