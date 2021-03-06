# Top level Makefile for TileAlignment source code
# # Copyright (c) 2016 Howard Hughes Medical Institute

SUBDIRS = src/pairwise src/multi-tile_pastix 

all: align2 align_multi

align2:
	cd src/pairwise; $(MAKE) all;

align_multi:
	cd src/multi-tile_pastix; $(MAKE) all;                    
	 
clean:
	cd src/pairwise; $(MAKE) clean
	cd src/multi-tile_pastix; $(MAKE) clean
