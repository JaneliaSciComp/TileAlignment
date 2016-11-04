all : compile compile_split

compile : align_tiles.c align_tiles.sh
	@mkdir -p ../../bin   
	./compile_align_tiles.sh
	cp align_tiles_pastix align_tiles.sh ../../bin

compile_split : align_tiles_split.c align_tiles_split.sh
	@mkdir -p ../../bin  
	./compile_align_tiles_split.sh
	cp align_tiles_split_pastix align_tiles_split.sh ../../bin

clean: 
	@rm -f *.o ../../bin/align_tiles*
	@rm -rf ../../bin      
.PHONY : compile compile_split clean

