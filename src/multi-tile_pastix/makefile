all : compile compile_split

compile : align_tiles.c align_tiles.sh
	@mkdir -p ../../bin   
	./compile_align_multi_tiles.sh

compile_split : align_tiles_split.c align_tiles_split.sh
	@mkdir -p ../../bin  
	./compile_align_multi_tiles_split.sh

clean: 
	@rm -f *.o ../../bin/align_multi_tiles*
	@rm -rf ../../bin      
.PHONY : compile compile_split clean

