all : compile compile_split copy_scripts

compile : align_tiles.c align_tiles.sh
	@mkdir -p ../../bin   
	./compile_align_tiles.sh

compile_split : align_tiles_split.c align_tiles_split.sh
	@mkdir -p ../../bin  
	./compile_align_tiles_split.sh

copy_scripts: align_tiles.sh align_tiles_split.sh run_align_tiles.sh run_align_tiles_split.sh
	@mkdir -p ../../bin
	cp     align_tiles.sh       ../../bin
	cp run_align_tiles.sh       ../../bin
	cp     align_tiles_split.sh ../../bin
	cp run_align_tiles_split.sh ../../bin
clean: 
	@rm -f *.o ../../bin/align_tiles* ../../bin/run_align_tiles*
	@rm -rf ../../bin      
.PHONY : compile compile_split copy_scripts clean

