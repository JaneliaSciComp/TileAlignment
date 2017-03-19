all : compile copy_scripts

compile: align_multi.c align_multi_nosplit.c main.c main.h compile_all.sh
	@mkdir -p ../../bin   
	./compile_all.sh
copy_scripts: align_multi.sh align_multi_nosplit.sh run_align_multi.sh run_align_multi_nosplit.sh
	@mkdir -p ../../bin
	cp     align_multi.sh         ../../bin
	cp run_align_multi.sh         ../../bin
	cp     align_multi_nosplit.sh ../../bin
	cp run_align_multi_nosplit.sh ../../bin
clean: 
	@rm -f *.o ../../bin/align_multi* ../../bin/run_align_multi*
	@rm -rf ../../bin      
.PHONY : compile copy_scripts clean

