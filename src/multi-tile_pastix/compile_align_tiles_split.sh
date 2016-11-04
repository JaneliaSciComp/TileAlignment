#!/bin/bash
echo "-------------------------------------------"
echo "-------------------------------------------"
echo "compiling pastix test example"
echo "-------------------------------------------"
echo "-------------------------------------------"

export SCOTCH_HOME=/usr/local/scotch_6.0.4
export export PASTIX_ROOT=/usr/local/pastix_5.2.2.22
LIBS="-L/home/denisovg/lib64 -L$PASTIX_ROOT/lib -lpastix -lifcore -lm -lrt -L/usr/local/scotch_6.0.4/lib -lscotch -lscotcherr -lscotcherrexit -lptscotch  -lptscotcherr  -lptscotcherrexit  -L/usr/local/hwloc-1.11.3/lib -lhwloc -lpthread -L$PASTIX_ROOT/install -L/usr/local/INTEL-2016/compilers_and_libraries_2016.2.181/linux/mpi/intel64/lib -L/usr/local/matlab-2016a/extern/lib/glnxa64 -L/usr/local/matlab-2016a/bin/glnxa64 -lstdc++ -lmatrix_driver -lmat -lmx -lmex";
INC="-I$PASTIX_ROOT/include -I$PASTIX_ROOT/src/example/src -I/usr/local/matlab-2016a/extern/include";
CC="mpicc -cc=icc -Wall -Wl,-rpath,'$ORIGIN";
CXX="mpicxx -cxx=icpc -Wall -Wl,-rpath,'$ORIGIN"
#CCOPT="-O3   -DCUDA_SM_VERSION=20 -I/usr/local/scotch_6.0.4/include -DDISTRIBUTED -DWITH_SCOTCH -I/usr/local/hwloc-1.11.3//include -DWITH_HWLOC  -DVERSION='' -DX_ARCHi686_pc_linux -DDOF_CONSTANT";
CCOPT="-g  -DCUDA_SM_VERSION=20 -I/usr/local/scotch_6.0.4/include -DDISTRIBUTED -DWITH_SCOTCH -I/usr/local/hwloc-1.11.3//include -DWITH_HWLOC  -DVERSION='' -DX_ARCHi686_pc_linux -DDOF_CONSTANT";
CXXOPT=$CCOPT;
#"_CXXOPT_";
CL="mpicc -cc=icc -Wall";
FC="mpif90 -f90=ifort -fpp";
FCOPT="";
FL="mpif90 -f90=ifort ";
OPTS="  -DPREC_DOUBLE  -DCUDA_SM_VERSION=20 -I/usr/local/scotch_6.0.4/include -DDISTRIBUTED -DWITH_SCOTCH -I/usr/local/hwloc-1.11.3/include -DWITH_HWLOC";
BLAS="-L/usr/local/BLAS-3.6.0 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core";
VERSION="exported";
LIBS_MURGE=`echo $LIBS | sed -e 's/-lpastix/-lpastix_murge -lpastix/g'`

usage="usage : $0 [options] - Shows PaStiX libs, includes and compiler\n";
usage="$usage    options : \n";
usage="$usage        --libs               - prints librairies\n";
usage="$usage        --libs_murge         - prints librairies\n";
usage="$usage        --incs               - prints includes\n";
usage="$usage        --cc                 - prints C compiler\n";
usage="$usage        --ccopts             - prints C compiler options\n";
usage="$usage        --cxx                - prints C++ compiler\n";
usage="$usage        --cxxcopts           - prints C++ compiler options\n";
usage="$usage        --cl                 - prints C linker\n";
usage="$usage        --fc                 - prints fortran compiler\n";
usage="$usage        --fcopts             - prints fortran compiler options\n";
usage="$usage        --fl                 - prints fortran linker\n";
usage="$usage        --opts               - prints PaStiX compiling options\n";
usage="$usage        --vers               - prints PaStiX version\n";
usage="$usage        --blas               - prints blas choosen in config.in\n";

if [ $# = 0 ]
then
    echo "Librairies               : $LIBS" ;
    echo "Librairies with murge    : $LIBS_MURGE";
    echo "Incs                     : $INC" ;
    echo "C Compiler               : $CC" ;
    echo "C Compiler options       : $CCOPT" ;
    echo "C++ Compiler             : $CXX" ;
    echo "C++ Compiler options     : $CXXOPT" ;
    echo "Fortran Compiler         : $FC" ;
    echo "Fortran Compiler options : $FCOPT" ;
    echo "C Linker                 : $CL" ;
    echo "Fortran Linker           : $FL" ;
    echo "Options                  : $OPTS" ;
    echo "Version                  : $VERSION" ;
    echo "Blas                     : $BLAS" ;
elif [ $# = 1 ]
then
    case $1 in
        --libs)
            echo $LIBS;;
        --libs_murge)
            echo $LIBS_MURGE;;
        --incs)
            echo $INC;;
        --cc)
            echo $CC;;
        --ccopts)
            echo $CCOPT;;
        --cxx)
            echo $CXX;;
        --cxxopts)
            echo $CXXOPT;;
        --fc)
            echo $FC;;
        --fcopts)
            echo $FCOPT;;
        --cl)
            echo $CL;;
        --fl)
            echo $FL;;
        --opts)
            echo $OPTS;;
        --blas)
            echo $BLAS;;
        --vers)
            echo $VERSION;;

        *)
            echo -e $usage
    esac;
else
    echo -e $usage
fi;

echo "$CC align_tiles_split.c $INC $LIBS $CCOPT $OPTS $BLAS  -o align_tiles_split_pastix"

mpicc align_tiles_split.c /home/denisovg/lib64/libstdc++.so.6 $INC $LIBS $CCOPT $OPTS $BLAS  -o align_tiles_split_pastix



