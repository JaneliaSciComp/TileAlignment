

#include	"CGBL_dmesh.h"

#include	"Cmdline.h"
#include	"File.h"
#include	"CAffineLens.h"
#include	"Debug.h"

#include        <string.h>
#include        <curl/curl.h>
#include        <sstream>

#include        "png.h"
#include        <iostream>

/* --------------------------------------------------------------- */
/* Constants ----------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

CGBL_dmesh	GBL;

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */






/* --------------------------------------------------------------- */
/* Object management --------------------------------------------- */
/* --------------------------------------------------------------- */

CGBL_dmesh::CGBL_dmesh()
{
	_arg.SCALE			= 999.0;
	_arg.XSCALE			= 999.0;
	_arg.YSCALE			= 999.0;
	_arg.SKEW			= 999.0;
	_arg.ima			= NULL;
	_arg.imb			= NULL;
	_arg.FLD			= 0;
	_arg.MODE			= 0;

	arg.CTR				= 999.0;
	arg.fma				= NULL;
	arg.fmb				= NULL;
	arg.Transpose		= false;
	arg.WithinSection	= false;
	arg.Verbose			= false;
	arg.SingleFold		= false;
	arg.Heatmap			= false;

	A.z	=  0;
	A.id	= -1;

	B.z	=  0;
	B.id	= -1;

        options.use_idb           = true;
        options.pts_file          = NULL;
}

/* --------------------------------------------------------------- */
/* ShowUsage  ---------------------------------------------------- */
/* --------------------------------------------------------------- */

void CGBL_dmesh::ShowUsage(int argc, char *argv[])
{
    fprintf( stderr,
    "Usage: %s argument(s) [ options ], where\n"
    "    Where argument(s) == one of\n"
    "           <imageA_file> <imageB_file> // read inputs from local files\n"
    "           <imageA_URL>  <imageB_URL>  // read inputs from URLs\n"
    "                                       \n"
    "    If argument(s) ==  <imageA_file> <imageB_file>\n"
    "    then options == \n"
    "      -h \n"
    "      -Tab=<six comma-separated values or URL>\n"
    "      -matchparams_file=<path>\n"
    "      -pts_file=<path>\n"
    "      -comp_file=<path>\n"
    "      -registered_file=<path>\n"
    "                             \n"
    "    If argument(s) ==  <imageA_URL> <imageB_URL>\n"
    "    then options == \n"
    "      -h \n"
    "      -Ta=<URL for affine coefficients for imageA>\n"
    "      -Tb=<URL for affine coefficients for imageB>\n"
    "      -matchparams_file=<path>\n"
    "      -pts_file=<path>\n"
    "      -comp_file=<path>\n"
    "      -registered_file=<path>\n"
             , argv[0]);
}

/* --------------------------------------------------------------- */

string extract_PNG_from_URL(string my_url)
{
    // Handling image extracted from URL and then stored in PNG file on a local drive
    // (using libcurl and libpng)
    // http://www.labbookpages.co.uk/software/imgProc/libPNG.html
    // http://stackoverflow.com/questions/12728524/save-an-image-from-jpeg-to-png-using-libcurl-and-gdkpixbuff
    //
    CURL *easy_handle;
    CURLcode imgresult;

    string suffix(".png");
    string output_path = tmpnam(NULL) + suffix;

    int code = 0;
    FILE *fp = NULL;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    png_bytep row = NULL;

    struct context
    {
        unsigned char *data;
        int allocation_size;
        int length;
    };
    struct context ctx;

    // Initialize write structure
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL) {
        fprintf(stderr, "Could not allocate write struct\n");
        code = 1;
        goto finalise;
    }

    // Initialize info structure
    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL) {
        fprintf(stderr, "Could not allocate info struct\n");
        code = 1;
        goto finalise;
    }

    // Setup Exception handling
    if (setjmp(png_jmpbuf(png_ptr))) {
        fprintf(stderr, "Error during png creation\n");
        code = 1;
        goto finalise;
    }

    easy_handle = curl_easy_init();
    if ( easy_handle ){
        // Open file
        fp = fopen((const char *)output_path.c_str(), "wb");
            if (fp == NULL) {
                (void)perror("The following error occurred");
                return NULL;
            }

        curl_easy_setopt(easy_handle, CURLOPT_URL, my_url.c_str());
        curl_easy_setopt(easy_handle, CURLOPT_WRITEFUNCTION, NULL);
        curl_easy_setopt(easy_handle, CURLOPT_WRITEDATA, fp);

        // Grab image
        imgresult = curl_easy_perform(easy_handle);
        if( imgresult ){
            cout << "Cannot grab the image!\n";
        }
    }
    finalise:
    curl_easy_cleanup(easy_handle);
    // Close the file
    fclose(fp);
    return output_path;
}


/* --------------------------------------------------------------- */
/* SetCmdLine ---------------------------------------------------- */
/* --------------------------------------------------------------- */

bool CGBL_dmesh::SetCmdLine( int argc, char* argv[] )
{
    if (argc == 1) {
        ShowUsage(argc, argv);
        return false;
    }

// Parse args

        vector<string>   key;
	vector<double>   vD;
        vector<double>   vD_url;
        vector<TAffine>  Ta;
        vector<TAffine>  Tb;
        const char      *matchparams_path = NULL;
        const char      *pts_file  = NULL;
        const char      *comp_file = NULL;
        const char      *registered_file = NULL;

        // Get the keys
        for( int i = 1; i < argc; ++i ) 
            if( argv[i][0] != '-' )
                        key.push_back(argv[i]);

        // Make sure the input files exist
        if (key.size() == 2) {
            FILE *fk;
            if (!strstr( key[0].c_str(), "http" )) {
                FILE *fk = fopen(key[0].c_str(), "r");
                if (fk == NULL) {
                    fprintf(stderr, "Input image file %s does not exist\n", key[0].c_str());
                    return false;
                }
                else
                    fclose(fk);
            }

            if (!strstr( key[1].c_str(), "http" )) {
                fk = fopen(key[1].c_str(), "r");
                if (fk == NULL) {
                    fprintf(stderr, "Input image file %s does not exist\n", key[1].c_str());
                    return false;
                }
                else
                    fclose(fk);
            }
        }

        FILE *flog;
        if (key.size() < 2)           
            flog = stdout;
        else {
            flog = stderr;
        }

        // Parse command line options
	for( int i = 1; i < argc; ++i ) {
                if( argv[i][0] != '-' )
                    ;
		else if( GetArgList( vD, "-Tdfm=", argv[i] ) ) {

			if( 6 == vD.size() )
				_arg.Tdfm.push_back( TAffine( &vD[0] ) );
			else {
				fprintf(flog,  
				"main: WARNING: Bad format in -Tdfm [%s].\n",
				argv[i] );
			}
		}
		else if( GetArgList( vD, "-Tab=", argv[i] ) ) {

			if( 6 == vD.size() )
				_arg.Tab.push_back( TAffine( &vD[0] ) );
			else {
				fprintf(flog,
				"main: WARNING: Bad format in -Tab [%s].\n",
				argv[i] );
			}
		}
                else if( GetArgList( vD, "-Ta=", argv[i] ) ) {

                        if( 6 == vD.size() )
                                Ta.push_back( TAffine( &vD[0] ) );
                        else {
                                fprintf(flog,
                                "main: WARNING: Bad format in -Ta [%s].\n",
                                argv[i] );
                        }
                }
                else if( GetArgList( vD, "-Tb=", argv[i] ) ) {

                        if( 6 == vD.size() )
                                Tb.push_back( TAffine( &vD[0] ) );
                        else {
                                fprintf(flog,
                                "main: WARNING: Bad format in -Tb [%s].\n",
                                argv[i] );
                        }
                }
                else if( GetArgListFromURL( vD, "-Tb=", argv[i] ) ) {
                        if( 6 == vD.size() )
                                Tb.push_back( TAffine( &vD[0] ) );
                        else {
                                fprintf(flog,
                                "main: WARNING: Bad format in -Tb [%s].\n",
                                argv[i] );
                        }
                }
                else if( GetArgListFromURL( vD, "-Ta=", argv[i] ) ) {
                        if( 6 == vD.size() )
                                Ta.push_back( TAffine( &vD[0] ) );
                        else {
                                fprintf(flog,
                                "main: WARNING: Bad format in -Ta [%s].\n",
                                argv[i] );
                        }
                }
                else if( GetArgList( vD, "-z=", argv[i] ) ) {
                        if( 2 == vD.size() ) {
                                A.z = vD[0];
                                B.z = vD[1];
                        } else {
                                fprintf(flog,
                                "main: option -z should specify two argument values [%s].\n",
                                argv[i] );
                        }
                }
		else if( GetArg( &_arg.SCALE, "-SCALE=%lf", argv[i] ) )
			;
		else if( GetArg( &_arg.XSCALE, "-XSCALE=%lf", argv[i] ) )
			;
		else if( GetArg( &_arg.YSCALE, "-YSCALE=%lf", argv[i] ) )
			;
		else if( GetArg( &_arg.SKEW, "-SKEW=%lf", argv[i] ) )
			;
		else if( GetArgStr( _arg.ima, "-ima=", argv[i] ) )
			;
		else if( GetArgStr( _arg.imb, "-imb=", argv[i] ) )
			;
                else if( GetArgStr(matchparams_path, "-matchparams_file=", argv[i] ) )
                        ;                   
                else if( GetArgStr(pts_file, "-pts_file=", argv[i] ) )
                        ;
                else if( GetArgStr(comp_file,"-comp_file=", argv[i] ) )
                        ;
                else if( GetArgStr(registered_file, "-registered_file=", argv[i] ) )
                        ;
		else if( GetArg( &_arg.FLD, "-FLD=%c", argv[i] ) )
			;
		else if( GetArg( &_arg.MODE, "-MODE=%c", argv[i] ) )
			;
		else if( GetArg( &arg.CTR, "-CTR=%lf", argv[i] ) )
			;
		else if( GetArgStr( arg.fma, "-fma=", argv[i] ) )
			;
		else if( GetArgStr( arg.fmb, "-fmb=", argv[i] ) )
			;
		else if( IsArg( "-tr", argv[i] ) )
			arg.Transpose = true;
		else if( IsArg( "-ws", argv[i] ) )
			arg.WithinSection = true;
		else if( IsArg( "-nf", argv[i] ) )
			_arg.FLD = 'N';
		else if( IsArg( "-sf", argv[i] ) )
			arg.SingleFold = true;
		else if( IsArg( "-v", argv[i] ) )
			arg.Verbose = true;
		else if( IsArg( "-heatmap", argv[i] ) )
			arg.Heatmap = true;
		else if( IsArg( "-dbgcor", argv[i] ) )
			dbgCor = true;
		else if( GetArgList( vD, "-Tmsh=", argv[i] ) ) {

			if( 6 == vD.size() )
				Tmsh.push_back( TAffine( &vD[0] ) );
			else {
				fprintf(flog,
				"main: WARNING: Bad format in -Tmsh [%s].\n",
				argv[i] );
			}
		}
		else if( GetArgList( vD, "-XYexp=", argv[i] ) ) {

			if( 2 == vD.size() )
				XYexp.push_back( Point( vD[0], vD[1] ) );
			else {
				fprintf(flog,
				"main: WARNING: Bad format in -XYexp [%s].\n",
				argv[i] );
			}
		}
		else {
			fprintf(flog,"Did not understand option '%s'.\n", argv[i] );
			return false;
		}
	}

        // Make sure the input parameters file exists   
        if (key.size() == 2) { 
            if (matchparams_path == NULL) {
                fprintf(stderr, "Input file matchparams.txt has not been specified\n");
                return false;
            } else {
                FILE *fm = fopen(matchparams_path, "r");
                if (fm == NULL) {
                    fprintf(stderr, "Input parameter file %s does not exist\n", matchparams_path);
                    return false;
                } else
                    fclose(fm);
            }
        }

        if (A.z < 0 || B.z < 0) {
            fprintf(flog, "\nPlease, specify non-negative z-values for images A and B: -z=<A.z>,<B.z> \n");
            return false;
        }

        // Decode labels in key
	if ( key.size() == 0 ||
	    (key.size() == 1 &&  
             4 != sscanf( key[0].c_str(), "%d.%d^%d.%d", &A.z, &A.id, &B.z, &B.id ))) 
        {
                ShowUsage(argc, argv);
		return false;
	}
        else if (key.size() == 2) {
            // Image file or URL input
            options.use_idb = false;
            if (matchparams_path != NULL) {
                options.matchparams_path = matchparams_path;
            }
            if (pts_file != NULL) 
                options.pts_file = pts_file;

            if (comp_file != NULL) {
                options.comp_file = comp_file;
            }
            if (registered_file != NULL) {
                options.registered_file = registered_file;
            }
        } 
 
// Rename stdout using image labels
        if (A.id >= 0 && B.id >= 0) 
            OpenPairLog( A.z, A.id, B.z, B.id );
        
	fprintf(flog, "\n---- dmesh start ----\n" );

// Record start time

	time_t	t0 = time( NULL );
	fprintf(flog, "main: Start: %s\n", ctime(&t0) );

// Get default parameters

	if( !ReadMatchParams( mch, A.z, B.z, options.matchparams_path, flog) )
		return false;

// Which file params to use according to (same,cross) layer

	double	cSCALE=1, cXSCALE=1, cYSCALE=1, cSKEW=0;
	int		cDfmFromTab;

	ctx.FLD = mch.FLD;

	if( A.z == B.z ) {

		cDfmFromTab	= mch.TAB2DFM_SL;

		//ctx.Tdfm = identity (default)
		ctx.XYCONF	= mch.XYCONF_SL;
		ctx.NBMXHT	= mch.NBMXHT_SL;
		ctx.HFANGDN	= mch.HFANGDN_SL;
		ctx.HFANGPR	= mch.HFANGPR_SL;
		ctx.RTRSH	= mch.RTRSH_SL;
		ctx.RIT		= mch.RIT_SL;
		ctx.RFA		= mch.RFA_SL;
		ctx.RFT		= mch.RFT_SL;
		ctx.OLAP2D	= mch.OLAP2D_SL;
		ctx.MODE	= mch.MODE_SL;
		ctx.THMDEC	= mch.THMDEC_SL;
		ctx.OLAP1D	= mch.OLAP1D_SL;
		ctx.LIMXY	= mch.LIMXY_SL;
		ctx.OPT		= mch.OPT_SL;
	}
	else {

		cSCALE	= mch.SCALE;
		cXSCALE	= mch.XSCALE;
		cYSCALE	= mch.YSCALE;
		cSKEW	= mch.SKEW;

		ctx.Tdfm.ComposeDfm( cSCALE, cXSCALE, cYSCALE, 0, cSKEW );

		cDfmFromTab	= mch.TAB2DFM_XL;

		ctx.XYCONF	= mch.XYCONF_XL;
		ctx.NBMXHT	= mch.NBMXHT_XL;
		ctx.HFANGDN	= mch.HFANGDN_XL;
		ctx.HFANGPR	= mch.HFANGPR_XL;
		ctx.RTRSH	= mch.RTRSH_XL;
		ctx.RIT		= mch.RIT_XL;
		ctx.RFA		= mch.RFA_XL;
		ctx.RFT		= mch.RFT_XL;
		ctx.OLAP2D	= mch.OLAP2D_XL;
		ctx.MODE	= mch.MODE_XL;
		ctx.THMDEC	= mch.THMDEC_XL;
		ctx.OLAP1D	= mch.OLAP1D_XL;
		ctx.LIMXY	= mch.LIMXY_XL;
		ctx.OPT		= true;
	}

// Fetch Til2Img entries
	fprintf(flog, "\n---- Input images ----\n" );

        if (A.id >= 0 && B.id >= 0) { 
             // Read inputs from an image database
             options.use_idb = true;
             IDBFromTemp( idb, "../../");  // get a path to idb

	    if( !IDBT2IGet1( A.t2i, idb, A.z, A.id, _arg.ima, flog ) ||
		!IDBT2IGet1( B.t2i, idb, B.z, B.id, _arg.imb, flog ) ) {

		return false;
	    }
        } else {
            // Read the inputs from the specified files or URLs
            options.use_idb = false;
            A.t2i.id   = A.id; 
            A.t2i.T    = TAffine( 1,0,Ta[0].t[2],0,1,Ta[0].t[5] ); // 1st image is considered as a reference
            A.t2i.col  = -999;
            A.t2i.row  = -999;            
            A.t2i.cam  = 0;
            if (!strstr( key[0].c_str(), "http" ))  // inputs from files
                A.t2i.path = key[0];
            else
                A.t2i.path = extract_PNG_from_URL(key[0]);

            B.t2i.id   = A.id;
            B.t2i.T    = TAffine( 1,0,Tb[0].t[2],0,1,Tb[0].t[5]);
            B.t2i.col  = -999;            
            B.t2i.row  = -999;            
            B.t2i.cam  = 0;
            if (!strstr( key[1].c_str(), "http" ))  // inputs from files
                B.t2i.path = key[1];
            else
                B.t2i.path = extract_PNG_from_URL(key[1]);


            // Inputs from files ( path is type "string")
//          std::cout << "Translation coefficients= 1,0," << Tb[0].t[2]-Ta[0].t[2] << ",0,1," << Tb[0].t[5]-Ta[0].t[5] << "\n";
        }
	PrintTil2Img( flog, 'A', A.t2i );
	PrintTil2Img( flog, 'B', B.t2i );

	fprintf(flog, "\n" );

// Commandline overrides

	fprintf(flog, "\n---- Command-line overrides ----\n" );


	if( _arg.Tab.size() ) {

		Tab = _arg.Tab[0];

		// remove lens parts of Tab coming from cross_thisblock

		if( mch.PXLENS && A.z != B.z ) {

			CAffineLens	LN;

			if( !LN.ReadIDB( idb ) )
				return false;

			LN.UpdateTFormRHS( Tab, A.t2i.cam, true );
			LN.UpdateTFormLHS( Tab, B.t2i.cam, false );
		}

		Tab.TPrint( flog, "Tab= " );
	}
	else
		Tab.FromAToB( A.t2i.T, B.t2i.T );

	int	altTdfm = false;

	if( _arg.Tdfm.size() ) {

		ctx.Tdfm	= _arg.Tdfm[0];
		altTdfm		= true;
	}
	else {

		if( _arg.SCALE != 999.0 ) {
			cSCALE	= _arg.SCALE;
			altTdfm	= true;
			fprintf(flog, "SCALE=%g\n", _arg.SCALE );
		}

		if( _arg.XSCALE != 999.0 ) {
			cXSCALE	= _arg.XSCALE;
			altTdfm	= true;
			fprintf(flog, "XSCALE=%g\n", _arg.XSCALE );
		}

		if( _arg.YSCALE != 999.0 ) {
			cYSCALE	= _arg.YSCALE;
			altTdfm	= true;
			fprintf(flog, "YSCALE=%g\n", _arg.YSCALE );
		}

		if( _arg.SKEW != 999.0 ) {
			cSKEW	= _arg.SKEW;
			altTdfm	= true;
			fprintf(flog, "SKEW=%g\n", _arg.SKEW );
		}

		if( altTdfm )
			ctx.Tdfm.ComposeDfm( cSCALE, cXSCALE, cYSCALE, 0, cSKEW );
	}

	if( !altTdfm && cDfmFromTab ) {

		TAffine	R;
		R.NUSetRot( -Tab.GetRadians() );

		ctx.Tdfm = Tab;
		ctx.Tdfm.SetXY( 0, 0 );
		ctx.Tdfm = R * ctx.Tdfm;
	}

	ctx.Tdfm.TPrint( flog, "Tdfm=" );

	if( _arg.FLD ) {
		ctx.FLD = _arg.FLD;
		fprintf(flog, "FLD=%c\n", _arg.FLD );
	}

	if( ctx.FLD == 'X' ) {
		ctx.FLD = (GBL.A.z == GBL.B.z ? 'N' : 'Y');
		fprintf(flog, "FLD=%c (was X)\n", ctx.FLD );
	}

	if( _arg.MODE ) {
		ctx.MODE = _arg.MODE;
		fprintf(flog, "MODE=%c\n", _arg.MODE );
	}

	if( ctx.MODE == 'Z' ) {
		ctx.MODE = 'C';
		arg.CTR = 0.0;
		fprintf(flog, "MODE=C (was Z)\n" );
	}
	else if( ctx.MODE == 'M' ) {
		ctx.MODE = 'N';
		arg.CTR = 0.0;
		fprintf(flog, "MODE=N (was M)\n" );
	}

	if( arg.CTR != 999.0 )
		fprintf(flog, "CTR=%g\n", arg.CTR );

	fprintf(flog, "\n" );

	return true;
}


