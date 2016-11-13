

#include	"CGBL_dmesh.h"
#include	"dmesh.h"
#include	"ApproximateMatch.h"
#include	"RegionToRegionMap.h"

#include	"Disk.h"
#include	"Inspect.h"
#include	"LinEqu.h"

#include	<stdlib.h>
#include	<string.h>

#define CALLOC(type,num) (type *) calloc((unsigned)(num),(sizeof (type)))

#define REALLOC(ptr,type,num) \
        (type *) realloc((type *)(ptr), (unsigned)((num)*sizeof(type)))

/* --------------------------------------------------------------- */
/* Macros -------------------------------------------------------- */
/* --------------------------------------------------------------- */

#define	FITAFF		0
#define	FITHMG		0
#define	FITDRAW		0
#define	FITTAB		0
#define	FITPOINTS	(FITAFF || FITHMG)

/* --------------------------------------------------------------- */
/* Types --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CStatus {
public:
	int	argn,
		brgn,
		thmok,
		ntri;
public:
	CStatus( int argn, int brgn )
	: argn(argn), brgn(brgn), thmok(0), ntri(0) {};
};


#if FITPOINTS

class PPair {
public:
	Point	A, B;
public:
	PPair( const Point& A, const Point& B ) : A(A), B(B) {};
};

#endif

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* Statics ------------------------------------------------------- */
/* --------------------------------------------------------------- */






/* --------------------------------------------------------------- */
/* RoughMatch ---------------------------------------------------- */
/* --------------------------------------------------------------- */

static bool RoughMatch(
	vector<TAffine>		&guesses,
	const PixPair		&px,
	CCropMask			&CM,
	const ConnRegion	&acr,
	const ConnRegion	&bcr,
	FILE*				flog )
{
	if( guesses.size() > 0 )
		return true;

	if( GBL.ctx.FLD == 'N'
		&& !GBL.mch.PXRESMSK
		&& !CM.IsFile( GBL.idb ) ) {

		// Call NoCR at most once. Possible states
		// are {0=never called, 1=failed, 2=success}.

		static vector<TAffine>	T;
		static int				state = 0;

		int	calledthistime = false;

		if( !state ) {
			state = 1 + ApproximateMatch_NoCR( T, px, flog );
			calledthistime = true;
		}

		if( state == 2 ) {

			if( !calledthistime ) {
				fprintf( flog, "\n---- Thumbnail matching ----\n" );
				T[0].TPrint( flog, "Reuse Approx: Best transform " );
			}

			guesses.push_back( T[0] );
			return true;
		}

		if( !calledthistime ) {
			fprintf( flog, "\n---- Thumbnail matching ----\n" );
			fprintf( flog, "FAIL: Approx: Case already failed.\n" );
		}

		return false;
	}
	else
		return ApproximateMatch( guesses, px, acr, bcr, flog );
}

/* --------------------------------------------------------------- */
/* UpscaleCoords ------------------------------------------------- */
/* --------------------------------------------------------------- */

static void UpscaleCoords( ffmap &maps, int scale )
{
	int	nT;

	if( scale > 1 && (nT = maps.transforms.size()) ) {

		for( int i = 0; i < nT; ++i ) {

			maps.transforms[i].MulXY( scale );

			maps.centers[i].x *= scale;
			maps.centers[i].y *= scale;
		}
	}
}


/* --------------------------------------------------------------- */
/* FitAffine ----------------------------------------------------- */
/* --------------------------------------------------------------- */

#if FITAFF

static void FitAffine(
	const PixPair		&px,
	const vector<PPair>	&pair,
	int					argn,
	int					brgn,
	FILE				*flog )
{
	int	np = pair.size();

	if( np < 3 ) {
		fprintf( flog,
		"Pipe: Too few points to fit affine [%d].\n", np );
		return;
	}

// Create system of normal equations

	double	RHS[6];
	double	LHS[6*6];
	int		i1[3] = { 0, 1, 2 },
			i2[3] = { 3, 4, 5 };

	memset( RHS, 0, 6   * sizeof(double) );
	memset( LHS, 0, 6*6 * sizeof(double) );

	for( int i = 0; i < np; ++i ) {

		const Point&	A = pair[i].A;
		const Point&	B = pair[i].B;

		double	v[3] = { A.x, A.y, 1.0 };

		AddConstraint_Quick( LHS, RHS, 6, 3, i1, v, B.x );
		AddConstraint_Quick( LHS, RHS, 6, 3, i2, v, B.y );
	}

// Solve

	printf( "Pipe: Aff solver returns: %d\n",
		Solve_Quick( LHS, RHS, 6 ) );

	TAffine	T( &RHS[0] );

// Report

	T.TPrint( flog, "Pipe: FitAffine: " );

#if FITTAB
	{
		CMutex		M;
		char		name[256];
		const char	*sud;

		// set file type and layer

		if( GBL.A.z < GBL.B.z )
			sud = "up";
		else if( GBL.A.z == GBL.B.z )
			sud = "same";
		else
			sud = "down";

		sprintf( name, "%s_%d_A", sud, GBL.A.z );

		if( M.Get( name ) ) {

			sprintf( name, "aff.%s", sud );
			FILE *f = fopen( name, "a" );

			if( f ) {

				// Write entry

				fprintf( f,
				"AFFINE"
				" %d.%d-%d %d.%d-%d"
				" %f %f %f %f %f %f\n",
				GBL.A.z, GBL.A.id, argn,
				GBL.B.z, GBL.B.id, brgn,
				RHS[0], RHS[1], RHS[2],
				RHS[3], RHS[4], RHS[5] );

				fflush( f );
				fclose( f );
			}

			M.Release();
		}
	}
#endif

// RMS error

	double	E = 0;

	for( int i = 0; i < np; ++i ) {

		Point	a = pair[i].A;

		T.Transform( a );

		double	err = a.DistSqr( pair[i].B );

		E += err;
	}

	E = sqrt( E / np );

	fprintf( flog, "Pipe: FitAffineRMSerr: %g\n", E );

// Paint

#if FITDRAW
	YellowView( px, T );
#endif
}

#endif

/* --------------------------------------------------------------- */
/* FitHmgphy ----------------------------------------------------- */
/* --------------------------------------------------------------- */

#if FITHMG

static void FitHmgphy(
	const PixPair		&px,
	const vector<PPair>	&pair,
	int					argn,
	int					brgn,
	FILE				*flog )
{
	int	np = pair.size();

	if( np < 4 ) {
		fprintf( flog,
		"Pipe: Too few points to fit homography [%d].\n", np );
		return;
	}

// Create system of normal equations

	double	RHS[8];
	double	LHS[8*8];
	int		i1[5] = { 0, 1, 2, 6, 7 },
			i2[5] = { 3, 4, 5, 6, 7 };

	for( int i = 0; i < np; ++i ) {

		const Point&	A = pair[i].A;
		const Point&	B = pair[i].B;

		double	v[5] = { A.x, A.y, 1.0, -A.x*B.x, -A.y*B.x };

		AddConstraint_Quick( LHS, RHS, 8, 5, i1, v, B.x );

		v[3] = -A.x*B.y;
		v[4] = -A.y*B.y;

		AddConstraint_Quick( LHS, RHS, 8, 5, i2, v, B.y );
	}

// Solve

	printf( "Pipe: Hmg solver returns: %d\n",
		Solve_Quick( LHS, RHS, 8 ) );

	THmgphy	T( &RHS[0] );

// Report

	T.TPrint( flog, "Pipe: FitHmgphy: " );

#if FITTAB
	{
		CMutex		M;
		char		name[256];
		const char	*sud;

		// set file type and layer

		if( GBL.A.z < GBL.B.z )
			sud = "up";
		else if( GBL.A.z == GBL.B.z )
			sud = "same";
		else
			sud = "down";

		sprintf( name, "%s_%d_H", sud, GBL.A.z );

		if( M.Get( name ) ) {

			sprintf( name, "hmg.%s", sud );
			FILE *f = fopen( name, "a" );

			if( f ) {

				// Write entry

				fprintf( f,
				"HMGPHY"
				" %d.%d-%d %d.%d-%d"
				" %f %f %f %f %f %f %.12g %.12g\n",
				GBL.A.z, GBL.A.id, argn,
				GBL.B.z, GBL.B.id, brgn,
				RHS[0], RHS[1], RHS[2],
				RHS[3], RHS[4], RHS[5],
				RHS[6], RHS[7] );

				fflush( f );
				fclose( f );
			}

			M.Release();
		}
	}
#endif

// RMS error

	double	E = 0;

	for( int i = 0; i < np; ++i ) {

		Point	a = pair[i].A;

		T.Transform( a );

		double	err = a.DistSqr( pair[i].B );

		E += err;
	}

	E = sqrt( E / np );

	fprintf( flog, "Pipe: FitAHmgphyRMSerr: %g\n", E );

// Paint

#if FITDRAW
	YellowView( px, T );
#endif
}

#endif

/* --------------------------------------------------------------- */
/* StoreEntry ---------------------------------------------------- */
/* --------------------------------------------------------------- */
static void StoreEntry(float **pGroupId, int **pId, float **p_x, float **p_y, 
                       float **qGroupId, int **qId, float **q_x, float **q_y,
                       int **w, int *num_entries, int *max_num_entries,
                       float GBL_A_z, int GBL_A_id, int ma, float pa_x, float pa_y,
                       float GBL_B_z, int GBL_B_id, int mb, float pb_x, float pb_y )
{
    if (*num_entries >= *max_num_entries)
    {
       *max_num_entries = *max_num_entries + 10;
       *pGroupId = REALLOC(*pGroupId, float, *max_num_entries);
       *qGroupId = REALLOC(*qGroupId, float, *max_num_entries);
       *pId      = REALLOC(*pId,      int  , *max_num_entries);
       *qId      = REALLOC(*qId,      int  , *max_num_entries);
       *p_x      = REALLOC(*p_x,      float, *max_num_entries);
       *p_y      = REALLOC(*p_y,      float, *max_num_entries);
       *q_x      = REALLOC(*q_x,      float, *max_num_entries);
       *q_y      = REALLOC(*q_y,      float, *max_num_entries);
       *w        = REALLOC(*w  ,      int  , *max_num_entries);
    }
   (*pGroupId)[*num_entries] = GBL_A_z;
   (*qGroupId)[*num_entries] = GBL_B_z;
   (*pId)[     *num_entries] = GBL_A_id;
   (*qId)[     *num_entries] = GBL_B_id;
   (*p_x)[     *num_entries] = pa_x;
   (*p_y)[     *num_entries] = pa_y;
   (*q_x)[     *num_entries] = pb_x;
   (*q_y)[     *num_entries] = pb_y;
   (*w)[       *num_entries] = 1;       
    (*num_entries)++;               
}

/* --------------------------------------------------------------- */
/* WritePtsJSON -------------------------------------------------- */
/* --------------------------------------------------------------- */
static void WritePtsJSON(float *pGroupId, int *pId, float *px, float *py,
                         float *qGroupId, int *qId, float *qx, float *qy,
                         int *w, int num_entries, FILE *f)
{
    int n = 0;   // number of 'good' entries
    for (int i=0; i<num_entries; i++)
        if (px[i] > 0)
            n++;
    fprintf(f, "{\n");
    fprintf(f, "    \"matches\": {\n");
    fprintf(f, "        \"p\": ");
    fprintf(f, "[[");
    if (n > 0) {
        for (int i=0; i<n; i++) {
            if (i < n-1)
                fprintf(f, "%.4f,", px[i]);
            else
                fprintf(f, "%.4f],[", px[i]);
        }
        for (int i=0; i<n; i++) {
            if (i < n-1)
                fprintf(f, "%.4f,", py[i]);
            else
                fprintf(f, "%.4f]],\n", py[i]);
        }
    } else {
        fprintf(f, "]],\n");
    }
    fprintf(f, "        \"q\": ");
    fprintf(f, "[[");
    if (n > 0) {
        for (int i=0; i<n; i++) {
            if (i < n-1)
                fprintf(f, "%.4f,", qx[i]);
            else
                fprintf(f, "%.4f],[", qx[i]);
        }
        for (int i=0; i<n; i++) {
            if (i < n-1)
                fprintf(f, "%.4f,", qy[i]);
            else
                fprintf(f, "%.4f]],\n", qy[i]);
        }
    } else {
        fprintf(f, "]],\n");
    }
    fprintf(f, "        \"w\": ");
    fprintf(f, "[");
    if (n > 0) {
        for (int i=0; i<n; i++) {
            if (i < n-1)
                fprintf(f, "%d,", 1);
            else
                fprintf(f, "%d]\n", 1);
        }
    } else {
        fprintf(f, "]\n");
    }
    fprintf(f, "    }\n");
    fprintf(f, "}\n");
}

/* --------------------------------------------------------------- */
/* WritePOINTEntries --------------------------------------------- */
/* --------------------------------------------------------------- */

// CPOINT2 entries are reported at full size.
//
static void WritePOINTEntries(
	const PixPair	      &px,
	const vector<CStatus> &vstat,
	const ffmap    	      &maps,
	const uint8           *fold_mask_a,
	const uint8           *fold_mask_b,
	int		       wf,
	int		       hf,
        const char            *pts_file,
	FILE                  *flog )
{
	fprintf( flog,
	"\n---- Tabulating point matches ----\n" );

	const char *sud;
	CMutex	    M;
	char	    name[256];
        int         num_entries, max_num_entries;
        float      *pGroupId, *qGroupId;
        int        *pId, *qId, *w, res;
        float      *p_x, *p_y, *q_x, *q_y;
        FILE       *f;

        // set pts file type and layer
        if( GBL.A.z < GBL.B.z )
                sud = "up";
        else if( GBL.A.z == GBL.B.z )
                sud = "same";
        else
                sud = "down";

        sprintf( name, "%s_%d", sud, GBL.A.z );

        if (pts_file != NULL) {
            num_entries = 0;
            max_num_entries = 20;
            pGroupId = CALLOC(float, max_num_entries);
            qGroupId = CALLOC(float, max_num_entries);
            pId      = CALLOC(int,   max_num_entries);
            qId      = CALLOC(int,   max_num_entries);
            p_x      = CALLOC(float, max_num_entries);
            p_y      = CALLOC(float, max_num_entries);
            q_x      = CALLOC(float, max_num_entries);
            q_y      = CALLOC(float, max_num_entries);
            w        = CALLOC(int,   max_num_entries);
            sprintf( name, "pts.%s", sud );
        }

	if ( M.Get( name ) ) 
	{

                if ( pts_file == NULL) {
                    f = fopen( name, "a" );
                } else {
                    if ( strstr(pts_file, "stdout")) 
                        f = stdout;
                    else
		        f = fopen( pts_file, "w" );
                }

		if( f ) {

                        // Write a "regular" pts file
                        //
			int	nstat = vstat.size(), tr0 = 0;

			for( int istat = 0; istat < nstat;
				tr0 += vstat[istat++].ntri ) 
                        {

				int	trlim = vstat[istat].ntri;

				if( !trlim )
					continue;

#if FITPOINTS
				vector<PPair>	pair;
#endif

				trlim += tr0;

				for( int itr = tr0; itr < trlim; ++itr ) {

					Point	pa = maps.centers[itr],
							pb = pa;
					int		ma, mb, ix, iy;

					// Lookup for A-point

					ix = int(pa.x);
					iy = int(pa.y);

					if( ix >= 0 && ix < wf && iy >= 0 && iy < hf ) {

						ma = fold_mask_a[ix + wf*iy];

						if( ma <= 0 ) {

							fprintf( flog,
							"CPOINT2: A-centroid has bad mask value:"
							" mask=%d @ (%d, %d).\n", ma, ix, iy );

							continue;
						}

						if( ma != vstat[istat].argn ) {

							fprintf( flog,
							"CPOINT2: A-centroid not in stat block:"
							" mask=%d, stat.argn=%d.\n",
							ma, vstat[istat].argn );

							continue;
						}
					}
					else {

						fprintf( flog,
						"CPOINT2: A-centroid out of A-image bounds"
						" (%d, %d).\n", ix, iy );

						continue;
					}

					// Lookup for B-point

					maps.transforms[itr].Transform( pb );

					ix = int(pb.x);
					iy = int(pb.y);

					if( ix >= 0 && ix < wf && iy >= 0 && iy < hf ) {

						mb = fold_mask_b[ix + wf*iy];

						if( mb <= 0 ) {

							fprintf( flog,
							"CPOINT2: B-centroid has bad mask value:"
							" mask=%d @ (%d, %d).\n", mb, ix, iy );

							continue;
						}

						if( mb != vstat[istat].brgn ) {

							fprintf( flog,
							"CPOINT2: B-centroid not in stat block:"
							" mask=%d, stat.brgn=%d.\n",
							mb, vstat[istat].brgn );

							continue;
						}
					}
					else {

						fprintf( flog,
						"CPOINT2: B-centroid out of B-image bounds"
						" (%d, %d).\n", ix, iy );

						continue;
					}

                                        if (pts_file != NULL) {
                                            StoreEntry(&pGroupId, &pId, &p_x, &p_y,
                                                       &qGroupId, &qId, &q_x, &q_y,
                                                       &w, &num_entries, &max_num_entries,
                                                       (float)GBL.A.z, (int)GBL.A.id, (int)1, (float)pa.x, (float)pa.y,
                                                       (float)GBL.B.z, (int)GBL.B.id, (int)1, (float)pb.x, (float)pb.y );

                                        } else {
					    // Write a "regular" entry
		    			    fprintf( f,
					    "CPOINT2"
					    " %d.%d-%d %f %f"
					    " %d.%d-%d %f %f\n",
					    GBL.A.z, GBL.A.id, ma, pa.x, pa.y,
					    GBL.B.z, GBL.B.id, mb, pb.x, pb.y );
                                        }
#if FITPOINTS
					// Accumulate point pair

					pair.push_back( PPair( pa, pb ) );
#endif
				}

				// Model transforms from point pairs

#if FITAFF
				FitAffine( px, pair,
					vstat[istat].argn,
					vstat[istat].brgn, flog );
#endif

#if FITHMG
				FitHmgphy( px, pair,
					vstat[istat].argn,
					vstat[istat].brgn, flog );
#endif
			}
                        if (pts_file != NULL) {
                            WritePtsJSON(pGroupId, pId, p_x, p_y,
                                        qGroupId, qId, q_x, q_y,
                                        w, num_entries, f);
                        } 
                        if (strcmp(pts_file, "stdout"))
                        {
			    fflush( f );
			    fclose( f );
                        }
		}
	}
    	M.Release();
}

/* --------------------------------------------------------------- */
/* PipelineDeformableMap ----------------------------------------- */
/* --------------------------------------------------------------- */

// The main routine the pipeline should call.
//
// Ntrans		- returned tform count
// tr_array		- array of these values
// map_mask		- <10=no map, else tform[pix-10]
// px			- a and b image pixels
// fold_mask_a	- 0=fold, 1,2,3...=region #
// fold_mask_b	- 0=fold, 1,2,3...=region #
// flog			- detailed output
//
void PipelineDeformableMap(
	int	       &Ntrans,
	double*	       &tr_array,
	uint16*		map_mask,
	const PixPair  &px,
	const uint8    *fold_mask_a,
	const uint8    *fold_mask_b,
        const char     *pts_file,
	FILE           *flog )
{
	int	wf = px.wf, hf = px.hf;

/* --------------- */
/* Default results */
/* --------------- */

	Ntrans		= 0;
	tr_array	= NULL;

	memset( map_mask, 0, wf * hf * sizeof(uint16) );

/* --------------------------------- */
/* Create the connected region lists */
/* --------------------------------- */

// Note that the connected region lists are always
// at the reduced resolution, if this is used.

	fprintf( flog, "\n---- Connected regions ----\n" );

	vector<ConnRegion>	Acr, Bcr;
	CCropMask			CM;

	if( GBL.ctx.FLD == 'N'
		&& !GBL.mch.PXRESMSK
		&& !CM.IsFile( GBL.idb ) ) {

		fprintf( flog, "Forcing single connected region.\n" );

		ConnRgnForce1( Acr, px.ws, px.hs );
		Bcr = Acr;
	}
	else {

		ConnRgnsFromFoldMask( Acr, fold_mask_a,
			wf, hf, px.scl, uint32(0.9 * GBL.mch.MMA), flog );

		ConnRgnsFromFoldMask( Bcr, fold_mask_b,
			wf, hf, px.scl, uint32(0.9 * GBL.mch.MMA), flog );
	}

/* ----------------------------------------- */
/* Find mappings for each region-region pair */
/* ----------------------------------------- */

	vector<CStatus>	vstat;
	ffmap			maps;  // transforms and centers
	FILE			*ftri	= NULL;

	//ftri = fopen( "Triangles.txt", "w" );

	for( int i = 0; i < Acr.size(); ++i ) {

		for( int j = 0; j < Bcr.size(); ++j ) {

			fprintf( flog, "\n---- Begin A-%d to B-%d ----\n",
			Acr[i].id, Bcr[j].id );

			// start list with user's transform arguments

			CStatus			stat( Acr[i].id, Bcr[j].id );
			vector<TAffine>	guesses = GBL.Tmsh;

			if( RoughMatch( guesses,
					px, CM, Acr[i], Bcr[j], flog ) ) {

				stat.thmok = true;

				// Try to get detailed mesh solution from each
				// guess {user + all returned from RoughMatch}.
				// The first to be successful (diff count > 0)
				// causes break.

				for( int k = 0; k < guesses.size(); ++k ) {

					int	count = maps.transforms.size();

					// Downscale coordinates
					guesses[k].MulXY( 1.0 / px.scl );

					RegionToRegionMap( maps, map_mask,
						px, Acr[i], Bcr[j],
						guesses[k], flog, ftri );

					count = maps.transforms.size() - count;

					stat.ntri = count;

					if( count )
						break;
				}
			}

			vstat.push_back( stat );
		}
	}

	//if( ftri )
	//	fclose( ftri );

/* ------------ */
/* Report total */
/* ------------ */

	Ntrans = maps.transforms.size();

	fprintf( flog,
	"Pipe: Returning %d total transforms.\n", Ntrans );

/* ---------------- */
/* Report by region */
/* ---------------- */

	fprintf( flog,
	"\n---- Summary Region-region results ----\n" );

	fprintf( flog,
	"Pipe: Table Headers {Az t r Bz t r thmok ntri}\n" );

	int	nstat = vstat.size();

	for( int i = 0; i < nstat; ++i ) {

		const CStatus& S = vstat[i];

		fprintf( flog,
		"FOUND: %5d %4d %3d  %5d %4d %3d %3d %3d\n",
		GBL.A.z, GBL.A.id, S.argn,
		GBL.B.z, GBL.B.id, S.brgn,
		S.thmok, S.ntri );
	}

/* ---------------- */
/* Feature matches  */
/* ---------------- */
	if( Ntrans || pts_file != NULL) {
		UpscaleCoords( maps, px.scl );

                WritePOINTEntries( px, vstat, maps,
                	fold_mask_a, fold_mask_b, wf, hf, pts_file, flog );
		tr_array = (double*)malloc( Ntrans * 6 * sizeof(double) );

		for( int i = 0; i < Ntrans; ++i ) {

			maps.transforms[i].ToMatlab();
			maps.transforms[i].CopyOut( tr_array + i*6 );
		}

		fprintf( flog, "\n" );
	}
	else
		memset( map_mask, 0, wf * hf * sizeof(uint16) );
}


