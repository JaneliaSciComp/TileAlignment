

#pragma once


#include	"CThmUtil.h"


/* --------------------------------------------------------------- */
/* Class --------------------------------------------------------- */
/* --------------------------------------------------------------- */

class CGBL_dmesh {

// =====
// Types
// =====

private:
	typedef struct {
		vector<TAffine>	Tdfm,
						Tab;
		double			SCALE,
						XSCALE,
						YSCALE,
						SKEW;
		const char		*ima,	// override idb paths
						*imb;
		int				FLD,
						MODE;
	} PrvDrvArgs;

public:
	typedef struct {
		double		CTR;
		const char	*fma,			// override idb paths
					*fmb;
		bool		Transpose,		// transpose all images
					WithinSection,	// overlap within a section
					Verbose,		// run inspect diagnostics
					SingleFold,		// assign id=1 to all non-fold rgns
					Heatmap;		// run CorrView
	} DriverArgs;

	typedef struct {
		TAffine	Tdfm;
		double	XYCONF,
				NBMXHT,
				HFANGDN,
				HFANGPR,
				RTRSH,
				RIT,
				RFA,
				RFT;
		long	OLAP2D;
		int		FLD,
				MODE,
				THMDEC,
				OLAP1D,
				LIMXY,
				OPT;
	} CntxtDep;

        typedef struct {
                bool        use_idb;
                const char *matchparams_path;
                const char *pts_file;
                const char *comp_file;
                const char *registered_file;
        } Options;
// ============
// Data members
// ============

private:
	PrvDrvArgs		_arg;

public:
	DriverArgs		arg;
	TAffine			Tab;	// start thumbs here
	vector<TAffine>	Tmsh;	// bypass thumbs, start mesh here
	vector<Point>	XYexp;	// command line expected XY
	MatchParams		mch;
	CntxtDep		ctx;
	string			idb;
	PicSpecs		A, B;
        Options                 options;

// =================
// Object management
// =================

public:
	CGBL_dmesh();

// =========
// Interface
// =========

public:
	bool SetCmdLine( int argc, char* argv[] );
        void ShowUsage(int argc, char *argv[]);
};

/* --------------------------------------------------------------- */
/* Globals ------------------------------------------------------- */
/* --------------------------------------------------------------- */

extern CGBL_dmesh	GBL;


