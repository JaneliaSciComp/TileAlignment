

#include	"Cmdline.h"

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include        <curl/curl.h>
#include        <sstream>

#include        <iostream>


/* --------------------------------------------------------------- */
/* IsArg --------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Return true if matched command line parameter.
//
// Example usage:
//
//	for( int i = 1; i < argc; ++i ) {
//		if( argv[i][0] != '-' )
//			noa.push_back( argv[i] );
//		else if( IsArg( "-nf", argv[i] ) )
//			NoFolds = true;
//	}
//
bool IsArg( const char *pat, const char *argv )
{
	return !strcmp( argv, pat );
}

/* --------------------------------------------------------------- */
/* GetArg -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Read argument from command line.
//
// Example usage:
//
//	for( int i = 1; i < argc; ++i ) {
//		if( argv[i][0] != '-' )
//			noa.push_back( argv[i] );
//		else if( GetArg( &ApproxScale, "-SCALE=%lf", argv[i] ) )
//			;
//		else if( GetArg( &Order, "-ORDER=%d", argv[i] ) )
//			;
//	}
//
bool GetArg( void *v, const char *pat, const char *argv )
{
	return 1 == sscanf( argv, pat, v );
}

/* --------------------------------------------------------------- */
/* GetArgStr -------------------------------------------------------- */
/* --------------------------------------------------------------- */

// Point at string argument on command line.
//
// Example usage:
//
//	char	*dirptr;
//
//	for( int i = 1; i < argc; ++i ) {
//		if( argv[i][0] != '-' )
//			noa.push_back( argv[i] );
//		else if( GetArgStr( dirptr, "-d=", argv[i] ) )
//			;
//	}
//
bool GetArgStr( const char* &s, const char *pat, char *argv )
{
	int		len = strlen( pat );
	bool	ok = false;

	if( !strncmp( argv, pat, len ) ) {

		s  = argv + len;
		ok = true;
	}

	return ok;
}

/* --------------------------------------------------------------- */
/* GetArgList ---------------------------------------------------- */
/* --------------------------------------------------------------- */

// Read integer argument list from command line.
//
// Example usage: ... -List=2,5,7 ...
//
//	vector<int>	I;
//
//	for( int i = 1; i < argc; ++i ) {
//		if( argv[i][0] != '-' )
//			noa.push_back( argv[i] );
//		else if( GetArgList( I, "-List=", argv[i] ) )
//			;
//	}
//
bool GetArgList( vector<int> &v, const char *pat, char *argv )
{
	int		len = strlen( pat );
	bool	ok = false;

	if( !strncmp( argv, pat, len ) ) {

		char	*s = strtok( argv + len, ":;, " );

		v.clear();

		while( s ) {
			v.push_back( atoi( s ) );
			s = strtok( NULL, ":;, " );
		}

		ok = true;
	}

	return ok;
}

// Read double argument list from command line.
//
// Example usage: ... -List=2.7,5,1.8e7 ...
//
//	vector<double>	D;
//
//	for( int i = 1; i < argc; ++i ) {
//		if( argv[i][0] != '-' )
//			noa.push_back( argv[i] );
//		else if( GetArgList( D, "-List=", argv[i] ) )
//			;
//	}
//
bool GetArgList( vector<double> &v, const char *pat, char *argv )
{
	int		len = strlen( pat );
	bool	ok = false;

	if( strstr(argv, "http:") == NULL && !strncmp( argv, pat, len ) ) {

		char	*s = strtok( argv + len, ":;, " );

		v.clear();

		while( s ) {
			v.push_back( atof( s ) );
			s = strtok( NULL, ":;, " );
		}

		ok = true;
	}

	return ok;
}

void Tokenize(const string& str, vector<string>& tokens,
              const string& delimiters = "\n")
{
     // Skip delimiters at beginning.
     string::size_type lastPos = str.find_first_not_of(delimiters, 0);
     // Find first "non-delimiter".
     string::size_type pos     = str.find_first_of(delimiters, lastPos);

     while (string::npos != pos || string::npos != lastPos)
     {
         // Found a token, add it to the vector.
         tokens.push_back(str.substr(lastPos, pos - lastPos));
         // Skip delimiters.  Note the "not_of"
         lastPos = str.find_first_not_of(delimiters, pos);
         // Find next "non-delimiter"
         pos = str.find_first_of(delimiters, lastPos);
     }
}

static size_t WriteCallback(void *contents, size_t size, size_t nmemb, void *userp)
{
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}


bool GetArgListFromURL( vector<double> &v, const char *pat, char *argv )
{
    int     len = strlen( pat );
    bool    ok = false;

    if ( strstr(argv, "http:") != NULL && !strncmp( argv, pat, len ) ) 
    {

        char *url = argv + len;
        static std::string readBuffer;

        if (strstr(url, "http:") == NULL)
            return false;

        CURL *easy_handle;
        CURLcode res;

        easy_handle = curl_easy_init();

//      curl_easy_setopt(easy_handle, CURLOPT_VERBOSE,    1L);
        curl_easy_setopt(easy_handle, CURLOPT_URL, url);
        readBuffer.clear();
        curl_easy_setopt(easy_handle, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(easy_handle, CURLOPT_WRITEDATA, &readBuffer);
        res = curl_easy_perform(easy_handle);
        curl_easy_cleanup(easy_handle);

        vector<string> tokens;

        int read_images_from_URLs = 0;
        double coeffs[] = {0., 0., 0., 0., 0, 0.};

        Tokenize(std::string(readBuffer), tokens, "\n");
        for (int i=0; i<tokens.size(); i++)
        {
            vector<string> tokens2;
            Tokenize(tokens[i], tokens2, " ");

            if ( tokens2[0].find("imageRow") != std::string::npos ) {
                const char *my_str    = tokens2[2].c_str();
                char *pEnd;
                double value = strtod(my_str, &pEnd);
                coeffs[3] = value;
            }

            if ( tokens2[0].find("imageCol") != std::string::npos ) {
                const char *my_str    = tokens2[2].c_str();
                char *pEnd;
                double value = strtod(my_str, &pEnd);
                coeffs[0] = value;
            }

            if ( tokens2[0].find("minX") != std::string::npos ) {
                const char *my_str    = tokens2[2].c_str();
                char *pEnd;
                double value = strtod(my_str, &pEnd);
                coeffs[1] = value;
            }

            if ( tokens2[0].find("minY") != std::string::npos ) {
                const char *my_str    = tokens2[2].c_str();
                char *pEnd;
                double value = strtod(my_str, &pEnd);
                coeffs[4] = value;
            }

            if ( tokens2[0].find("stageX") != std::string::npos ) {
                const char *my_str    = tokens2[2].c_str();
                char *pEnd;
                double value = strtod(my_str, &pEnd);
                coeffs[2] = value;
            }

            if ( tokens2[0].find("stageY") != std::string::npos ) {
                const char *my_str    = tokens2[2].c_str();
                char *pEnd;
                double value = strtod(my_str, &pEnd);
                coeffs[5] = value;
            }
        }
        v.assign (coeffs,coeffs+6);
        readBuffer.clear();
        ok = true;
    }
    return ok;
}

