// CommonC.cpp: Common classes and functions needed by most files.
//
//////////////////////////////////////////////////////////////////////

#include "CommonC.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

int this_node = -1;
int n_nodes = -1;

/******************* exported variables **********************/
/** buffer for error messages during the integration process. NULL if no errors occured. */
//char *error_msg;
int n_error_msg = 0;

void warning (char *warning_msg)
{
	static int n_warn=0;

	cout << "WARNING: " << warning_msg << endl;

	n_warn++;
	if (n_warn >= MAX_W)
	{
		throw ComException("maximum number of warnings exceeded.");
	}
}

/*  WCLOCK: Get time in secs  */
double wclock()
{
	double  wclock;
	
	wclock  = (double) clock()/CLOCKS_PER_SEC;
	return (wclock);
}