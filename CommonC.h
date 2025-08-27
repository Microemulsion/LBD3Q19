////////////////////////////////////////////////////////////////////////////////
// 09-05-07 Kxsong
// CommonC.h
// Ver. 0.1 09/05/2007
// Common class needed by most files.
////////////////////////////////////////////////////////////////////////////////

#ifndef COMMONC_H_
#define COMMONC_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include<map>
#include<vector>
#include<string>
#include<set>
#include<fstream>
#include<iostream>
#include<math.h>
#include<time.h>

#ifndef CXX
#ifndef ISSP

#endif
#endif

using namespace std;

/*  Define hard limits */
#define  MAX_W           10          /* Max # warnings */
#define  MAX_Y         2048          /* Max Y dimension */
#define  MAX_Z         1024          /* Max Z dimension */

/** Precision for capture of round off errors. */
#define ROUND_ERROR_PREC 1.0e-14


/** index of \ref box_l in \ref #fields */
#define FIELD_BOXL                0  
/** index of \ref DomainDecomposition::cell_grid in  \ref #fields */
#define FIELD_CELLGRID            1
/** index of \ref DomainDecomposition::cell_size in  \ref #fields */
#define FIELD_CELLSIZE            2
/** index of \ref dpd_gamma in  \ref #fields */
#define FIELD_DPD_GAMMA          3
/** index of \ref dpd_r_cut in  \ref #fields */
#define FIELD_DPD_RCUT            4
/** index of \ref langevin_gamma in  \ref #fields */
#define FIELD_LANGEVIN_GAMMA      5
/** index of \ref integ_switch in \ref #fields */
#define FIELD_INTEG_SWITCH        6
/** index of \ref local_box_l in \ref #fields */
#define FIELD_LBOXL               7
/** index of \ref max_cut in \ref #fields */
#define FIELD_MCUT                8
/** index of \ref max_num_cells  in \ref #fields */
#define FIELD_MAXNUMCELLS         9
/** index of \ref max_seen_particle in \ref #fields */
#define FIELD_MAXPART             10
/** index of \ref max_range in \ref #fields */
#define FIELD_MAXRANGE            11
/** index of \ref max_skin in  \ref #fields */
#define FIELD_MAXSKIN             12
/** index of \ref min_num_cells  in \ref #fields */
#define FIELD_MINNUMCELLS         13
/** index of \ref n_layers in  \ref #fields */
#define FIELD_NLAYERS             14
/** index of \ref n_nodes in \ref #fields */
#define FIELD_NNODES              15
/** index of \ref n_total_particles in  \ref #fields */
#define FIELD_NPART               16
/** index of \ref n_particle_types in \ref #fields */
#define FIELD_NPARTTYPE           17
/** index of \ref n_rigidbonds in \ref #fields */
#define FIELD_RIGIDBONDS          18
/** index of \ref node_grid in \ref #fields */
#define FIELD_NODEGRID            19
/** index of \ref nptiso_gamma0 in \ref #fields */
#define FIELD_NPTISO_G0           20
/** index of \ref nptiso_gammav in \ref #fields */
#define FIELD_NPTISO_GV           21
/** index of \ref nptiso_struct::p_ext in \ref #fields */
#define FIELD_NPTISO_PEXT         22      
/** index of \ref nptiso_struct::p_inst in \ref #fields */
#define FIELD_NPTISO_PINST        23
/** index of \ref nptiso_struct::p_inst_av in \ref #fields */
#define FIELD_NPTISO_PINSTAV      24
/** index of \ref nptiso_struct::p_diff in \ref #fields */
#define FIELD_NPTISO_PDIFF        25
/** index of \ref nptiso_struct::piston in \ref #fields */
#define FIELD_NPTISO_PISTON       26
/** index of \ref #periodic in \ref #fields */
#define FIELD_PERIODIC            27
/** index of \ref #skin in \ref #fields */
#define FIELD_SKIN                28
/** index of \ref #temperature in \ref #fields */
#define FIELD_TEMPERATURE         29
/** index of \ref thermo_switch in \ref #fields */
#define FIELD_THERMO_SWITCH       30
/** index of \ref sim_time in  \ref #fields */
#define FIELD_SIMTIME             31
/** index of \ref time_step in \ref #fields */
#define FIELD_TIMESTEP            32
/** index of \ref timing_samples in  \ref #fields */
#define FIELD_TIMINGSAMP          33
/** index of \ref transfer_rate  in \ref #fields */
#define FIELD_TRANSFERRATE        34
/** index of \ref rebuild_verletlist in \ref #fields */
#define FIELD_VERLETFLAG          35
/** index of \ref verlet_reuse in  \ref #fields */
#define FIELD_VERLETREUSE         36
/** index of \ref lattice_switch in \ref #fields */
#define FIELD_LATTICE_SWITCH      37
/** index of \ref dpd_tgamma in \ref #fields */
#define FIELD_DPD_TGAMMA          38
/** index of \ref dpd_tr_cut in \ref #fields */
#define FIELD_DPD_TRCUT          39
/** index of \ref dpd_twf in \ref #fields */
#define FIELD_DPD_TWF          40
/** index of \ref dpd_wf in \ref #fields */
#define FIELD_DPD_WF           41
/** index of \ref address_var in \ref #fields */
#define FIELD_ADRESS           42
/*@}*/


/** \name request codes in asynchronous mode
    Keep this in sync with \ref communication::callbacks and \ref communication::names. */
/*@{*/
/** Action number for \ref mpi_stop. */
#define REQ_TERM      0
/** Action number for \ref mpi_bcast_parameter. */
#define REQ_BCAST_PAR 1
/** Action number for \ref mpi_who_has. */
#define REQ_WHO_HAS   2
/** Action number for \ref mpi_bcast_event. */
#define REQ_EVENT     3
/** Action number for \ref mpi_place_particle. */
#define REQ_PLACE     4

/** Action number for \ref mpi_send_v. */
#define REQ_SET_V     5
/** Action number for \ref mpi_send_f. */
#define REQ_SET_F     6
/** Action number for \ref mpi_send_q. */
#define REQ_SET_Q     7
/** Action number for \ref mpi_send_type. */
#define REQ_SET_TYPE  8
/** Action number for \ref mpi_send_bond. */
#define REQ_SET_BOND  9

/** Action number for \ref mpi_recv_part. */
#define REQ_GET_PART  10
/** Action number for \ref mpi_integrate. */
#define REQ_INTEGRATE 11
/** Action number for \ref mpi_bcast_ia_params. */
#define REQ_BCAST_IA  12
/** Action number for \ref mpi_bcast_n_particle_types. */
#define REQ_BCAST_IA_SIZE  13
/** Action number for \ref mpi_gather_stats. */
#define REQ_GATHER    14

/** Action number for \ref mpi_set_time_step. */
#define REQ_SET_TIME_STEP  15
/** Action number for \ref mpi_get_particles. */
#define REQ_GETPARTS  16
/** Action number for \ref mpi_bcast_coulomb_params. */
#define REQ_BCAST_COULOMB 17
/** Action number for \ref mpi_send_ext. */
#define REQ_SET_EXT 18
/** Action number for \ref mpi_place_particle. */
#define REQ_PLACE_NEW   19

/** Action number for \ref mpi_remove_particle */
#define REQ_REM_PART   20
/** Action number for \ref mpi_bcast_constraint */
#define REQ_BCAST_CONSTR 21
/** Action number for \ref mpi_random_seed */
#define REQ_RANDOM_SEED 22
/** Action number for \ref mpi_random_stat */
#define REQ_RANDOM_STAT 23
/** Action number for \ref mpi_lj_cap_forces. */
#define REQ_BCAST_LFC 24

/** Action number for \ref mpi_tab_cap_forces. */
#define REQ_BCAST_TFC 25
/** Action number for \ref mpi_random_seed */
#define REQ_BIT_RANDOM_SEED 26
/** Action number for \ref mpi_random_stat */
#define REQ_BIT_RANDOM_STAT 27
/** Action number for \ref mpi_get_constraint_force */
#define REQ_GET_CONSFOR 28
/** Action number for \ref mpi_rescale_particles */
#define REQ_RESCALE_PART 29

/** Action number for \ref mpi_bcast_cell_structure */
#define REQ_BCAST_CS  30
/** Action number for \ref mpi_send_quat. */
#define REQ_SET_QUAT  31
/** Action number for \ref mpi_send_omega. */
#define REQ_SET_OMEGA 32
/** Action number for \ref mpi_send_torque. */
#define REQ_SET_TORQUE 33
/** Action number for \ref mpi_send_mol_id. */
#define REQ_SET_MOLID 34

/** Action number for \ref mpi_bcast_nptiso_geom */
#define REQ_BCAST_NPTISO_GEOM 35
/** Action number for \ref mpi_update_mol_ids  */
#define REQ_UPDATE_MOL_IDS 36
/** Action number for \ref mpi_sync_topo_part_info  */
#define REQ_SYNC_TOPO 37
/** Action number for \ref mpi_send_mass. */
#define REQ_SET_MASS  38
/** Action number for \ref mpi_buck_cap_forces. */
#define REQ_BCAST_BFC 39

/** Action number for \ref mpi_gather_runtime_errors  */
#define REQ_GET_ERRS  40
/** Action number for \ref mpi_send_exclusion. */
#define REQ_SET_EXCL  41
/** Action number for \ref mpi_morse_cap_forces. */
#define REQ_BCAST_MFC 42
/** Action number for \ref mpi_bcast_lb_params. */
#define REQ_BCAST_LBPAR 43
/** Action number for \ref mpi_send_dip. */
#define REQ_SET_DIP 44
/** Action number for \ref mpi_send_dipm. */
#define REQ_SET_DIPM 45
/** Action number for \ref mpi_send_fluid. */
#define REQ_SET_FLUID 46
/** Action number for \ref mpi_recv_fluid. */
#define REQ_GET_FLUID 47
/** Action number for \ref mpi_local_stress_tensor*/
#define REQ_GET_LOCAL_STRESS_TENSOR 48
/** Action number for \ref mpi_ljangle_cap_forces. */
#define REQ_BCAST_LAFC 49
/** Action number for \ref mpi_send_isVirtual. */
#define REQ_SET_ISVI 50
/** Total number of action numbers. */
#define REQ_MAXIMUM 51

/*@}*/

/** Tiny angle cutoff for sinus calculations */
#define TINY_SIN_VALUE 1e-10
/** Tiny angle cutoff for cosine calculations */
#define TINY_COS_VALUE 0.9999999999

/*************************************************************/
/** \name Mathematical, physical and chemical constants.     */
/*************************************************************/
/*@{*/
/** Pi. */
#define PI     3.14159265358979323846264338328 
/** Square root of Pi */
#define wupi   1.77245385090551602729816748334 
/** One over square root of Pi. */
#define wupii  0.56418958354775627928034964498 
/** Pi to the power 1/3. */
#define driwu2 1.25992104989487316476721060728 
/*@}*/


#ifdef COMM_DEBUG
#define COMM_TRACE(cmd) { cmd; }
#else
/** Equals { cmd } iff COMM_DEBUG is set. */
#define COMM_TRACE(cmd)
#endif

/*************************************************************/
/* mass helper macro                                         */
/*************************************************************/

//#ifdef MASS
///** macro for easy use of mass. If masses are not switched on, the particle mass is defined to 1,
//    so it should be compiled out in most cases. */
//#define PMASS(pt) (pt).p.mass
//#else
//#define PMASS(pt) 1
//#endif

/** \name Exported Variables */
/*@{*/
/** The number of this node. */
extern int this_node;
/** The total number of nodes. */
extern int n_nodes;
/*@}*/

/** buffer for error messages during the integration process. */
//extern char *error_msg;
extern int n_error_msg;

// A Exception class, which can show the position of the errors.
class ComException {
public:
	ComException(  const string& message  = "NormalTerminate"
		, const string& fileName = ""
		, int   numLine          = 0 )
		:	  p_message(  message )
		, p_fileName( fileName )
		, p_numLine(  numLine )
	{ ; }
	virtual ~ComException() { ; }
	string message() const { return p_message; }
    virtual void display(ostream& out ) const {
		out << p_message ;
		if( p_fileName.size() > 0 ) {
			out << " Error at " << p_fileName;
			if( p_numLine > 0 ) {
				out << " lineNumber:"   << p_numLine;
			}
		}
		out << endl;
	}
protected:
    string p_message;   // Error Message
    string p_fileName;  // file name
    int    p_numLine;   // line number
};

//realloc in C++
template< class T >
T* realloc_p(T *a,int size)
{
	T   *t;
	t  = new T[ size ];
	memcpy( t, a, size*sizeof(T) );
	delete[] a;
	return t;
} 


/*  WARNING: Prints limitted warning messages.*/
void warning (char *warning_msg);

/** Mathematically rounds 'double'-typed x, returning 'double'. */
inline double dround(double x) { return floor(x+0.5); }

/** Calculates the SQuaRe of 'double' x, returning 'double'. */
inline double SQR(double x) { return x*x; }

/** calculates the scalar product of two vectors a nd b */
inline double scalar(double a[3], double b[3]) {
  double d2 = 0.0;
  int i;
  for(i=0;i<3;i++)
    d2 += a[i]*b[i];
  return d2;
}

/************************************************
 * data types
 ************************************************/

/** check for runtime errors on all nodes. This has to be called on all nodes synchronously.
    @return the number of characters in the error messages of all nodes together. */
int check_runtime_errors();

void mpi_issue(int reqcode, int node, int param);

/** Calculates the maximum of 'double'-typed a and b, returning 'double'. */
inline double dmax(double a, double b) { return (a>b) ? a : b; }

/** Calculates the minimum of 'double'-typed a and b, returning 'double'. */
inline double dmin(double a, double b) { return (a<b) ? a : b; }

/** Calculates the maximum of 'int'-typed a and b, returning 'int'. */
inline int imax(int a, int b) { return (a>b) ? a : b; }

/** Calculates the minimum of 'int'-typed a and b, returning 'int'. */
inline int imin(int a, int b) { return (a<b) ? a : b; }

/** Calculates the remainder of a division */
inline double drem_down(double a, double b) { return a - floor(a/b)*b; }

/*************************************************************/
/** \name Three dimensional grid operations                  */
/*************************************************************/
/*@{*/

/** get the linear index from the position (a,b,c) in a 3D grid
 *  of dimensions adim[]. returns linear index.
 *
 * @return        the linear index
 * @param a       x position 
 * @param b       y position 
 * @param c       z position 
 * @param adim    dimensions of the underlying grid  
 */
inline int get_linear_index(int a, int b, int c, int adim[3])
{
  return (a + adim[0]*(b + adim[1]*c));
}


/** get the position a[] from the linear index in a 3D grid
 *  of dimensions adim[].
 *
 * @param i       linear index
 * @param a       x position (return value) 
 * @param b       y position (return value) 
 * @param c       z position (return value) 
 * @param adim    dimensions of the underlying grid  
 */
inline void get_grid_pos(int i, int *a, int *b, int *c, int adim[3])
{
  *a = i % adim[0];
  i /= adim[0];
  *b = i % adim[1];
  i /= adim[1];
  *c = i;
}

/*************************************************************/
/** \name List operations .                                  */
/*************************************************************/
/*@{*/

extern int this_node;





/** calculates the squared length of a vector */
inline double sqrlen(double v[3]) {
  double d2 = 0.0;
  int i;
  for(i=0;i<3;i++)
    d2 += SQR(v[i]);
  return d2;
}

/** Returns the distance between two positions squared and stores the
    distance vector pos1-pos2 in vec.
 *  \param pos1 Position one.
 *  \param pos2 Position two.
 *  \param vec  vecotr pos1-pos2.
 *  \return distance squared
*/
inline double distance2vec(double pos1[3], double pos2[3], double vec[3])
{
  vec[0] = pos1[0]-pos2[0];
  vec[1] = pos1[1]-pos2[1];
  vec[2] = pos1[2]-pos2[2];
  return SQR(vec[0]) + SQR(vec[1]) + SQR(vec[2]);
}

/** calculates the vector product c of two vectors a and b */
inline void vector_product(double a[3], double b[3], double c[3]) {
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
  return ;
}

/** returns the distance between two positions squared.
 *  \param pos1 Position one.
 *  \param pos2 Position two.
*/
inline double distance2(double pos1[3], double pos2[3])
{
  return SQR(pos1[0]-pos2[0]) + SQR(pos1[1]-pos2[1]) + SQR(pos1[2]-pos2[2]);
}

/** returns the distance between two positions squared.
 *  \param pos1 Position one.
 *  \param pos2 Position two.
*/
inline double atanh(double value)
{
  return log((1.0 / value + 1.0) / (1.0 / value - 1.0)) / 2.0;
}

/** Subtracts vector v2 from vector v1 and stores resuld in vector dv */
inline void vecsub(double v1[3], double v2[3], double dv[3])
{
  dv[0] = v1[0] - v2[0];
  dv[1] = v1[1] - v2[1];
  dv[2] = v1[2] - v2[2];
}

/*  WCLOCK: Get time in secs  */
double wclock();

#endif  /* COMMONC_H_ */
