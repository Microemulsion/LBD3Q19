////////////////////////////////////////////////////////////////////////////////
// 07-06-08 Kxsong
// Particles.h
// Ver. 0.1 07/06/2008
// The definitions of particles.
////////////////////////////////////////////////////////////////////////////////
// Particles.h: interface for the Particles class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_PARTICLES_H_)
#define _PARTICLES_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

//#define _CRT_SECURE_NO_WARNINGS
// Disable warning messages 4996.
#pragma warning( disable : 4996 )

/*The cell system (\ref Cell Structure) describes how particles are
distributed on the cells and how particles of different cells
(regardless if they reside on the same or different nodes)
interact with each other. The following cell systems are implemented:
<li> domain decomposition: The simulation box is divided spatially
ino cells (see \ref domain_decomposition.h). This is suitable for
short range interctions.

<li> All cells, real cells as well as ghost cells, are stored in the array 
\ref cells::cells with size \ref n_cells. The size of this array has to be 
changed with \ref realloc_cells. <li> Their are two lists of cell pointers 
to acces particles and ghost particles on a node: \ref local_cells contains
pointers to all cells containing the particles physically residing on that
node. \ref ghost_cells contains pointers to all cells containing
the ghost particles of that node. The size of these lists has to be
changed with \ref realloc_cellplist
<li> An example using the cell pointer lists to access particle
data can be found in the function \ref 
print_local_particle_positions. DO NOT INVENT YOUR OWN WAY!!! */

#include <limits.h>
#include <string.h>
#include <vector>
#include "Mesh.h"

/** CELLS: Default value for the maximal number of cells per node. */
#define CELLS_MAX_NUM_CELLS 216

/** half the number of cell neighbors in 3 Dimensions. */
#define CELLS_MAX_NEIGHBORS 14

/** Flag for exchange_and_sort_particles : Do neighbor exchange. */
#define CELL_NEIGHBOR_EXCHANGE 0

class sPart
{
public:
	sPart(void);
	~sPart(void);

	//ParticleProperties p;
	/* Properties of a particle which are not supposed to
	change during the integration, but have to be known
	for all ghosts. Ghosts are particles which are
	needed in the interaction calculation, but are just copies of
	particles stored on different nodes.*/
	/** unique identifier for the particle. */
	int    identity;
	/* Molecule identifier. */
	int    mol_id;
	/** particle type, used for non bonded interactions. */
	int    type;
	//#ifdef MASS
	/** particle mass */
	double mass;
	//#endif

	//ParticlePosition p;
	/* Positional information on a particle. Information that is
	communicated to calculate interactions with ghost particles. */
	/* periodically folded position. */
	double p[3];

	/** check whether a particle is a ghost or not */
	//bool g_flag;

	// ParticleMomentum v;
	/* Momentum information on a particle. Information not contained in
	communication of ghost particles so far, but a communication would
	be necessary for velocity dependend potentials. */
	/* velocity. */
	double v[3];

	// ParticleForce f;
	/* Force information on a particle. Forces of ghost particles are
	collected and added up to the force of the original particle. */
	/* force. */
	double f[3];

	//  ParticleLocal l;
	/** Information on a particle that is needed only on the
	node the particle belongs to */
	/** position in the last time step befor last Verlet list update. */
	double p_old[3];
	/** index of the simulation box image where the particle really sits. */
	int   i[3];

	//#ifdef EXTERNAL_FORCES
	/** flag whether to fix a particle in space.
	Values:
	<ul> <li> 0 no external influence
	<li> 1 apply external force \ref ParticleLocal::ext_force
	<li> 2,3,4 fix particle coordinate 0,1,2
	</ul>*/
	int ext_flag;
	/** External force, apply if \ref ParticleLocal::ext_flag == 1. */
	double ext_force[3];
	//#endif


	//#ifdef LB  ParticleLatticeCoupling lc;
	//#endif
	//#ifdef LB
	/** Data related to the Lattice Boltzmann hydrodynamic coupling */
	/** fluctuating part of the coupling force */
	double f_random[3];

	/** bonded interactions list. The format is pretty simple:
	Just the bond type, and then the particle ids. The number of particle ids can be determined
	easily from the bonded_ia_params entry for the type. */
	vector<int> bl;//IntList* bl;
};

typedef struct {
	int neighborCell;
	/** Verlet list for non bonded interactions of a cell with a neighbor cell. */
	vector< pair< sPart* , sPart* > > vList;
} VerletList;

/** A cell is a \ref ParticleList representing a particle group with
respect to the integration algorithm.*/
class Cell
{
public:
	Cell(void);
	~Cell(void);
	vector<sPart*>     p_objects;

	/** number of interacting neighbor cells . 
	A word about the interacting neighbor cells:
	In a 3D lattice each cell has 27 neighbors (including
	itself!). Since we deal with pair forces, it is sufficient to
	calculate only half of the interactions (Newtons law: actio =
	reactio). For each cell 13+1=14 neighbors. This has only to be
	done for the inner cells. 
	Caution: This implementation needs double sided ghost
	communication! For single sided ghost communication one would
	need some ghost-ghost cell interaction as well, which we do not
	need!
	It follows: inner cells: n_neighbors = 14
	ghost cells:             n_neighbors = 0
	*/
	int n_neighbors;
	VerletList cell_inter[14];
	/** bool to mark the cell containing ghosts */
	bool p_contGhosts;

	//Ghost Position gp;
	/* periodically folded position. */
	double gPMove[3];

	//Corresponded no ghost cell num;
	int noGCell;

	/** Add a particle pair to a verlet pair list.
	Checks verlet pair list size and reallocates memory if necessary.
	*  \param p1 Pointer to paricle one.
	*  \param p2 Pointer to paricle two.
	*  \param pl Pointer to the verlet pair list.
	*/
	void add_pair(int n, sPart* p1, sPart* p2);

	void clear_pairs(int n);
};

class Particles
{
public:
	Particles(const BaseMesh& pmesh, double max_cut, double skin, double time_step, char* outputPFileName);
	~Particles(void);

	// Get the detail of particles from input file.
	void readInputF(char* inputFileName);

	/*  Particles: output  */
	void outPut (int t);

	/*  Read checkpoint file  */
	void readOutputF(int t);

	//initialize particles .
	void initObjects(int n_particles);

	/** Move a particle to a new position.
	If it does not exist, it is created. the position must
	be on the local node!
	@param part the identity of the particle to move
	@param p    its new position
	@param new  if true, the particle is allocated, else has to exists already */
	void place_particle(int part, double p[3], bool bnew);

	/** fold particle coordinates to primary simulation box.
	\param pos the position...
	\param image_box and the box
	Both pos and image_box are I/O,
	i. e. a previously folded position will be folded correctly. */
	void fold_position(const BaseMesh& p_mesh, double pos[3],int image_box[3]);

	/** unfold coordinates to physical position.
	\param pos the position...
	\param image_box and the box

	Both pos and image_box are I/O, i.e. image_box will be (0,0,0)
	afterwards.
	*/
	void unfold_position(double pos[3],int image_box[3]);

	/** Update the entries in \ref local_particles for all particles in the list pl.
	@param pl the list to put in.
	*/
	void update_local_particles(Cell *pl);

	/** Allocate storage for local particles and ghosts. This version
	does \em not care for the bond information to be freed if necessary.
	\param plist the list on which to operate
	\param size the size to provide at least. It is rounded
	up to multiples of \ref PART_INCREMENT.
	\return true iff particle adresses have changed */
	int realloc_particlelist(Cell *plist, int size);

	/** implements \ref CellStructure::position_to_cell. */
	Cell* position_to_cell(double pos[3]);

	/************************************************
	* exported variables
	************************************************/

	/** total number of particles on all nodes. */
	int  n_total_particles;

	/** id->particle mapping on all nodes. This is used to find partners
	of bonded interactions. */
	vector<sPart*> local_particles;

	/************************************************************/
	/** \name Exported Variables */
	/************************************************************/
	/*@{*/

	/** list of all cells. */
	/** Initialize the cell structure on program start with the default
	cell structure of type domain decomposition. */
	void cells_pre_init();

	/** Calculate maximal interaction range. 
	Uses \ref calc_maximal_cutoff.
	\ref max_range  = \ref max_cut + \ref #skin; */
	void integrate_vv_recalc_maxrange(double max_cut, double skin);

	/** Calculate cell grid dimensions, cell sizes and number of cells.
	*  Calculates the cell grid, based on \ref local_box_l and \ref
	*  max_range. If the number of cells is larger than \ref
	*  max_num_cells, it increases max_range until the number of cells is
	*  smaller or equal \ref max_num_cells. It sets: \ref
	*  DomainDecomposition::cell_grid, \ref
	*  DomainDecomposition::ghost_cell_grid, \ref
	*  DomainDecomposition::cell_size, \ref
	*  DomainDecomposition::inv_cell_size, and \ref n_cells.
	*/
	void dd_create_cell_grid(double max_cut);

	/** get the minimal distance vector of two vectors in the current bc.
	@param a the vector to subtract from
	@param b the vector to subtract
	@param res where to store the result
	*/
	void get_mi_vector(double res[3], double a[3], double b[3]);

	/** update ghost information. If \ref rebuild_verletlist == 1, for some cell structures also a
	resorting of the particles takes place. */
	void cells_update_ghosts();

	/** Init cell interactions for cell system domain decomposition.
	* initializes the interacting neighbor cell list of a cell The
	* created list of interacting neighbor cells is used by the verlet
	* algorithm (see verlet.c) to build the verlet lists.
	*/
	void dd_init_cell_interactions();

	/** Add velocity to all particles. */
	void addVelocity(double* vel);
	
	/** Get total moment. */
	void predict_momentum_particles(double* momentum);

	/** Get kinetic energy. */
	double calc_kinetic_energy();

	/** function that rescales all velocities on one node according to a
    new time step. */
	void rescale_velocities(double scale);

	/*  Particles: output  properties*/
	void outProp(int t, FILE *file_ptr);

	/*  Read properties file  */
	void inProp(int t, FILE	*file_ptr);

	/** update ghost cells. */
	void update_ghosts();

	/** list of all cells. */
	vector<Cell*> local_cells;

	/** list of all cells containing ghosts */
	vector<Cell*> ghost_cells;

	/** list of all cells. */
	vector<Cell*> cells;

	/** size of \ref cells::cells */
	int n_cells;

	/** If non-zero, the verlet list has to be rebuilt. */
	bool rebuild_verletlist;

	/** flag for using Verlet List */
	bool use_vList;

	double max_skin;

private:
	const BaseMesh&       p_mesh;

	/** Used by \ref mpi_place_particle, should not be used elsewhere.
	Called if on a different node a new particle was added.
	@param part the identity of the particle added */
	sPart* add_one_particle();

	/** The number of nodes in each spatial dimension. */
	int node_grid[3];
	/** position of node in node grid */
	int node_pos[3];

	/** Dimensions of the box a single node is responsible for. */ 
	double local_box_l[3];

	/** linked cell grid in nodes spatial domain. */
	int cell_grid[3];

	/** linked cell grid with ghost frame. */
	int ghost_cell_grid[3];

	/** cell size. 
	Def: \verbatim cell_grid[i] = (int)(local_box_l[i]/max_range); \endverbatim */
	double cell_size[3];

	/** Left (bottom, front) corner of this nodes local box. */ 
	double my_left[3];
	/** Right (top, back) corner of this nodes local box. */ 
	double my_right[3];

	/** Maximal interaction range (max_cut + skin). */
	double max_range;
	/** Square of \ref max_range. It's initial value is -1.0 which is
	used to determine wether max_range/max_range2 has been set
	properly by \ref integrate_vv_recalc_maxrange or not. */
	double max_range2;


	/** Maximal number of cells per node. In order to avoid memory
	*  problems due to the cell grid one has to specify the maximal
	*  number of \ref cells::cells . The corresponding callback function
	*  is \ref max_num_cells_callback. If the number of cells \ref
	*  n_cells, is larger than max_num_cells the cell grid is
	*  reduced. max_num_cells has to be larger than 27, e.g one inner
	*  cell.  max_num_cells is initialized with the default value
	*  specified in \ref config.h: \ref CELLS_MAX_NUM_CELLS.
	*/
	int max_num_cells;
	int min_num_cells;

	double p_timestep;

	char* p_outputPFileName;

	/** Reallocate the list of all cells (\ref cells::cells). */
	void realloc_cells(vector<Cell*>& icells, int size);

	/** sort the particles into the cells and initialize the ghost particle structures.
	@param global_flag if this is CELLS_GLOBAL_EXCHANGE, particle positions can have changed
	arbitrarly, otherwise the change should have been smaller then skin.  */
	void cells_resort_particles(int global_flag);

	/** Go through \ref ghost_cells and remove the ghost entries from \ref
	local_particles. Part of \ref dd_exchange_and_sort_particles.*/
	void invalidate_ghosts();

	/** Just resort the particles. Used during integration. The particles
	are stored in the cell structure.

	@param global_flag Use DD_GLOBAL_EXCHANGE for global exchange and
	DD_NEIGHBOR_EXCHANGE for neighbor exchange (recommended for use within
	Molecular dynamics, or any other integration scheme using only local
	particle moves) 
	*/
	void dd_exchange_and_sort_particles(int global_flag);

	/** fold a coordinate to primary simulation box.
	\param pos         the position...
	\param image_box   and the box
	\param dir         the coordinate to fold: dir = 0,1,2 for x, y and z coordinate.

	Both pos and image_box are I/O,
	i. e. a previously folded position will be folded correctly.
	*/
	void fold_coordinate(double pos[3], int image_box[3], int dir);

	/** Remove a particle from one particle List and append it to another.
	Refill the sourceList with last particle and update its entry in
	local_particles. Reallocates particles if necessary.  This
	procedure cares for \ref local_particles.
	\param destList   List where the particle is appended.
	\param sourceList List where the particle will be removed.
	\param ind        Index of the particle in the sourceList.
	\return Pointer to new location of the particle.
	*/
	void move_indexed_particle(Cell *destList, Cell *sourceList, int ind);
};

#endif // !defined(_PARTICLES_H_)
