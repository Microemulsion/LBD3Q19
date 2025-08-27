////////////////////////////////////////////////////////////////////////////////
// Jan 29, 2009 Kxsong
// LBD3Q19.h
// Ver. 0.1 1/29/2009
// The lattice Boltzmann D3Q19 classes.
////////////////////////////////////////////////////////////////////////////////

#ifndef LBD3Q19_H_
#define LBD3Q19_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>
#include "Mesh.h"
#include "random.h"
#include "Particles.h"

class LBD3Q19
{
public:
	LBD3Q19(const BaseMesh& pmesh, const double ptau, const double ptauv, const double pagrid, double pviscosity, double bulk_viscosity, double pfriction, double ptotal_rho, double ptemperature, double ptime_step, char* OFFileName);
	~LBD3Q19(void);

	/** The Lattice Boltzmann streaming step.
	* The populations are moved to the neighbouring lattice sites
	* according to the velocity sublattice. This can be done in two ways:
	* First, one can use a temporary field to store the updated configuration.
	* Second, one can order the updates such that only populations are
	* overwritten which have already been propagated. The halo region
	* serves as a buffer. This requires two sweeps through the lattice,
	* one bottom up and one bottom down. One has to be careful if the
	* velocities are upgoing or downgoing. This can be a real bugfest!
	*/
	void stream();

	/** Apply external forces to the fluid.
	*
	* Eq. (28) Ladd and Verberg, J. Stat. Phys. 104(5/6):1191 (2001).
	* Note that the second moment of the force is neglected.
	*/
	void lb_external_forces(double ext_force[3]);

	/** The Lattice Boltzmann collision step.
	* Loop over all lattice sites and perform the collision update.
	* If fluctuations are present, the fluctuating part of the stress tensor
	* is added. The update is only accepted then, if no negative populations
	* occur.
	*/
	void collision();

	/** Propagates the Lattice Boltzmann system for one time step.
	* This function performs the collision step and the streaming step.
	* If external forces are present, they are applied prior to the collisions.
	* If boundaries are present, it also applies the boundary conditions.
	*/
	void lb_propagate();

	/** Performs a full initialization of
	*  the Lattice Boltzmann system. All derived parameters
	*  and the fluid are reset to their default values. */
	void lb_init(double ptemperature, double ptime_step);

	/*  LBE_BCONDS: Sets up extra pointers for periodic boundaries.  */
	void bconds();

	/** Calculates the coupling of MD particles to the LB fluid.
	* This function  is called from \ref force_calc. The force is added
	* to the particle force and the corresponding momentum exchange is
	* applied to the fluid.
	* Note that this function changes the state of the fluid!
	*/
	void calc_particle_lattice_ia(Particles& parts);

	/*  Output data.  */
	void outPut(int t);

	/*  Output data.  */
	void outTxt(int t);

	/*  Output data.  */
	void outN(int t);

	/*  Input data.  */
	void inPut(int t);

	/* Get total mass*/
	double lb_calc_fluid_mass();

	/* Get total movement*/
	void lb_calc_fluid_momentum(double* momentum);

	/** Calculate temperature of the LB fluid.
	* \param result Fluid temperature
	*/
	double lb_calc_fluid_temp();

	/*  LB: output  properties*/
	void outProp(int t, FILE *file_ptr);

	/*  Read properties file  */
	void inProp(int t, FILE	*file_ptr);

private:

	/** local populations of the velocity directions */
	double**  f;

	/** index of the local populations of the velocity direction*/
	int**  f_index;

	/** index of the local populations of the velocity directions 1,3,5,7,10,11,14,15,18*/
	int*  fa_index;
	
	/** index of the local populations of the velocity directions 2,4,6,8,9,12,13,16,17*/
	int*  fb_index;

	/** local momentum */
	double** j;

	/** local density */
	double* rho;

	/** local stress tensor */
	double** pi;

	const BaseMesh&       p_mesh;

	//Size of the mesh.
	int    max_x, max_y, max_z, volume, offset;

	//size of the f_index;
	int    imax_x, imax_y, imax_z, ivolume, ioffset;

	//Size of the fmesh.
	int    fmax_x, fmax_y, fmax_z, fvolume, foffset;

	//surface = halo_grid_volume - grid_volume ;
	int surface;

	//Output file name.
	char* outputFFileName;

	/** lattice constant */
	double agrid;

	/*time_step comes from the discretization.*/
	double time_step;

	/** number of velocities */
	int n_veloc;

	/** speed of sound squared */
	double c_sound_sq;

	//tau
	double tau_v;

	/** Lattice Boltzmann time step
	* This variable is used for convenience instead of having to type lbpar.tau everywhere. */
	double tau;

	/** measures the MD time since the last fluid update */
	double fluidstep;

	/** number density (LJ units) */
	double total_rho;

	/** eigenvalue of the collision operator for relaxation of shear modes */
	double lblambda;
	/** eigenvalue of the collision operator for relaxation of bulk modes */
	double lblambda_bulk;

	/** Flag indicating whether fluctuations are present. */
	int fluct;
	/** amplitude of the fluctuations in the fluid stress tensor */
	double lb_fluct_pref;
	/** amplitude of the bulk fluctuations of the stress tensor */
	double lb_fluct_pref_bulk;
	/** amplitude of the fluctuations in the viscous coupling */
	double lb_coupl_pref;

	/** kinematic viscosity (LJ units) */
	double viscosity;
	/** bulk viscosity (LJ units) */
	double bulk_viscosity;

	/** friction coefficient for viscous coupling (LJ units)
	* Note that the friction coefficient is quite high and may
	* lead to numerical artifacts with low order integrators */
	double friction;

	/** temperature. */
	double temperature;

	/** Calculate local populations from hydrodynamic fields.
	*
	* The mapping is given in terms of the equilibrium distribution.
	*
	* Eq. (2.15) Ladd, J. Fluid Mech. 271, 295-309 (1994)
	* Eq. (4) in Berk Usta, Ladd and Butler, JCP 122, 094902 (2005)
	*
	* @param local_node Pointer to the local lattice site (Input).
	* @param trace      Trace of the local stress tensor (Input).
	* @param trace_eq   Trace of equilibriumd part of local stress tensor (Input).
	*/
	void lb_calc_local_n(int index);

	/** Add fluctuating part to the stress tensor and update the populations.
	*
	* Ladd, J. Fluid Mech. 271, 285-309 (1994).<br>
	* Berk Usta, Ladd and Butler, JCP 122, 094902 (2005).<br>
	* Ahlrichs, PhD-Thesis (2000).
	*   
	* @param local_node Pointer to the local lattice site.
	*/
	void lb_add_fluct_pi(int index);

	/** Collision update of the stress tensor.
	* The stress tensor is relaxed towards the equilibrium.
	*
	* See Eq. (5) in Berk Usta, Ladd and Butler, JCP 122, 094902 (2005)
	*
	* @param local_node Pointer to the local lattice site (Input).
	* @param trace      Trace of local stress tensor (Output).
	* @param trace_eq   Trace of equilibrium part of local stress tensor (Output).
	*/
	void lb_update_local_pi(int index);

	/** Calculate the local fluid fields.
	* The calculation is implemented explicitly for the special case of D3Q19.
	*
	* Original Author: Ahlrichs 06/11/97, 29/03/98
	*
	* @param local_node   The local lattice site.
	* @param calc_pi_flag Flag indicating whether stress tensor should be
	*                     computed.
	*/
	void lb_calc_local_fields(int index, bool calc_pi_flag) {

		double *local_n = f[index];
		double *local_j = j[index];
		double *local_pi = pi[index];
		double avg_rho = total_rho * (agrid * agrid * agrid);

		rho[index] = avg_rho + local_n[0] + local_n[1] + local_n[2]
		+ local_n[3] + local_n[4] + local_n[5] + local_n[6]
		+ local_n[7] + local_n[8] + local_n[9] + local_n[10]
		+ local_n[11] + local_n[12] + local_n[13] + local_n[14]
		+ local_n[15] + local_n[16] + local_n[17] + local_n[18];
		//if (rho[index] !=0.85 )cout << " " << local_n[1] <<" x " <<local_n[2] <<" y " <<local_n[3]<<" z " <<index;

		local_j[0] = local_n[1] - local_n[2] + local_n[7] - local_n[8]
		+ local_n[9] - local_n[10] + local_n[11] - local_n[12]
		+ local_n[13] - local_n[14];
		local_j[1] = local_n[3] - local_n[4] + local_n[7] - local_n[8]
		- local_n[9] + local_n[10] + local_n[15] - local_n[16]
		+ local_n[17] - local_n[18];
		local_j[2] = local_n[5] - local_n[6] + local_n[11] - local_n[12]
		- local_n[13] + local_n[14] + local_n[15] - local_n[16]
		- local_n[17] + local_n[18];

		if (calc_pi_flag) {
			local_pi[0] = avg_rho / 3.0 + local_n[1] + local_n[2]
			+ local_n[7] + local_n[8] + local_n[9] + local_n[10]
			+ local_n[11] + local_n[12] + local_n[13] + local_n[14];
			local_pi[2] = avg_rho / 3.0 + local_n[3] + local_n[4]
			+ local_n[7] + local_n[8] + local_n[9] + local_n[10]
			+ local_n[15] + local_n[16] + local_n[17] + local_n[18];
			local_pi[5] = avg_rho / 3.0 + local_n[5] + local_n[6]
			+ local_n[11] + local_n[12] + local_n[13] + local_n[14]
			+ local_n[15] + local_n[16] + local_n[17] + local_n[18];
			local_pi[1] = local_n[7] - local_n[9] + local_n[8]
			- local_n[10];
			local_pi[3] = local_n[11] + local_n[12] - local_n[13]
			- local_n[14];
			local_pi[4] = local_n[15] + local_n[16] - local_n[17]
			- local_n[18];
		}

	}

	/** Sets the density and momentum on a local lattice site.
	* @param index The index of the lattice site within the local domain (Input)
	* @param rho   Local density of the fluid (Input)
	* @param j     Local momentum of the fluid (Input)
	*/
	void lb_set_local_fields(int index, const double irho, const double *iv,
		const double *ipi);

	/** (Re-)initializes the fluid. */
	void lb_reinit_fluid();

	/** (Re-)allocate memory for the fluid and initialize pointers. */
	void lb_create_fluid();

	/** Initialize lattice.
	*
	* This function initializes the variables describing the lattice
	* layout. Important: The lattice data is <em>not</em> allocated here!
	*
	* \param lattice pointer to the lattice
	* \param agrid   lattice spacing
	* \param tau     time step for lattice dynamics
	*/
	void init_lattice();

	/** (Re-)initializes the derived parameters
	*  for the Lattice Boltzmann system.
	*  The current state of the fluid is unchanged. */
	void lb_reinit_parameters(double ptemperature, double ptime_step);

	/** Coupling of a particle to viscous fluid with Stokesian friction.
	* 
	* Section II.C. Ahlrichs and Duenweg, JCP 111(17):8225 (1999)
	*
	* @param p          The coupled particle (Input).
	* @param force      Coupling force between particle and fluid (Output).
	*/
	void lb_viscous_momentum_exchange(sPart *p, double force[3]);

	/** Map a spatial position to the surrounding lattice sites.
	*
	* This function takes a global spatial position and determines the
	* surrounding elementary cell of the lattice for this position.
	* The distance fraction in each direction is also calculated.
	* <br><em>Remarks:</em>
	* <ul>
	* <li>The spatial position has to be in the local domain.</li>
	* <li>The lattice sites of the elementary cell are returned as local indices</li>
	* </ul>
	* \param lattice    pointer to the lattice (Input)
	* \param pos        spatial position (Input)
	* \param node_index local indices of the surrounding lattice sites (Output)
	* \param delta      distance fraction of pos from the surrounding
	*                   elementary cell, 6 directions (Output)
	*/
	void map_position_to_lattice(const double pos[3], int node_index[8], double delta[6]);

	/** Transfer a certain amount of momentum to a elementray cell of fluid.
	* 
	* Eq. (14) Ahlrichs and Duenweg, JCP 111(17):8225 (1999).
	*
	* @param momentum   Momentum to be transfered to the fluid (lattice
	*                   units) (Input).
	* @param node_index Indices of the sites of the elementary lattice
	*                   cell (Input).
	* @param delta      Weights for the assignment to the single lattice
	*                   sites (Input).
	* @param badrandoms Flag/Counter for the occurrence negative
	*                   populations (Output).
	*/
	void lb_transfer_momentum(const double momentum[3], const int node_index[8], const double delta[6]);
};

#endif /* LBD3Q19_H_ */