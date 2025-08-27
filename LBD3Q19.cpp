// LBD3Q19.cpp: implementation of the D3Q19 lattice Boltzmann class.
//
//////////////////////////////////////////////////////////////////////

#include "LBD3Q19.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

LBD3Q19::LBD3Q19(const BaseMesh& pmesh, const double ptau, const double ptauv, const double pagrid, double pviscosity, double pbulk_viscosity, double pfriction, double ptotal_rho, double ptemperature, double ptime_step, char* OFFileName)
:p_mesh( pmesh )
,tau( ptau )
,tau_v( ptauv )
,agrid( pagrid )
,viscosity( pviscosity )
,bulk_viscosity( pbulk_viscosity )
,friction( pfriction )
,total_rho( ptotal_rho )
,outputFFileName(OFFileName)
{
	n_veloc = 19;

	c_sound_sq = 1./3.;

	if (viscosity < 0.0) {
		//errtext = runtime_error(128);
		//ERROR_SPRINTF(errtext,"{101 Lattice Boltzmann fluid viscosity not set} ");
		throw ComException(" Lattice Boltzmann fluid viscosity not set. ", __FILE__, __LINE__ );
	}
	if (total_rho < 0.0) {
		//errtext = runtime_error(128);
		//ERROR_SPRINTF(errtext,"{100 Lattice Boltzmann fluid density not set} ");
		throw ComException(" Lattice Boltzmann fluid density not set. ", __FILE__, __LINE__ );
	}
	if (agrid < 0.0) {
		//errtext = runtime_error(128);
		//ERROR_SPRINTF(errtext,"{098 Lattice Boltzmann agrid not set} ");
		throw ComException(" Lattice Boltzmann agrid not set. ", __FILE__, __LINE__ );
	}
	if (tau < 0.0) {
		//errtext = runtime_error(128);
		//ERROR_SPRINTF(errtext,"{099 Lattice Boltzmann time step not set} ");
		throw ComException(" Lattice Boltzmann time step not set. ", __FILE__, __LINE__ );
	}

	/** amplitude of the fluctuations in the fluid stress tensor */
	lb_fluct_pref = 0.0;
	/** amplitude of the bulk fluctuations of the stress tensor */
	lb_fluct_pref_bulk = 0.0;

	/** eigenvalue of the collision operator for relaxation of shear modes */
	lblambda = -1;
	/** eigenvalue of the collision operator for relaxation of bulk modes */
	lblambda_bulk = -1;

	/** measures the MD time since the last fluid update */
	fluidstep=0.0;

	lb_init(ptemperature, ptime_step);
}

LBD3Q19::~LBD3Q19(void)
{
	int i;
	for(i = 0;i < fvolume;i++){
		delete[] f[i], pi[i], j[i];
	}
	for(i = 0;i < ivolume;i++){
		delete[] f_index[i];
	}
	delete[] f, pi, rho, j, fa_index, fb_index;
}

/***********************************************************************/
/** \name Streaming step */
/***********************************************************************/
/*@{*/

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
void LBD3Q19::stream()
{
	int k;
	int indexA,indexB;
	int* rindex;

	double* local_fA;
	double* local_fB;
	double* remot_f;

	/* bottom up sweep */
	for (k=0;k<ivolume;k++) {
		indexA = fa_index[k];
		indexB = fb_index[k];
		local_fA = f[indexA];
		local_fB = f[indexB];

		rindex = f_index[k];

		/* propagation to higher indices */		
		remot_f = f[rindex[1]];
		local_fA[1] = remot_f[1];

		remot_f = f[rindex[3]];
		local_fA[3] = remot_f[3];

		remot_f = f[rindex[5]];
		local_fA[5] = remot_f[5];

		remot_f = f[rindex[7]];
		local_fA[7] = remot_f[7];

		remot_f = f[rindex[10]];
		local_fA[10] = remot_f[10];

		remot_f = f[rindex[11]];
		local_fA[11] = remot_f[11];

		remot_f = f[rindex[14]];
		local_fA[14] = remot_f[14];

		remot_f = f[rindex[15]];
		local_fA[15] = remot_f[15];

		remot_f = f[rindex[18]];
		local_fA[18] = remot_f[18];

		/* propagation to lower indices */
		remot_f = f[rindex[2]];
		local_fB[2] = remot_f[2];

		remot_f = f[rindex[4]];
		local_fB[4] = remot_f[4];

		remot_f = f[rindex[6]];
		local_fB[6] = remot_f[6];

		remot_f = f[rindex[8]];
		local_fB[8] = remot_f[8];

		remot_f = f[rindex[9]];
		local_fB[9] = remot_f[9];

		remot_f = f[rindex[12]];
		local_fB[12] = remot_f[12];

		remot_f = f[rindex[13]];
		local_fB[13] = remot_f[13];

		remot_f = f[rindex[16]];
		local_fB[16] = remot_f[16];

		remot_f = f[rindex[17]];
		local_fB[17] = remot_f[17];
	}
}

/***********************************************************************/
/** \name External forces */
/***********************************************************************/
/*@{*/
/** Apply external forces to the fluid.
*
* Eq. (28) Ladd and Verberg, J. Stat. Phys. 104(5/6):1191 (2001).
* Note that the second moment of the force is neglected.
*/
void LBD3Q19::lb_external_forces(double ext_force[3]) {

	int x, y, z, index;
	double* local_f;
	double* local_j;
	double	delta_j[3];

	for (z=1; z<=max_z; z++) {
		for (y=1; y<=max_y; y++) {
			for (x=1; x<=max_x; x++) {
				index = (x + fmax_x*(y + fmax_y*z));
				/* calculate momentum due to ext_force in lattice units */
				/* ext_force is the force per volume in LJ units */
				delta_j[0] = ext_force[0]*tau*tau*agrid*agrid;
				delta_j[1] = ext_force[1]*tau*tau*agrid*agrid;
				delta_j[2] = ext_force[2]*tau*tau*agrid*agrid;

				//cout << delta_j[0] <<" delta_j[0]"<< delta_j[1] <<" delta_j[1]"<< delta_j[2] <<" delta_j[2]"<< index<<endl;

				local_j = j[index];
				local_j[0] += delta_j[0];
				local_j[1] += delta_j[1];
				local_j[2] += delta_j[2];

				local_f = f[index];

				local_f[1]  +=   1./6. * delta_j[0];
				local_f[2]  += - 1./6. * delta_j[0];
				local_f[3]  +=   1./6. * delta_j[1];
				local_f[4]  += - 1./6. * delta_j[1];
				local_f[5]  +=   1./6. * delta_j[2];
				local_f[6]  += - 1./6. * delta_j[2];
				local_f[7]  +=   1./12. * (delta_j[0]+delta_j[1]);
				local_f[8]  += - 1./12. * (delta_j[0]+delta_j[1]);
				local_f[9]  +=   1./12. * (delta_j[0]-delta_j[1]);
				local_f[10] += - 1./12. * (delta_j[0]-delta_j[1]);
				local_f[11] +=   1./12. * (delta_j[0]+delta_j[2]);
				local_f[12] += - 1./12. * (delta_j[0]+delta_j[2]);
				local_f[13] +=   1./12. * (delta_j[0]-delta_j[1]);
				local_f[14] += - 1./12. * (delta_j[0]-delta_j[1]);
				local_f[15] +=   1./12. * (delta_j[1]+delta_j[2]);
				local_f[16] += - 1./12. * (delta_j[1]+delta_j[2]);
				local_f[17] +=   1./12. * (delta_j[1]-delta_j[2]);
				local_f[18] += - 1./12. * (delta_j[1]-delta_j[2]);
			}
		}
	}
}
/*@}*/

/** The Lattice Boltzmann collision step.
* Loop over all lattice sites and perform the collision update.
* If fluctuations are present, the fluctuating part of the stress tensor
* is added. The update is only accepted then, if no negative populations
* occur.
*/
void LBD3Q19::collision() {

	int index, x, y, z;
	/* loop over all nodes (halo excluded) */
	index = foffset;
	for (z=1;z<=max_z;z++) {
		for (y=1;y<=max_y;y++) {
			for (x=1;x<=max_x;x++) {

				lb_calc_local_fields(index,true);

				lb_update_local_pi(index);

				if (fluct) lb_add_fluct_pi(index);

				lb_calc_local_n(index);

				++index;
			}
			index += 2;
		}
		index += 2*fmax_x;
	}
}
/*@}*/

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
void LBD3Q19::lb_calc_local_n(int index) {

	double *local_n   = f[index];
	double local_rho = rho[index];
	double *local_j   = j[index];
	double *local_pi  = pi[index];
	double trace;
	const double rhoc_sq = local_rho*c_sound_sq;
	const double avg_rho = total_rho*(agrid*agrid*agrid);

	/* see Eq. (4) in Berk Usta, Ladd and Butler, JCP 122, 094902 (2005) */

	/* reduce the pressure tensor to the part needed here */
	local_pi[0] -= rhoc_sq;
	local_pi[2] -= rhoc_sq;
	local_pi[5] -= rhoc_sq;

	trace = local_pi[0] + local_pi[2] + local_pi[5];

	double rho_times_coeff;
	double tmp1,tmp2;

	/* update the q=0 sublattice */
	local_n[0] = 1./3. * (local_rho-avg_rho) - 1./2.*trace;

	/* update the q=1 sublattice */
	rho_times_coeff = 1./18. * (local_rho-avg_rho);

	local_n[1] = rho_times_coeff + 1./6.*local_j[0] + 1./4.*local_pi[0] - 1./12.*trace;
	local_n[2] = rho_times_coeff - 1./6.*local_j[0] + 1./4.*local_pi[0] - 1./12.*trace;
	local_n[3] = rho_times_coeff + 1./6.*local_j[1] + 1./4.*local_pi[2] - 1./12.*trace;
	local_n[4] = rho_times_coeff - 1./6.*local_j[1] + 1./4.*local_pi[2] - 1./12.*trace;
	local_n[5] = rho_times_coeff + 1./6.*local_j[2] + 1./4.*local_pi[5] - 1./12.*trace;
	local_n[6] = rho_times_coeff - 1./6.*local_j[2] + 1./4.*local_pi[5] - 1./12.*trace;

	/* update the q=2 sublattice */
	rho_times_coeff = 1./36. * (local_rho-avg_rho);

	tmp1 = local_pi[0] + local_pi[2];
	tmp2 = 2.0*local_pi[1];

	local_n[7]  = rho_times_coeff + 1./12.*(local_j[0]+local_j[1]) + 1./8.*(tmp1+tmp2) - 1./24.*trace;
	local_n[8]  = rho_times_coeff - 1./12.*(local_j[0]+local_j[1]) + 1./8.*(tmp1+tmp2) - 1./24.*trace;
	local_n[9]  = rho_times_coeff + 1./12.*(local_j[0]-local_j[1]) + 1./8.*(tmp1-tmp2) - 1./24.*trace;
	local_n[10] = rho_times_coeff - 1./12.*(local_j[0]-local_j[1]) + 1./8.*(tmp1-tmp2) - 1./24.*trace;

	tmp1 = local_pi[0] + local_pi[5];
	tmp2 = 2.0*local_pi[3];

	local_n[11] = rho_times_coeff + 1./12.*(local_j[0]+local_j[2]) + 1./8.*(tmp1+tmp2) - 1./24.*trace;
	local_n[12] = rho_times_coeff - 1./12.*(local_j[0]+local_j[2]) + 1./8.*(tmp1+tmp2) - 1./24.*trace;
	local_n[13] = rho_times_coeff + 1./12.*(local_j[0]-local_j[2]) + 1./8.*(tmp1-tmp2) - 1./24.*trace;
	local_n[14] = rho_times_coeff - 1./12.*(local_j[0]-local_j[2]) + 1./8.*(tmp1-tmp2) - 1./24.*trace;

	tmp1 = local_pi[2] + local_pi[5];
	tmp2 = 2.0*local_pi[4];

	local_n[15] = rho_times_coeff + 1./12.*(local_j[1]+local_j[2]) + 1./8.*(tmp1+tmp2) - 1./24.*trace;
	local_n[16] = rho_times_coeff - 1./12.*(local_j[1]+local_j[2]) + 1./8.*(tmp1+tmp2) - 1./24.*trace;
	local_n[17] = rho_times_coeff + 1./12.*(local_j[1]-local_j[2]) + 1./8.*(tmp1-tmp2) - 1./24.*trace;
	local_n[18] = rho_times_coeff - 1./12.*(local_j[1]-local_j[2]) + 1./8.*(tmp1-tmp2) - 1./24.*trace;

	/* restore the pressure tensor to the full part */
	local_pi[0] += rhoc_sq;
	local_pi[2] += rhoc_sq;
	local_pi[5] += rhoc_sq;
}

/** Add fluctuating part to the stress tensor and update the populations.
*
* Ladd, J. Fluid Mech. 271, 285-309 (1994).<br>
* Berk Usta, Ladd and Butler, JCP 122, 094902 (2005).<br>
* Ahlrichs, PhD-Thesis (2000).
*   
* @param local_node Pointer to the local lattice site.
*/
void LBD3Q19::lb_add_fluct_pi(int index) {

	double *local_pi = pi[index];
	double tmp, sum=0.0;

	const double pref1 = sqrt(2.0) * lb_fluct_pref;

	/* off-diagonal components */
	local_pi[1] += lb_fluct_pref * (d_random()-0.5);
	local_pi[3] += lb_fluct_pref * (d_random()-0.5);
	local_pi[4] += lb_fluct_pref * (d_random()-0.5);

	/* diagonal components */
	tmp = (d_random()-0.5);
	sum += tmp;
	local_pi[0] += pref1 * tmp;
	tmp = (d_random()-0.5);
	sum += tmp;
	local_pi[2] += pref1 * tmp;
	tmp = (d_random()-0.5);
	sum += tmp;
	local_pi[5] += pref1 * tmp;

	/* make shear modes traceless and add bulk fluctuations on the trace */
	sum *= (lb_fluct_pref_bulk/sqrt(3.0) - pref1/3.0);
	local_pi[0] += sum;
	local_pi[2] += sum;
	local_pi[5] += sum;
}

/** Collision update of the stress tensor.
* The stress tensor is relaxed towards the equilibrium.
*
* See Eq. (5) in Berk Usta, Ladd and Butler, JCP 122, 094902 (2005)
*
* @param local_node Pointer to the local lattice site (Input).
* @param trace      Trace of local stress tensor (Output).
* @param trace_eq   Trace of equilibrium part of local stress tensor (Output).
*/
void LBD3Q19::lb_update_local_pi(int index) {

	const double local_rho = rho[index];
	double *local_j  = j[index];
	double *local_pi = pi[index];
	double local_pi_eq[6];
	double trace, trace_eq;
	double tmp;

	const double rhoc_sq = local_rho*c_sound_sq;
	const double onepluslambda = 1.0 + lblambda;

	/* calculate the equilibrium part of the pressure tensor */
	local_pi_eq[0] = rhoc_sq + local_j[0]*local_j[0]/local_rho;
	tmp = local_j[1]/local_rho;
	local_pi_eq[1] = local_j[0]*tmp;
	local_pi_eq[2] = rhoc_sq + local_j[1]*tmp;
	tmp = local_j[2]/local_rho;
	local_pi_eq[3] = local_j[0]*tmp;
	local_pi_eq[4] = local_j[1]*tmp;
	local_pi_eq[5] = rhoc_sq + local_j[2]*tmp;

	/* calculate the traces */
	trace_eq = local_pi_eq[0] + local_pi_eq[2] + local_pi_eq[5];
	trace = local_pi[0] + local_pi[2] + local_pi[5];

	/* relax the local pressure tensor */
	local_pi[0] = local_pi_eq[0] + onepluslambda*(local_pi[0] - local_pi_eq[0]);
	local_pi[1] = local_pi_eq[1] + onepluslambda*(local_pi[1] - local_pi_eq[1]);
	local_pi[2] = local_pi_eq[2] + onepluslambda*(local_pi[2] - local_pi_eq[2]);
	local_pi[3] = local_pi_eq[3] + onepluslambda*(local_pi[3] - local_pi_eq[3]);
	local_pi[4] = local_pi_eq[4] + onepluslambda*(local_pi[4] - local_pi_eq[4]);
	local_pi[5] = local_pi_eq[5] + onepluslambda*(local_pi[5] - local_pi_eq[5]);  
	tmp = 1./3.*(lblambda_bulk-lblambda)*(trace - trace_eq);
	local_pi[0] += tmp;
	local_pi[2] += tmp;
	local_pi[5] += tmp;
}

/***********************************************************************/
/** \name Integration step for the lattice Boltzmann fluid             */
/***********************************************************************/
/*@{*/
/** Propagate the Lattice Boltzmann dynamics.
* This function is called from the integrator. Since the time step
* for the lattice dynamics can be coarser than the MD time step,
* we monitor the time since the last lattice update.
*/
void LBD3Q19::lb_propagate() {

	fluidstep+=time_step ;
	//cout <<fluidstep<<" fluidstep "<<tau<< " tau "<<endl;
	if (fluidstep>=tau) {
		//cout << "test"<<endl;
		fluidstep=0.0 ;
		//outN(200+fluidstep);
		//outN(100+fluidstep);

		/* collision step */
		collision();

		/* exchange halo regions */
		//update_halo_comm->halo_communication();
		//outN(100+fluidstep);
		/* streaming step */
		stream();
		//outN(200+fluidstep);

	}

}

/** Performs a full initialization of
*  the Lattice Boltzmann system. All derived parameters
*  and the fluid are reset to their default values. */
void LBD3Q19::lb_init(double ptemperature, double ptime_step) {

	/* initialize derived parameters */
	lb_reinit_parameters(ptemperature, ptime_step);

	/* initialize the local lattice domain */
	init_lattice();

	/* allocate memory for data structures */
	lb_create_fluid();

	/* setup the initial particle velocity distribution */
	lb_reinit_fluid();

}

/** Sets the hydrodynamic fields on a local lattice site.
* @param index The index of the lattice site within the local domain (Input)
* @param rho   Local density of the fluid (Input)
* @param v     Local velocity of the fluid (Input)
*/
void LBD3Q19::lb_set_local_fields(int index, const double irho, const double *iv, const double *ipi) {

	double *local_j = j[index];
	double *local_pi = pi[index];

	rho[index] = irho;

	local_j[0] = irho * iv[0];
	local_j[1] = irho * iv[1];
	local_j[2] = irho * iv[2];

	local_pi[0] = ipi[0];
	local_pi[1] = ipi[1];
	local_pi[2] = ipi[2];
	local_pi[3] = ipi[3];
	local_pi[4] = ipi[4];
	local_pi[5] = ipi[5];

	/* calculate populations according to equilibrium distribution */
	lb_calc_local_n(index);

}
/*@}*/

/** (Re-)initializes the fluid according to the given value of rho. */
void LBD3Q19::lb_reinit_fluid() {

	int k ;

	/* default values for fields in lattice units */
	double prho = total_rho*(agrid*agrid*agrid) ;
	double pv[3] = { 0., 0., 0. };
	double ppi[6] = { prho*c_sound_sq, 0., prho*c_sound_sq, 0., 0., prho*c_sound_sq };

	for (k=0;k<fvolume;k++) {
		lb_set_local_fields(k,prho,pv,ppi);
		//cout << k << " k " << rho << " rho " << v[0]+v[1]+v[2] << " v " << pi[0] << " "<< pi[1] << " "<< pi[2] << " "<< pi[3] << " "<< pi[4]<< " "<<pi[5]<< " "<<endl;
	}
}

/** (Re-)allocate memory for the fluid and initialize pointers. */
void LBD3Q19::lb_create_fluid() {

	int i;

	f =  new double* [fvolume];
	f_index =  new int* [ivolume];
	fa_index =  new int [ivolume];
	fb_index =  new int [ivolume];
	j =  new double* [fvolume];
	rho =  new double [fvolume];
	pi = new double* [fvolume];
	for(i = 0;i < fvolume;i++){
		f[i] = new double [n_veloc];
		f_index[i] =  new int [n_veloc];
		pi[i] = new double [6];
		j[i] = new double [3];
	}
}

/** Initialize lattice.
*
* This function initializes the variables describing the lattice
* layout. Important: The lattice data is <em>not</em> allocated here!
*
* \param lattice pointer to the lattice
* \param agrid   lattice spacing
* \param tau     time step for lattice dynamics
*/
void LBD3Q19::init_lattice() {
	max_x = p_mesh.getNx()/agrid;
	max_y = p_mesh.getNy()/agrid;
	max_z = p_mesh.getNz()/agrid;
	volume = max_x * max_y * max_z;
	offset = (1 + max_x*(1 + max_y*1));

	imax_x = max_x + 1;
	imax_y = max_y + 1;
	imax_z = max_z + 1;
	ivolume = imax_x * imax_y * imax_z;
	ioffset = (1 + imax_x*(1 + imax_y*1));

	fmax_x = max_x + 2;
	fmax_y = max_y + 2;
	fmax_z = max_z + 2;
	fvolume = fmax_x * fmax_y * fmax_z;
	foffset = (1 + fmax_x*(1 + fmax_y*1));

	surface = fvolume - volume;
}

/** (Re-)initializes the fluid. */
void LBD3Q19::lb_reinit_parameters(double ptemperature, double ptime_step) {

	time_step = ptime_step;
	temperature = ptemperature;

	/* Eq. (3) Ahlrichs and Duenweg, JCP 111(17):8225 (1999). */
	lblambda = -2./(6.*viscosity*tau/(agrid*agrid)+1.);

	if (bulk_viscosity > 0.0) {
		lblambda_bulk = -2./(9.*bulk_viscosity*tau/(agrid*agrid)+1.);
	}

	if (temperature > 0.0) {  /* fluctuating hydrodynamics ? */

		fluct = 1 ;

		/* lb_fluct_pref is stored in lattice units (pressure)
		* Eq. (7) Ahlrichs and Duenweg, JCP 111(17):8225 (1999).
		* The factor 12 comes from the fact that we use random numbers
		* from -0.5 to 0.5 (equally distributed) which have variance 1/12.
		*/
		lb_fluct_pref = sqrt(12.*2.*viscosity*total_rho*temperature*SQR(lblambda)*tau*tau*tau/agrid);
		cout <<" "<<lb_fluct_pref<< " lb_fluct_pref "<<viscosity<<" viscosity "<<total_rho<<" rho "<<temperature<<" temperature "<<endl;
		cout <<" "<<SQR(lblambda)<<" sqr(lblambda) "<<tau<<" tau "<<agrid<<" argid "<<endl;
		if (bulk_viscosity > 0.0) {
			lb_fluct_pref_bulk = sqrt(12.*2.*bulk_viscosity*total_rho*temperature*SQR(lblambda_bulk)*tau*tau*tau/agrid);
			cout <<" "<<lb_fluct_pref_bulk<<" lb_fluct_pref_bulk "<<bulk_viscosity<<" lbpar.bulk_viscosity "<<endl;
		}

		//LB_TRACE(fprintf(stderr,"%d: lb_fluct_pref=%f lb_fluct_pref_bulk=%f (temp=%f, lambda=%f, lambda_v=%f, tau=%f, agrid=%f)\n",this_node,lb_fluct_pref,lb_fluct_pref_bulk,temperature,lblambda,lblambda_bulk,tau,agrid));

		/* lb_coupl_pref is stored in MD units (force)
		* Eq. (16) Ahlrichs and Duenweg, JCP 111(17):8225 (1999).
		* The factor 12 comes from the fact that we use random numbers
		* from -0.5 to 0.5 (equally distributed) which have variance 1/12.
		* time_step comes from the discretization.
		*/
		lb_coupl_pref = sqrt(12.*2.*friction*temperature/time_step); 

		//LB_TRACE(fprintf(stderr,"%d: lb_coupl_pref=%f (temp=%f, friction=%f, time_step=%f)\n",this_node,lb_coupl_pref,temperature,lbpar.friction,time_step));

	} else {
		/* no fluctuations at zero temperature */
		fluct = 0 ;
		lb_fluct_pref = 0.0;
		lb_coupl_pref = 0.0;
	}

}


/*  Output data.  */
void LBD3Q19::outPut (int t)
{
	char	filename[1024];
	FILE	*file_ptr;

	int index;
	int i,m, k;

	sprintf(filename,"%s%03d.dat" , outputFFileName, t) ;
	file_ptr = fopen(filename, "w");
	if(!file_ptr) throw ComException("Cann not open the output data of fliuds. ", __FILE__, __LINE__ );

	fwrite (&(t), sizeof(int), 1, file_ptr);
	//cout << lblattice->halo_grid_volume << " lblattice->halo_grid_volume "<<endl;
	for (index=0;index<fvolume;index++) {

		//double *local_rho = rho[index];
		double *local_j   = j[index];
		double *local_pi  = pi[index];

		//Comment it to reduce the calculation.
		//lb_calc_local_fields(index, true);

		fwrite (&(rho[index]), sizeof(double), 1, file_ptr);
		//cout << *local_rho << " rho ";
		m = 0;
		for (i=0;i<3;i++) {
			fwrite (&(local_j[i]), sizeof(double), 1, file_ptr);
			//if (index < 2200) cout << " "<<local_j[i] << " j" << i << " ";
			for (k=0;k<=i;k++) {
				fwrite (&(local_pi[m]), sizeof(double), 1, file_ptr);
				//cout << "local_pi[m]" << local_pi[m] << " m " << m << " ";
				m++;
			}
		}
		//if (index == 32) cout << index << " index " << endl;
	}
	//cout << endl;
	fclose (file_ptr);  file_ptr = 0;
}

/*  Input data.  */
void LBD3Q19::inPut (int t)
{
	char	filename[1024];
	FILE	*file_ptr;
	int index; double prho; double pj[3]; double ppi[6];
	int m, i, k;
	int dt;

	//sprintf(filename,"%s.dat" , outputFFileName) ;
	sprintf(filename,"%s%03d.dat" , outputFFileName, t) ;
	file_ptr = fopen(filename, "rb");
	if(file_ptr != 0) 
	{
		cout << "open " << filename << " to read lbm details." << endl;

		fread ( &(dt), sizeof(int), 1, file_ptr);
		//cout << lblattice->halo_grid_volume << " lblattice->halo_grid_volume "<<endl;
		for (index=0;index<fvolume;index++) {

			fread (&prho, sizeof(double), 1, file_ptr);
			//cout << prho << " rho ";
			m = 0;
			for (i=0;i<3;i++) {
				fread (&(pj[i]), sizeof(double), 1, file_ptr);
				//cout << pj[i]<<" j"<< i <<" ";
				for (k=0;k<=i;k++) {
					//io_ptr = (void *) pi;
					fread (&(ppi[m]), sizeof(double), 1, file_ptr);
					//cout << ppi[m]<<" pi"<<m <<" ";
					m++;
				}
			}
			//cout <<" index "<< index;
			//lb_set_local_fields(index,prho,pj,ppi);
		}
		//cout <<endl;
		fclose (file_ptr);  file_ptr = 0;
	}
}

/*  LBE_BCONDS: Sets up extra pointers for periodic boundaries.  */
void LBD3Q19::bconds()
{
	int* local_f;
	int x, y, z;
	int index, bondary;
	int moveA,moveB;


	for (z = 0; z < imax_z; z++){
		for (y = 0; y <  imax_y; y++){
			for (x = 0; x <  imax_x; x++){
				bondary = 0;
				index = (x + imax_x*(y + imax_y*z));
				local_f = f_index[index];
				fb_index[index] = (x + fmax_x*(y + fmax_y*z));

				//next[2]  = -  1;                  // (-1, 0, 0)
				if (x == max_x){
					moveA = - max_x;
					bondary += 1;
				} else {
					moveA = 1;;
				}
				index = ((x+moveA) + fmax_x*(y + fmax_y*z));
				local_f[2]= index;


				//next[4]  = -  yperiod;            // ( 0,-1, 0)
				if (y == max_y){
					moveA = - max_y;
					bondary += 2;
				} else {
					moveA = 1;
				}
				index = (x + fmax_x*((y+moveA) + fmax_y*z));
				local_f[4]= index;

				//next[6]  = -  zperiod;            // ( 0, 0,-1)
				if (z == max_z){
					moveA = - max_z;
					bondary += 4;
				} else {
					moveA = 1;
				}
				index = (x + fmax_x*(y + fmax_y*(z+moveA)));
				local_f[6]= index;

				//next[8]  = -  (1+yperiod);        // (-1,-1, 0)
				if (x == max_x){
					moveA = - max_x+1;
				} else {
					moveA = 1;
				}
				if (y == max_y){
					moveB = - max_y;
					moveA = 0;
				} else {
					moveB = 1;
				}
				index = ((x+moveA) + fmax_x*((y+moveB) + fmax_y*z));
				local_f[8]= index;

				//next[9]  =    (1-yperiod);        // ( 1,-1, 0)
				if (x == 1){
					moveA = max_x-1;
				} else {
					moveA = -1;;
				}
				if (y == max_y){
					moveB = - max_y;
					moveA = 0;
				} else {
					moveB = 1;
				}
				index = ((x+moveA) + fmax_x*((y+moveB) + fmax_y*z));
				local_f[9]= index;

				//next[12] = -  (1+zperiod);       // (-1, 0,-1)
				if (x == max_x){
					moveA = - max_x+1;
				} else {
					moveA = 1;
				}
				if (z == max_z){
					moveB = - max_z;
					moveA = 0;
				} else {
					moveB = 1;
				}
				index = ((x+moveA) + fmax_x*(y + fmax_y*(z+moveB)));
				local_f[12]= index;

				//next[13] =    (1-zperiod);       // ( 1, 0,-1)
				if (x == 1){
					moveA = max_x-1;
				} else {
					moveA = -1;;
				}
				if (z == max_z){
					moveB = - max_z;
					moveA = 0;
				} else {
					moveB = 1;
				}
				index = ((x+moveA) + fmax_x*(y + fmax_y*(z+moveB)));
				local_f[13]= index;

				//next[16] = -  (yperiod+zperiod); // ( 0,-1,-1)
				if (y == max_y){
					moveA = - max_y+1;
				} else {
					moveA = 1;
				}
				if (z == max_z){
					moveB = - max_z;
					moveA = 0;
				} else {
					moveB = 1;
				}
				index = (x + fmax_x*((y+moveA) + fmax_y*(z+moveB)));
				local_f[16]= index;

				//next[17] =    (yperiod-zperiod); // ( 0, 1,-1)
				if (y == 1){
					moveA = max_y-1;
				} else {
					moveA = -1;;
				}
				if (z == max_z){
					moveB = - max_z;
					moveA = 0;
				} else {
					moveB = 1;
				}
				index = (x + fmax_x*((y+moveA) + fmax_y*(z+moveB)));
				local_f[17]= index;

				////next[19]  = -  (1+yperiod);        // (-1,-1,-1)
				//if (x == max_x){
				//	moveA = - max_x+1;
				//} else {
				//	moveA = 1;
				//}
				//if (y == max_y){
				//	moveB = - max_y+1;
				//	moveA = 0;
				//} else {
				//	moveB = 1;
				//}

				//if (z == max_z){
				//	moveB = - max_z+1;
				//	moveA = 0;
				//} else {
				//	moveC = 1;
				//}
				//index = ((x+moveA) + fmax_x*((y+moveB) + fmax_y*(z+moveC)));
				//local_f[19]= index;

				local_f[0]=bondary;
			}
		}
	}

	for (z = imax_z; z > 0; z--){
		for (y = imax_y; y > 0; y--){
			for (x = imax_x; x > 0; x--){
				index = ((imax_x-x) + imax_x*((imax_y-y) + imax_y*(imax_z-z)));
				local_f = f_index[index];
				fa_index[index] = (x + fmax_x*(y + fmax_y*z));

				//next[1]  =    1;                  // ( 1, 0, 0) +
				if (x == 1){
					moveA = max_x;
				} else {
					moveA = -1;;
				}
				index = ((x+moveA) + fmax_x*(y + fmax_y*z));
				local_f[1]= index;


				//next[3]  =    yperiod;            // ( 0, 1, 0) +
				if (y == 1){
					moveA = max_y;
				} else {
					moveA = -1;
				}
				index = (x + fmax_x*((y+moveA) + fmax_y*z));
				local_f[3]= index;

				//next[5]  =    zperiod;            // ( 0, 0, 1) +
				if (z == 1){
					moveA = max_z;
				} else {
					moveA = -1;
				}
				index = (x + fmax_x*(y + fmax_y*(z+moveA)));
				local_f[5]= index;

				//next[7]  =    (1+yperiod);        // ( 1, 1, 0) +
				if (x == 1){
					moveA = max_x-1;
				} else {
					moveA = -1;
				}
				if (y == 1){
					moveB = max_y;
					moveA = 0;
				} else {
					moveB = -1;
				}
				index = ((x+moveA) + fmax_x*((y+moveB) + fmax_y*z));
				local_f[7]= index;

				//next[10] = -  (1-yperiod);       // (-1, 1, 0) +
				if (x == max_x){
					moveA = -max_x+1;
				} else {
					moveA = 1;
				}
				if (y == 1){
					moveB = max_y;
					moveA = 0;
				} else {
					moveB = -1;
				}
				index = ((x+moveA) + fmax_x*((y+moveB) + fmax_y*z));
				local_f[10]= index;

				//next[11] =    (1+zperiod);       // ( 1, 0, 1) +
				if (x == 1){
					moveA = max_x-1;
				} else {
					moveA = -1;
				}
				if (z == 1){
					moveB = max_z;
					moveA = 0;
				} else {
					moveB = -1;
				}
				index = ((x+moveA) + fmax_x*(y + fmax_y*(z+moveB)));
				local_f[11]= index;

				//next[14] = -  (1-zperiod);       // (-1, 0, 1) +
				if (x == max_x){
					moveA = -max_x+1;
				} else {
					moveA = 1;
				}
				if (z == 1){
					moveB = max_z;
					moveA = 0;
				} else {
					moveB = -1;
				}
				index = ((x+moveA) + fmax_x*(y + fmax_y*(z+moveB)));
				local_f[14]= index;

				//next[15] =    (yperiod+zperiod); // ( 0, 1, 1) +
				if (y == 1){
					moveA = max_y-1;
				} else {
					moveA = -1;
				}
				if (z == 1){
					moveB = max_z;
					moveA = 0;
				} else {
					moveB = -1;
				}
				index = (x + fmax_x*((y+moveA) + fmax_y*(z+moveB)));
				local_f[15]= index;

				//next[18] = -  (yperiod-zperiod); // ( 0,-1, 1) +
				if (y == max_y){
					moveA = -max_y+1;
				} else {
					moveA = 1;
				}
				if (z == 1){
					moveB = max_z;
					moveA = 0;
				} else {
					moveB = -1;
				}
				index = (x + fmax_x*((y+moveA) + fmax_y*(z+moveB)));
				local_f[18]= index;
			}
		}
	}
}

/** Calculate particle lattice interactions.
* So far, only viscous coupling with Stokesian friction is
* implemented.
* Include all particle-lattice forces in this function.
* The function is called from \ref force_calc.
*
* Parallelizing the fluid particle coupling is not straightforward
* because drawing of random numbers makes the whole thing nonlocal.
* One way to do it is to treat every particle only on one node, i.e.
* the random numbers need not be communicated. The particles that are 
* not fully inside the local lattice are taken into account via their
* ghost images on the neighbouring nodes. But this requires that the 
* correct values of the surrounding lattice nodes are available on 
* the respective node, which means that we have to communicate the 
* halo regions before treating the ghost particles. Moreover, after 
* determining the ghost couplings, we have to communicate back the 
* halo region such that all local lattice nodes have the correct values.
* Thus two communication phases are involved which will most likely be 
* the bottleneck of the computation.
*
* Another way of dealing with the particle lattice coupling is to 
* treat a particle and all of it's images explicitly. This requires the
* communication of the random numbers used in the calculation of the 
* coupling force. The problem is now that, if random numbers have to 
* be redrawn, we cannot efficiently determine which particles and which 
* images have to be re-calculated. We therefore go back to the outset
* and go through the whole system again until no failure occurs during
* such a sweep. In the worst case, this is very inefficient because
* many things are recalculated although they actually don't need.
* But we can assume that this happens extremely rarely and then we have
* on average only one communication phase for the random numbers, which
* probably makes this method preferable compared to the above one.
*/
void LBD3Q19::calc_particle_lattice_ia(Particles& parts) {

	int i, k, c, np;
	Cell* cell ;
	sPart* p ;
	double force[3];
	int ncell_size = parts.local_cells.size();

	//sPart gp ;
	//Cell* gcell ;

	for (k=0;k<fvolume;k++) {
		lb_calc_local_fields(k, false);
	}

	/* draw random numbers for local particles */
	for (c=0;c<ncell_size;c++) {
		cell = parts.local_cells[c] ;
		np = cell->p_objects.size() ;
		for (i=0;i<np;i++) {
			p = cell->p_objects[i];
			p->f_random[0] = -lb_coupl_pref*(d_random()-0.5);
			p->f_random[1] = -lb_coupl_pref*(d_random()-0.5);
			p->f_random[2] = -lb_coupl_pref*(d_random()-0.5);

		}
	}

	/* local cells */
	for (c=0;c<ncell_size;c++) {
		cell = parts.local_cells[c] ;
		//p = cell->part ;
		np = cell->p_objects.size() ;

		for (i=0;i<np;i++) {
			p = cell->p_objects[i];

			lb_viscous_momentum_exchange(p,force) ;

			/* add force to the particle */
			p->f[0] += force[0];
			p->f[1] += force[1];
			p->f[2] += force[2];

			//ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LB f = (%.6e,%.3e,%.3e)\n",this_node,p->f.f[0],p->f.f[1],p->f.f[2]));
		}
	}

	//ncell_size = parts.ghost_cells.size();
	//	/* local cells */
	//for (c=0;c<ncell_size;c++) {
	//	gcell = parts.ghost_cells[c] ;
	//	cell = parts.local_cells[gcell->noGCell];
	//	//p = cell->part ;
	//	np = cell->p_objects.size() ;

	//	for (i=0;i<np;i++) {
	//		p = cell->p_objects[i];

	//		gp = *p;

	//		gp.p[0] += gcell->gPMove[0];
	//		gp.p[1] += gcell->gPMove[1];
	//		gp.p[2] += gcell->gPMove[2];

	//		if (gp.p[0] >= -agrid && gp.p[0] < max_x
 //           && gp.p[1] >= -agrid && gp.p[1] < max_y
 //           && gp.p[2] >= -agrid && gp.p[2] < max_z) {
	//			lb_viscous_momentum_exchange(&gp,force) ;
	//		}

	//		//ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LB f = (%.6e,%.3e,%.3e)\n",this_node,p->f.f[0],p->f.f[1],p->f.f[2]));
	//	}
	//}
}

/** Coupling of a particle to viscous fluid with Stokesian friction.
* 
* Section II.C. Ahlrichs and Duenweg, JCP 111(17):8225 (1999)
*
* @param p          The coupled particle (Input).
* @param force      Coupling force between particle and fluid (Output).
*/
void LBD3Q19::lb_viscous_momentum_exchange(sPart *p, double force[3]) {

	int x,y,z;
	int node_index[8];
	double delta[6];
	double local_rho, *local_j, interpolated_u[3], delta_j[3];

	//ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: f = (%.3e,%.3e,%.3e)\n",this_node,p->f.f[0],p->f.f[1],p->f.f[2]));

	/* determine elementary lattice cell surrounding the particle 
	and the relative position of the particle in this cell */ 
	map_position_to_lattice(p->p,node_index,delta) ;

	//ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LB delta=(%.3f,%.3f,%.3f,%.3f,%.3f,%.3f) pos=(%.3f,%.3f,%.3f)\n",this_node,delta[0],delta[1],delta[2],delta[3],delta[4],delta[5],p->r.p[0],p->r.p[1],p->r.p[2]));

	/* calculate fluid velocity at particle's position
	this is done by linear interpolation
	(Eq. (11) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
	interpolated_u[0] = interpolated_u[1] = interpolated_u[2] = 0.0 ;
	for (z=0;z<2;z++) {
		for (y=0;y<2;y++) {
			for (x=0;x<2;x++) {

				//local_node = &lbfluid[node_index[(z*2+y)*2+x]];
				local_rho  = rho[node_index[(z*2+y)*2+x]];
				local_j    = j[node_index[(z*2+y)*2+x]];

				interpolated_u[0] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*local_j[0]/(local_rho);
				interpolated_u[1] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*local_j[1]/(local_rho);	  
				interpolated_u[2] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*local_j[2]/(local_rho);
				//cout << " " << interpolated_u[0]<< " "<< interpolated_u[1]<<" "<<interpolated_u[2]<<endl; 

			}
		}
	}

	//ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LB u = (%.16e,%.3e,%.3e) v = (%.16e,%.3e,%.3e)\n",this_node,interpolated_u[0],interpolated_u[1],interpolated_u[2],p->m.v[0],p->m.v[1],p->m.v[2]));

	/* calculate viscous force
	* take care to rescale velocities with time_step and transform to MD units 
	* (Eq. (9) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
	force[0] = - friction * (p->v[0]/time_step - interpolated_u[0]*agrid/tau);
	force[1] = - friction * (p->v[1]/time_step - interpolated_u[1]*agrid/tau);
	force[2] = - friction * (p->v[2]/time_step - interpolated_u[2]*agrid/tau);

	//ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LB f_drag = (%.6e,%.3e,%.3e)\n",this_node,force[0],force[1],force[2]));

	//ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LB f_random = (%.6e,%.3e,%.3e)\n",this_node,p->lc.f_random[try][0],p->lc.f_random[try][1],p->lc.f_random[try][2]));

	force[0] = force[0] + p->f_random[0];
	force[1] = force[1] + p->f_random[1];
	force[2] = force[2] + p->f_random[2];

	//ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LB f_tot = (%.6e,%.3e,%.3e) (try=%d)\n",this_node,force[0],force[1],force[2],try));

	/* transform momentum transfer to lattice units
	(Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
	delta_j[0] = - force[0]*time_step*tau/agrid;
	delta_j[1] = - force[1]*time_step*tau/agrid;
	delta_j[2] = - force[2]*time_step*tau/agrid;

	//cout << " force[0] " << force[0]<< " time_step "<<time_step<<" tau "<< tau << " agrid "<<agrid<<endl;
	lb_transfer_momentum(delta_j,node_index,delta);
}

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
void LBD3Q19::map_position_to_lattice(const double pos[3], int node_index[8], double delta[6]) {

	int dir, ind[3], box_l[3];
	double lpos, rel;

	int indexB;
	int boundary;

	box_l[0] = max_x;
	box_l[1] = max_y;
	box_l[2] = max_z;

	/* determine the elementary lattice cell containing the particle
	and the relative position of the particle in this cell */
	for (dir = 0; dir < 3; dir++) {

		lpos = pos[dir];
		rel = lpos / agrid + 1.0; // +1 for halo offset
		ind[dir] = (int) floor(rel);
		//cout << " lpos " << lpos << " rel " << rel <<" ind[dir] "<< ind[dir]<<endl;

		/* surrounding elementary cell is not completely inside this box,
		adjust if this is due to round off errors */
		if (ind[dir] < 0) {
			/*fprintf(
			stderr,
			"%d: map_position_to_lattice: position (%f,%f,%f) not inside a local plaquette.\n",
			this_node, pos[0], pos[1], pos[2]);*/
			throw ComException("map_position_to_lattice: position not inside a local plaquette.", __FILE__, __LINE__ );

		} else if (ind[dir] > box_l[dir]) {
			if (lpos - box_l[dir] < ROUND_ERROR_PREC)
				ind[dir] = box_l[dir] ;
			else
				/*fprintf(
				stderr,
				"%d: map_position_to_lattice: position (%f,%f,%f) not inside a local plaquette.\n",
				this_node, pos[0], pos[1], pos[2]);*/
				throw ComException("map_position_to_lattice: position not inside a local plaquette.", __FILE__, __LINE__ );
		}

		delta[3 + dir] = rel - ind[dir]; // delta_x/a
		delta[dir] = 1.0 - delta[3 + dir];

	}

	//cout <<" x " <<ind[0]<<" y " <<ind[1]<<" z " <<ind[2]<<" n "<<endl;
	//cout <<" dx " <<delta[0]<<" dy " <<delta[1]<<" dz " <<delta[2]<<" n "<<endl;
	//cout <<" dx " <<delta[3]<<" dy " <<delta[4]<<" dz " <<delta[5]<<" n "<<endl;
	indexB = (ind[0]+ imax_x * (ind[1]+ ind[2]*imax_y));
	node_index[0] = fb_index[indexB];//ind[0]+ fmax_x * (ind[1]+ ind[2]*fmax_y);
	//cout <<" node_index[0] " <<node_index[0]<<endl;

	//indexB = fb_index[(ind[0]+ imax_x * (ind[1]+ ind[2]*imax_y))];
	//cout <<" indexB " <<indexB<<endl;
	//if (indexB != node_index[0])throw ComException("node_index[0] != indexB.", __FILE__, __LINE__ );
	
	boundary = f_index[indexB][0];
	//printf("  boundary:%d \n",boundary);
	switch (boundary){
		case 0:
			node_index[1] = node_index[0] + 1;
			node_index[2] = node_index[0] + fmax_x;
			node_index[3] = node_index[0] + fmax_x + 1;
			node_index[4] = node_index[0] + fmax_x * fmax_y;
			node_index[5] = node_index[4] + 1;
			node_index[6] = node_index[4] + fmax_x;
			node_index[7] = node_index[4] + fmax_x + 1;
			break;
		case 1:
			node_index[1] = node_index[0] + 1 - max_x;
			node_index[2] = node_index[0] + fmax_x;
			node_index[3] = node_index[0] + fmax_x  + 1 - max_x;
			node_index[4] = node_index[0] + fmax_x * fmax_y;
			node_index[5] = node_index[4] + 1 - max_x;
			node_index[6] = node_index[4] + fmax_x;
			node_index[7] = node_index[4] + fmax_x  + 1 - max_x;
			break;
		case 2:
			node_index[1] = node_index[0] + 1;
			node_index[2] = node_index[0] + fmax_x*( 1 - max_y);
			node_index[3] = node_index[0] + fmax_x*( 1 - max_y) + 1;
			node_index[4] = node_index[0] + fmax_x * fmax_y;
			node_index[5] = node_index[4] + 1;
			node_index[6] = node_index[4] + fmax_x*( 1 - max_y);
			node_index[7] = node_index[4] + fmax_x*( 1 - max_y) + 1;
			break;
		case 4:
			node_index[1] = node_index[0] + 1;
			node_index[2] = node_index[0] + fmax_x;
			node_index[3] = node_index[0] + fmax_x + 1;
			node_index[4] = node_index[0] + fmax_x * fmax_y*( 1 - max_z);
			node_index[5] = node_index[4] + 1;
			node_index[6] = node_index[4] + fmax_x;
			node_index[7] = node_index[4] + fmax_x + 1;
			break;
		case 3:
			node_index[1] = node_index[0] + 1 - max_x;
			node_index[2] = node_index[0] + fmax_x*( 1 - max_y);
			node_index[3] = node_index[0] + fmax_x*( 1 - max_y) + 1 - max_x;
			node_index[4] = node_index[0] + fmax_x * fmax_y;
			node_index[5] = node_index[4] + 1 - max_x;
			node_index[6] = node_index[4] + fmax_x*( 1 - max_y);
			node_index[7] = node_index[4] + fmax_x*( 1 - max_y) + 1 - max_x;
			break;
		case 5:
			node_index[1] = node_index[0] + 1 - max_x;
			node_index[2] = node_index[0] + fmax_x;
			node_index[3] = node_index[0] + fmax_x + 1 - max_x;
			node_index[4] = node_index[0] + fmax_x * fmax_y*( 1 - max_z);
			node_index[5] = node_index[4] + 1 - max_x;
			node_index[6] = node_index[4] + fmax_x;
			node_index[7] = node_index[4] + fmax_x + 1 - max_x;
			break;
		case 6:
			node_index[1] = node_index[0] + 1;
			node_index[2] = node_index[0] + fmax_x*( 1 - max_y);
			node_index[3] = node_index[0] + fmax_x*( 1 - max_y) + 1;
			node_index[4] = node_index[0] + fmax_x * fmax_y*( 1 - max_z);
			node_index[5] = node_index[4] + 1;
			node_index[6] = node_index[4] + fmax_x*( 1 - max_y);
			node_index[7] = node_index[4] + fmax_x*( 1 - max_y) + 1;
			break;
		case 7:
			node_index[1] = node_index[0] + 1 - max_x;
			node_index[2] = node_index[0] + fmax_x*( 1 - max_y);
			node_index[3] = node_index[0] + fmax_x*( 1 - max_y) + 1 - max_x;
			node_index[4] = node_index[0] + fmax_x * fmax_y*( 1 - max_z);
			node_index[5] = node_index[4] + 1 - max_x;
			node_index[6] = node_index[4] + fmax_x*( 1 - max_y);
			node_index[7] = node_index[4] + fmax_x*( 1 - max_y) + 1 - max_x;
			break;
		default:
			throw ComException("Wrong boundary Type.", __FILE__, __LINE__ );
			return;
	}

	//node_index[1] = node_index[0] + 1;
	//node_index[2] = node_index[0] + fmax_x;
	//node_index[3] = node_index[0] + fmax_x + 1;
	//node_index[4] = node_index[0] + fmax_x * fmax_y;
	//node_index[5] = node_index[4] + 1;
	//node_index[6] = node_index[4] + fmax_x;
	//node_index[7] = node_index[4] + fmax_x + 1;
}

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
void LBD3Q19::lb_transfer_momentum(const double momentum[3], const int node_index[8], const double delta[6]) {

	int x, y, z, index;
	double *local_n, *n_tmp;
	double *local_j, delta_j[3];

	/* We don't need to save the local populations because 
	* we use a trick for their restoration:
	* We substract the old random force from the new one,
	* hence the previous change in the local populations
	* is automatically revoked during the recalculation.
	* Note that this makes it necessary to actually apply 
	* all changes and forbids to return immediately when negative
	* populations occur.
	*/

	for (z=0;z<2;z++) {
		for (y=0;y<2;y++) {
			for (x=0;x<2;x++) {

				index = node_index[(z*2+y)*2+x];
				n_tmp = f[index];
				local_n = n_tmp;
				local_j = j[index];

				delta_j[0] = delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*momentum[0];
				delta_j[1] = delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*momentum[1];
				delta_j[2] = delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*momentum[2];

				
				local_n[1] = n_tmp[1]  + 1./6.*delta_j[0];
				local_n[2] = n_tmp[2]  - 1./6.*delta_j[0];

				local_n[3]  = n_tmp[3]  + 1./6.*delta_j[1];
				local_n[4]  = n_tmp[4]  - 1./6.*delta_j[1];
				local_n[5]  = n_tmp[5]  + 1./6.*delta_j[2];
				local_n[6]  = n_tmp[6]  - 1./6.*delta_j[2];
				local_n[7]  = n_tmp[7]  + 1./12.*(delta_j[0]+delta_j[1]);
				local_n[8]  = n_tmp[8]  - 1./12.*(delta_j[0]+delta_j[1]);
				local_n[9]  = n_tmp[9]  + 1./12.*(delta_j[0]-delta_j[1]);
				local_n[10] = n_tmp[10] - 1./12.*(delta_j[0]-delta_j[1]);
				local_n[11] = n_tmp[11] + 1./12.*(delta_j[0]+delta_j[2]);
				local_n[12] = n_tmp[12] - 1./12.*(delta_j[0]+delta_j[2]);
				local_n[13] = n_tmp[13] + 1./12.*(delta_j[0]-delta_j[2]);
				local_n[14] = n_tmp[14] - 1./12.*(delta_j[0]-delta_j[2]);
				local_n[15] = n_tmp[15] + 1./12.*(delta_j[1]+delta_j[2]);
				local_n[16] = n_tmp[16] - 1./12.*(delta_j[1]+delta_j[2]);
				local_n[17] = n_tmp[17] + 1./12.*(delta_j[1]-delta_j[2]);
				local_n[18] = n_tmp[18] - 1./12.*(delta_j[1]-delta_j[2]);
			}
		}
	}

}

/* Get total mass*/
double LBD3Q19::lb_calc_fluid_mass()
{
	int index;
	double t_mass;
	int x, y, z;

	t_mass = 0.0;

	for (z=1; z<=max_z; z++) {
		for (y=1; y<=max_y; y++) {
			for (x=1; x<=max_x; x++) {
				index = (x + fmax_x*(y + fmax_y*z));
				//Comment it to reduce the calculation.
				//lb_calc_local_fields(index, false);
				t_mass += rho[index];
			}
		}
	}

	return t_mass;
}

/*  Output data.  */
void LBD3Q19::outTxt (int t)
{
	char	filename[1024];
	FILE	*file_ptr;

	int index;
	int i,m, k;
	int x, y, z;

	sprintf(filename,"%s%03d.dat" , outputFFileName, t) ;
	file_ptr = fopen(filename, "w");
	if(!file_ptr) throw ComException("Cann not open the output data of fliuds. ", __FILE__, __LINE__ );

	//fwrite (&(fvolume), sizeof(int), 1, file_ptr);
	for (x = 0; x < fmax_x; x++)
	{
		for (y = 0; y <  fmax_y; y++){
			for (z = 0; z <  fmax_z; z++){
				index = (x + fmax_x*(y + fmax_y*z));

				double local_rho = rho[index];
				double *local_j   = j[index];
				double *local_pi  = pi[index];

				//Comment it to reduce the calculation.
				//lb_calc_local_fields(index, true);

				fprintf (file_ptr, " %d %d %d ", x, y, z);

				//fwrite (&(rho[index]), sizeof(double), 1, file_ptr);
				//cout << *local_rho << " rho ";
				fprintf (file_ptr, "%f ", local_rho);
				m = 0;
				for (i=0;i<3;i++) {
					fprintf (file_ptr, "%f ", local_j[i]);
					//fwrite (&(local_j[i]), sizeof(double), 1, file_ptr);
					for (k=0;k<=i;k++) {
						fprintf (file_ptr, "%f ", local_pi[m]);
						//fwrite (&(local_pi[m]), sizeof(double), 1, file_ptr);
						m++;
					}
				}
				fprintf (file_ptr, "\n");
			}
		}
	}
	fclose (file_ptr);  file_ptr = 0;
}

/*  Output data.  */
void LBD3Q19::outN(int t)
{
	char	filename[1024];
	FILE	*file_ptr;

	int index;
	int i;
	int x, y, z;
	double *local_n;

	sprintf(filename,"%sN%06d.txt" , outputFFileName, t) ;
	file_ptr = fopen(filename, "w");
	if(!file_ptr) throw ComException("Cann not open the output data of fliuds. ", __FILE__, __LINE__ );

	fwrite (&(t), sizeof(int), 1, file_ptr);
	for (x = 0; x < fmax_x; x++)
	{
		for (y = 0; y <  fmax_y; y++){
			for (z = 0; z <  fmax_z; z++){
				index = (x + fmax_x*(y + fmax_y*z));
				local_n = f[index];

				//double *local_rho = rho[index];
				double *local_j   = j[index];
				double *local_pi  = pi[index];

				//Comment it to reduce the calculation.
				//lb_calc_local_fields(index,true);

				fprintf (file_ptr, " %d %d %d ", x, y, z);

				for (i=0;i<19;i++) {
					fprintf (file_ptr, "%f ", local_n[i]);
					//cout << " " << local_n[i];
				}
				fprintf (file_ptr, "\n");
			}
		}
	}
	fclose (file_ptr);  file_ptr = 0;
}

/* Get total movement*/
void LBD3Q19::lb_calc_fluid_momentum(double* momentum)
{
	int index;
	int x, y, z;
	double* local_j;

	momentum[0] = 0.0;momentum[1] = 0.0;momentum[2] = 0.0;

	for (z=1; z<=max_z; z++) {
		for (y=1; y<=max_y; y++) {
			for (x=1; x<=max_x; x++) {
				index = (x + fmax_x*(y + fmax_y*z));
				local_j   = j[index];
				//Comment it to reduce the calculation.
				//lb_calc_local_fields(index, false);
				momentum[0] += local_j[0];
				momentum[1] += local_j[1];
				momentum[2] += local_j[2];
			}
		}
	}
	momentum[0] *= agrid/tau;
    momentum[1] *= agrid/tau;
    momentum[2] *= agrid/tau;
}

/** Calculate temperature of the LB fluid.
 * \param result Fluid temperature
 */
double LBD3Q19::lb_calc_fluid_temp() {
  int x, y, z, index;
  double local_rho, local_j2;
  double temp = 0.0;

  for (x=1; x<=max_x; x++) {
    for (y=1; y<=max_y; y++) {
      for (z=1; z<=max_z; z++) {
		  index = (x + fmax_x*(y + fmax_y*z));

		  //Comment it to reduce the calculation.
		  //lb_calc_local_fields(index, false);
		  //lb_calc_local_j(&lbfluid[index]);
		  //lb_calc_local_rho(&lbfluid[index]);

		  local_rho = rho[index];
		  local_j2  = scalar(j[index],j[index]);

		  temp += local_j2;
      }
    }
  }
  temp *= 1./(total_rho*volume*tau*tau*pow(agrid,4));
  //cout << "  temp "<<temp<<" total_rho "<<total_rho<<" volume "<<volume<<" tau "<<tau<<" pow(agrid,4) "<<pow(agrid,4)<<endl;
  return temp;
}

/*  LB: output  properties*/
void LBD3Q19::outProp(int t, FILE	*file_ptr)
{
	void   *io_ptr;

	int index;
	int i,m, k;


	//fwrite (&(fvolume), sizeof(int), 1, file_ptr);
	//cout << lblattice->halo_grid_volume << " lblattice->halo_grid_volume "<<endl;

	io_ptr = (void *) &t;
	fwrite (io_ptr, sizeof(int), 1, file_ptr);

	for (index=0;index<fvolume;index++) {

		//double *local_rho = rho[index];
		double *local_j   = j[index];
		double *local_pi  = pi[index];

		//Comment it to reduce the calculation.
		//lb_calc_local_fields(index, true);

		fwrite (&(rho[index]), sizeof(double), 1, file_ptr);
		//cout << *local_rho << " rho ";
		m = 0;
		for (i=0;i<3;i++) {
			fwrite (&(local_j[i]), sizeof(double), 1, file_ptr);
			//if (index < 2200) cout << " "<<local_j[i] << " j" << i << " ";
			for (k=0;k<=i;k++) {
				fwrite (&(local_pi[m]), sizeof(double), 1, file_ptr);
				//cout << "local_pi[m]" << local_pi[m] << " m " << m << " ";
				m++;
			}
		}
		//if (index == 32) cout << index << " index " << endl;
	}
	//cout << endl;
}

/*  Read properties file  */
void LBD3Q19::inProp(int t, FILE *file_ptr)
{
	int index; double prho; double pj[3]; double ppi[6];
	int m, i, k;
	int dt;

	if(file_ptr != 0) 
	{
		//cout << "open " << filename << " to read lbm details." << endl;

		fread ( &(dt), sizeof(int), 1, file_ptr);
		if (t!=dt)throw ComException("time in the fluid properties file does not match.", __FILE__, __LINE__ );

		//cout << lblattice->halo_grid_volume << " lblattice->halo_grid_volume "<<endl;
		for (index=0;index<fvolume;index++) {

			fread (&prho, sizeof(double), 1, file_ptr);
			//cout << prho << " rho ";
			m = 0;
			for (i=0;i<3;i++) {
				fread (&(pj[i]), sizeof(double), 1, file_ptr);
				//cout << pj[i]<<" j"<< i <<" ";
				for (k=0;k<=i;k++) {
					//io_ptr = (void *) pi;
					fread (&(ppi[m]), sizeof(double), 1, file_ptr);
					//cout << ppi[m]<<" pi"<<m <<" ";
					m++;
				}
			}
			//cout <<" index "<< index;
			//lb_set_local_fields(index,prho,pj,ppi);
		}
		//cout <<endl;
	}
}