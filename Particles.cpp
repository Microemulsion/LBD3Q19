// Particles.cpp: implementation of the Particles class.
//
//////////////////////////////////////////////////////////////////////

#include "Particles.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

sPart::sPart(void)
{
	identity = -1;
	type     = 0;
	mol_id   = -1;
	mass	= 1.0;
	//#endif

	p[0] = 0.0; p[1] = 0.0; p[2] = 0.0;

	v[0] = 0.0; v[1] = 0.0; v[2] = 0.0;

	f[0] = 0.0; f[1] = 0.0; f[2] = 0.0;

	p_old[0] = 0.0; p_old[1] = 0.0; p_old[2] = 0.0;
	/** index of the simulation box image where the particle really sits. */
	i[0] = 0; i[1] = 0; i[2] = 0;

	ext_flag = 0;
	/** External force, apply if \ref ParticleLocal::ext_flag == 1. */
	ext_force[0] = 0.0; ext_force[01] = 0.0; ext_force[2] = 0.0;
	//#endif

	f_random[0] = 0.0; f_random[1] = 0.0; f_random[2] = 0.0;

}

sPart::~sPart(void)
{
}

Cell::Cell(void)
{
	n_neighbors = 14;
	p_contGhosts = false;
	gPMove[0] = 0.0;
	gPMove[1] = 0.0;
	gPMove[2] = 0.0;
	noGCell = -1;
}

Cell::~Cell(void)
{
}

/* Add a particle pair to a verlet pair list. */
void Cell::add_pair(int n, sPart* p1, sPart* p2){
	pair< sPart* , sPart* > datax;
	datax.first  = p1; datax.second = p2;
	cell_inter[n].vList.push_back(datax);
}

void Cell::clear_pairs(int n){
	cell_inter[n].vList.erase(cell_inter[n].vList.begin(), cell_inter[n].vList.end());
}

Particles::Particles(const BaseMesh& pmesh, double max_cut, double skin, double time_step, char* outputPFileName)
:p_mesh( pmesh )
,p_outputPFileName( outputPFileName )
,p_timestep( time_step )
{
	int i;
	double box_l[3];

	max_num_cells = CELLS_MAX_NUM_CELLS;
	min_num_cells = 1;

	/** The number of nodes in each spatial dimension. */
	node_grid[0] = 1;
	node_grid[1] = 1;
	node_grid[2] = 1;
	/** position of node in node grid */
	node_pos[0] = 0;
	node_pos[1] = 0;
	node_pos[2] = 0;

	box_l[0] = p_mesh.getNx();
	box_l[1] = p_mesh.getNy();
	box_l[2] = p_mesh.getNz();

	integrate_vv_recalc_maxrange(max_cut, skin);

	for (i=0; i<3; i++){

		/** Dimensions of the box a single node is responsible for. */ 
		local_box_l[i] = box_l[i]/(double)node_grid[i];

		my_left[i]   = node_pos[i]    *local_box_l[i];
		my_right[i]  = (node_pos[i]+1)*local_box_l[i];
		//cout << local_box_l[i]<< " "<<cell_size[i]<<" "<<my_left[i]<<" "<<my_right[i]<<endl;
	}
	dd_create_cell_grid(max_cut);//for cell_grid

	for (i=0; i<3; i++){
		cell_size[i] = local_box_l[i]/(double)cell_grid[i];
		//cout << local_box_l[i]<< " "<<cell_size[i]<<" "<<my_left[i]<<" "<<my_right[i]<<endl;
	}


	n_total_particles = 0;

	n_cells = 0;

	use_vList = false;
}

Particles::~Particles(void)
{
	unsigned int i;

	for (i=0; i<local_particles.size(); i++ )
		delete local_particles[i];

	n_total_particles = 0;

	n_cells = 0;
	for (i=0; i<local_cells.size(); i++ )
		delete local_cells[i];

	for (i=0; i<ghost_cells.size(); i++ )
		delete ghost_cells[i];
}

// Get the detail of particles from input file.
void Particles::readInputF(char* inputFileName)
{
	int identity, type;
	double mass;
	double p[3], f[3], v[3];

	string txt;
	sPart* temp_particle;

	if (inputFileName == NULL) throw ComException("Input file name is NULL. ", __FILE__, __LINE__ );
	ifstream fp(inputFileName);
	if (fp.fail()) throw ComException("Can not open input file. ", __FILE__, __LINE__ );

	getline(fp,txt);
	while (txt.find("Object parameters")) {
		getline(fp,txt);
		if (fp.eof()) throw ComException("Object parameters not found in input file. ", __FILE__, __LINE__ );
	}

	while (txt.find("OBJECTDETAILS")) {
		getline(fp,txt);
		if (fp.eof()) throw ComException("Object details not found in input file. ", __FILE__, __LINE__ );
	}

	//cout <<"\n Object details: " << endl;

	n_total_particles = 0;
	fp	>> identity >> type >> mass >> p[0] >> p[1] >> p[2] >> v[0] >> v[1] >> v[2] >> f[0] >> f[1] >> f[2];
	while (!fp.fail())
	{


		//cout <<"  "<< identity <<"  "<< type <<"  "<< mass <<"  "<< p[0] <<"  "<< p[1] <<"  "<< p[2] <<"  "<< f[0] <<"  "<< f[1] <<"  "<< f[2] <<"  "<< v[0] <<"  "<< v[1] <<"  "<< v[2] << endl;

		if (identity!=(n_total_particles)) break;

		place_particle(identity, p, true);

		temp_particle = local_particles[identity];
		temp_particle->mass = mass;
		temp_particle->type = type;
		memcpy(temp_particle->v, v, 3*sizeof(double));
		memcpy(temp_particle->f, f, 3*sizeof(double));
		fp	>> identity >> type >> mass >> p[0] >> p[1] >> p[2] >> v[0] >> v[1] >> v[2]>> f[0] >> f[1] >> f[2];
	}
	fp.close();
}

void Particles::place_particle(int part, double p[3], bool bnew)
{
	Cell* pcell;
	int i[3];
	sPart* temp_particle;

	if (part < 0) throw ComException("Particle identity < 0. ", __FILE__, __LINE__ );

	i[0] = 0;
	i[1] = 0;
	i[2] = 0;
	fold_position(p_mesh, p, i);

	if (bnew) {

		temp_particle = add_one_particle();

		/* allocate particle anew */
		pcell = position_to_cell(p);

		//cout << " pcell "<< pcell->n_neighbors <<endl;
		if (!pcell) {
			//fprintf(stderr, "%d: INTERNAL ERROR: particle %d at %f(%f) %f(%f) %f(%f) does not belong on this node\n",
			//   this_node, part, p[0], pp[0], p[1], pp[1], p[2], pp[2]);
			//errexit();
			throw ComException("INTERNAL ERROR: particle does not belong on this node\n. ", __FILE__, __LINE__ );
		}
		pcell->p_objects.push_back(temp_particle);

		//cout << " size "<<pcell->p_objects.size();

		temp_particle->identity = part;
		//cout << " cell indenti " <<local_cells[0]->p_objects[pcell->p_objects.size()-1]->identity<<endl;
	}
	else
		temp_particle = local_particles[part];

	///*PART_TRACE(fprintf(stderr, "%d: local_place_particle: got particle id=%d @ %f %f %f\n",
	//this_node, part, p[0], p[1], p[2]));*/

	memcpy(temp_particle->p, p, 3*sizeof(double));
	memcpy(temp_particle->i, i, 3*sizeof(int));
}

sPart* Particles::add_one_particle()
{
	sPart* tempsPart  = new sPart();

	n_total_particles++;

	local_particles.push_back(tempsPart);

	return tempsPart;
}

/** fold particle coordinates to primary simulation box.
\param pos the position...
\param image_box and the box
Both pos and image_box are I/O,
i. e. a previously folded position will be folded correctly.*/
void Particles::fold_position(const BaseMesh& p_mesh, double pos[3],int image_box[3])
{
	int i;
	double box_l_i[3], box_l[3];
	box_l[0] = p_mesh.getNx();
	box_l[1] = p_mesh.getNy();
	box_l[2] = p_mesh.getNz();

	box_l_i[0] = 1.0/box_l[0];
	box_l_i[1] = 1.0/box_l[1];
	box_l_i[2] = 1.0/box_l[2];
	for(i=0;i<3;i++)
	{	/** fold a coordinate to primary simulation box.
		\param pos         the position...
		\param image_box   and the box
		\param dir         the coordinate to fold: dir = 0,1,2 for x, y and z coordinate.
		Both pos and image_box are I/O,
		i. e. a previously folded position will be folded correctly.*/

		int tmp;

		image_box[i] += (tmp = (int)floor(pos[i]*box_l_i[i]));
		pos[i]        = pos[i] - tmp*box_l[i];    
		if(pos[i] < 0 || pos[i] >= box_l[i]) {
			/* slow but safe */
			if (fabs(pos[i]*box_l_i[i]) >= INT_MAX/2) {
				//char *errtext = runtime_error(128 + TCL_INTEGER_SPACE);
				//ERROR_SPRINTF(errtext,"{086 particle coordinate out of range, pos = %f, image box = %d} ", pos[i], image_box[i]);
				throw ComException("Iparticle coordinate out of range.\n. ", __FILE__, __LINE__ );
				image_box[i] = 0;
				pos[i] = 0;
				return;
			}
		}
	}
}

/** unfold coordinates to physical position.
\param pos the position...
\param image_box and the box

Both pos and image_box are I/O, i.e. image_box will be (0,0,0)
afterwards.
*/
void Particles::unfold_position(double pos[3],int image_box[3])
{
	double box_l[3];
	box_l[0] = p_mesh.getNx();
	box_l[1] = p_mesh.getNy();
	box_l[2] = p_mesh.getNz();

	int i;
	for(i=0;i<3;i++) {
		pos[i]       = pos[i] + image_box[i]*box_l[i];    
		image_box[i] = 0;
	}
}

Cell* Particles::position_to_cell(double pos[3])
{
	int i,cpos[3];
	double lpos;

	for(i=0;i<3;i++) {
		lpos = pos[i] - my_left[i];

		//cpos[i] = (int)(lpos/cell_size[i])+1;
		cpos[i] = (int)(lpos/cell_size[i]);

		//kxsong use normal grid.
		//if (cpos[i] < 1) {
		//	cpos[i] = 1;
		//	if (lpos > -ROUND_ERROR_PREC) cpos[i] = 1;
		//	else	return NULL;
		//}
		//else if (cpos[i] > cell_grid[i]) {
		//	if (lpos < local_box_l[i] + ROUND_ERROR_PREC) cpos[i] = cell_grid[i];
		//	else	return NULL;
		//}
	}
	//cout << " cpos[0] "<< cpos[0]<< " cpos[1] " <<cpos[1]<< " cpos[2] " <<cpos[2]<<endl;
	//cout << " ghost_cell_grid "<< ghost_cell_grid[0]<< " cpos[1] " <<ghost_cell_grid[1]<< " cpos[2] " <<ghost_cell_grid[2]<<endl;
	i = get_linear_index(cpos[0],cpos[1],cpos[2], cell_grid);  
	//cout << " i "<< i << " local_cells.size " << local_cells.size() << endl;
	return local_cells[i];
}

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
void Particles::dd_create_cell_grid(double max_cut)
{
	int i, j, k, n_local_cells,new_cells,min_ind;
	double cell_range[3], min_box_l, min_size, scale, volume;
	bool is_ghost;
	int n_ghostcells, n_localcells;
	//CELL_TRACE(fprintf(stderr, "%d: dd_create_cell_grid: max_range %f\n",this_node,max_range));
	//CELL_TRACE(fprintf(stderr, "%d: dd_create_cell_grid: local_box %f-%f, %f-%f, %f-%f,\n",this_node,my_left[0],my_right[0],my_left[1],my_right[1],my_left[2],my_right[2]));

	/* initialize */
	min_box_l = dmin(dmin(local_box_l[0],local_box_l[1]),local_box_l[2]);
	cell_range[0]=cell_range[1]=cell_range[2] = max_range;

	if(max_range2 < 0.0) {
		/* this is the initialization case */
		n_local_cells = cell_grid[0] = cell_grid[1] = cell_grid[2]=1;
	}
	else {
		/* Calculate initial cell grid */
		volume = local_box_l[0];
		for(i=1;i<3;i++) volume *= local_box_l[i];
		scale = pow(max_num_cells/volume, 1./3.);
		for(i=0;i<3;i++) {
			/* this is at least 1 */
			cell_grid[i] = (int)ceil(local_box_l[i]*scale);
			cell_range[i] = local_box_l[i]/cell_grid[i];
			//cout << "local_box_l[i] "<< local_box_l[i] << " scale "<< scale<<endl; 
			//cout << " cell_grid[i] "<< cell_grid[i]<< " cell_range[i] "<< cell_range[i]<<endl; 
			if ( cell_range[i] < max_range ) {
				/* ok, too many cells for this direction, set to minimum */
				cell_grid[i] = (int)floor(local_box_l[i]/max_range);
				if ( cell_grid[i] < 1 ) {
					//char *error_msg = runtime_error(TCL_INTEGER_SPACE + 2*TCL_DOUBLE_SPACE + 128);
					//ERROR_SPRINTF(error_msg, "{002 interaction range %f in direction %d is larger than the local box size %f} ",
					//		max_range, i, local_box_l[i]);
					throw ComException("002 interaction range in direction is larger than the local box size.\n. ", __FILE__, __LINE__ );
					cell_grid[i] = 1;
				}
				cell_range[i] = local_box_l[i]/cell_grid[i];
			}
		}

		/* It may be necessary to asymmetrically assign the scaling to the coordinates, which the above approach will not do.
		For a symmetric box, it gives a symmetric result. Here we correct that. */
		for (;;) {
			n_local_cells = cell_grid[0];
			for (i = 1; i < 3; i++)
				n_local_cells *= cell_grid[i];

			/* done */
			if (n_local_cells <= max_num_cells)
				break;

			/* find coordinate with the smallest cell range */
			min_ind = 0;
			min_size = cell_range[0];
			for (i = 1; i < 3; i++)
				if (cell_grid[i] > 1 && cell_range[i] < min_size) {
					min_ind = i;
					min_size = cell_range[i];
				}
				//     CELL_TRACE(fprintf(stderr, "%d: minimal coordinate %d, size %f, grid %d\n", this_node,min_ind, min_size, dd.cell_grid[min_ind]));

				cell_grid[min_ind]--;
				cell_range[min_ind] = local_box_l[min_ind]/cell_grid[min_ind];
		}
		//   CELL_TRACE(fprintf(stderr, "%d: final %d %d %d\n", this_node, dd.cell_grid[0], dd.cell_grid[1], dd.cell_grid[2]));

		/* sanity check */
		if (n_local_cells < min_num_cells) {
			//     char *error_msg = runtime_error(TCL_INTEGER_SPACE + 2*TCL_DOUBLE_SPACE + 128);
			//     ERROR_SPRINTF(error_msg, "{001 number of cells %d is smaller than minimum %d (interaction range too large or max_num_cells too small)} ",
			//	    n_local_cells, min_num_cells);
			throw ComException("001 number of cells is smaller than minimum (interaction range too large or max_num_cells too small).\n. ", __FILE__, __LINE__ );
		}
	}

	/* quit program if unsuccesful */
	if(n_local_cells > max_num_cells) {
		//   char *error_msg = runtime_error(128);
		//   ERROR_SPRINTF(error_msg, "{003 no suitable cell grid found} ");
		throw ComException("003 no suitable cell grid found.\n. ", __FILE__, __LINE__ );
	}

	/* now set all dependent variables */
	new_cells=1;
	for(i=0;i<3;i++) {
		ghost_cell_grid[i] = cell_grid[i]+2;	
		new_cells              *= ghost_cell_grid[i];
		cell_size[i]       = local_box_l[i]/(double)cell_grid[i];
		//   dd.inv_cell_size[i]   = 1.0 / dd.cell_size[i];
	}
	cell_range[0] = dmin(dmin(cell_size[0],cell_size[1]),cell_size[2]);
	if (max_cut >= 0.0)
		max_skin = cell_range[0] - max_cut;
	else
		max_skin = 0.0;

	/* allocate cell array and cell pointer arrays */
	//cout << " n_local_cells " << n_local_cells<<endl;
	realloc_cells(local_cells, n_local_cells);
	n_cells = n_local_cells;
	//cout << " (new_cells-n_local_cells) "<<(new_cells-n_local_cells)<< " new_cells "<<new_cells<<endl;
	realloc_cells(ghost_cells, (new_cells-n_local_cells));
	for(i=0; i<(new_cells-n_local_cells); i++) {
		ghost_cells[i]->n_neighbors = 0;
		ghost_cells[i]->p_contGhosts = true;
	}
	n_cells = new_cells;

	//kxsong creat total cells.
	is_ghost = false;
	n_ghostcells = 0;
	n_localcells = 0;
	for(k=0; k<ghost_cell_grid[2]; k++)
	{
		for(j=0; j<ghost_cell_grid[1]; j++)
		{
			for(i=0; i<ghost_cell_grid[0]; i++)
			{
				if (i == 0 ||i == (ghost_cell_grid[0]-1)) is_ghost = true;
				if (j == 0 ||j == (ghost_cell_grid[1]-1)) is_ghost = true;
				if (k == 0 ||k == (ghost_cell_grid[2]-1)) is_ghost = true;
				//cout << " i "<<i<<" j "<<j<<" k "<< k<<" index "<< get_linear_index(i,j,k,ghost_cell_grid)<<endl;
				if (is_ghost) { 
					//cout << " n_ghostcells "<< n_ghostcells<<endl;
					cells.push_back(ghost_cells[n_ghostcells]);
					//cout << " n_ghostcells "<< n_ghostcells<<endl;
					n_ghostcells++;

					is_ghost = false;
				} else {
					//cout << " n_localcells "<< n_localcells<<endl;
					cells.push_back(local_cells[n_localcells]);
					//cout << " n_localcells "<< n_localcells<<endl;
					n_localcells++;

				}
			}
		}
	}

	// CELL_TRACE(fprintf(stderr, "%d: dd_create_cell_grid, n_cells=%d, local_cells.n=%d, ghost_cells.n=%d, dd.ghost_cell_grid=(%d,%d,%d)\n", this_node, n_cells,local_cells.n,ghost_cells.n,dd.ghost_cell_grid[0],dd.ghost_cell_grid[1],dd.ghost_cell_grid[2]));
}

void Particles::integrate_vv_recalc_maxrange(double max_cut, double skin)
{
	//INTEG_TRACE(fprintf(stderr,"%d: integrate_vv_recalc_maxrange:\n",this_node));

	/* maximal interaction cutoff */
	//calc_maximal_cutoff();

	if (max_cut <= 0.0) {
		max_range  = -1.0;
		max_range2 = -1.0;
		return;
	}
	max_range            = max_cut;

	/* at beginning be nice */
	if (skin > 0.0) {
		max_range            += skin;
	}

	max_range2            = SQR(max_range);
}

void Particles::realloc_cells(vector<Cell*>& icells, int size)
{
	int i;
	//CELL_TRACE(fprintf(stderr, "%d: realloc_cells %d\n", this_node, size));
	int old_size = icells.size();

	if (size < old_size){
		/* free all memory associated with cells to be deleted. */
		for(i=size; i<old_size; i++) {
			delete icells[i];
		}
		icells.resize(size);
	}

	if (size > old_size){

		/* initialize new cells */
		for(i=old_size; i<size; i++) {
			Cell* tempCell  = new Cell();
			icells.push_back(tempCell);
		}
	}
	//cout << "icells.size() "<<icells.size()<< endl;
}

/*  Read checkpoint file  */
void Particles::readOutputF(int t)
{
	char	filename[1024];
	FILE	*file_ptr;
	void   *io_ptr;
	int n_sph;

	sPart data;

	sprintf(filename,"%s%03d.dat" , p_outputPFileName, t) ;
	file_ptr = fopen (filename, "rb");

	if (file_ptr != 0)
	{
		cout << "open " << filename << " to read particle details." << endl;

		io_ptr = (void *) &t;
		fread (io_ptr, sizeof(int), 1, file_ptr);

		for (n_sph = 0; n_sph < n_total_particles; n_sph++)
		{
			io_ptr = (void *)    &data.p;
			fread (io_ptr, sizeof(data.p), 1, file_ptr);
			io_ptr = (void *)    &data.v;
			fread (io_ptr, sizeof(data.v), 1, file_ptr);
			io_ptr = (void *)    &data.f;
			fread (io_ptr, sizeof(data.f), 1, file_ptr);
			io_ptr = (void *)    &data.mass;
			fread (io_ptr, sizeof(data.mass), 1, file_ptr);
			io_ptr = (void *)    &data.type;
			fread (io_ptr, sizeof(data.type), 1, file_ptr);
			io_ptr = (void *)     &data.identity;
			fread (io_ptr, sizeof(data.identity), 1, file_ptr);
			//cout << " ind "<<data.identity<<endl;

			//place_particle(data.identity, data.p, 1);
			sPart* temp_particle = local_particles[n_sph];
			temp_particle->identity = data.identity;
			temp_particle->mass = data.mass;
			temp_particle->type = data.type;
			memcpy(temp_particle->p, data.p, 3*sizeof(double));
			memcpy(temp_particle->v, data.v, 3*sizeof(double));
			memcpy(temp_particle->f, data.f, 3*sizeof(double));
		}
		fclose (file_ptr);  file_ptr = 0;
	}
}

/*  Particles: outPut  */
void Particles::outPut (int t)
{
	char	filename[1024];
	FILE	*file_ptr;
	int n_sph;
	void   *io_ptr;

	sprintf(filename,"%s%03d.dat" , p_outputPFileName, t) ;
	file_ptr = fopen(filename, "wb");
	if(!file_ptr) throw ComException("Can not open the output data of particles. ", __FILE__, __LINE__ );

	io_ptr = (void *) &t;
	fwrite (io_ptr, sizeof(int), 1, file_ptr);
	for (n_sph = 0; n_sph < n_total_particles; n_sph++)
	{
		sPart* data = local_particles[n_sph];

		if ( data!= NULL) {
			unfold_position(data->p, data->i);
			io_ptr = (void *)     &data->p;
			fwrite (io_ptr, sizeof(data->p), 1, file_ptr);
			io_ptr = (void *)     &data->v;
			fwrite (io_ptr, sizeof(data->v), 1, file_ptr);
			io_ptr = (void *)     &data->f;
			fwrite (io_ptr, sizeof(data->f), 1, file_ptr);
			io_ptr = (void *)     &data->mass;
			fwrite (io_ptr, sizeof(data->mass), 1, file_ptr);
			io_ptr = (void *)     &data->type;
			fwrite (io_ptr, sizeof(data->type), 1, file_ptr);
			io_ptr = (void *)     &data->identity;
			fwrite (io_ptr, sizeof(data->identity), 1, file_ptr);
			//cout<<" data->p "<<data->p[0] << " data->i "<<data->identity<<" data->v "<<data->v[0]<<" data->f "<< data->f[0]<< " data->TYPE "<< data->type <<" data->mass " <<data->mass<<endl;
		}
	}
	fclose (file_ptr);  file_ptr = 0;
}

void Particles::get_mi_vector(double res[3], double a[3], double b[3])
{
	int i;

	for(i=0;i<3;i++) {
		res[i] = a[i] - b[i];
		res[i] -= dround(res[i]/local_box_l[i])*local_box_l[i];
	}
}

void Particles::cells_update_ghosts()
{

	/* methods using skin rsp. rebuild_verletlist */
	if(use_vList) {
		if (rebuild_verletlist)
			/* Communication step:  number of ghosts and ghost information */
			cells_resort_particles(CELL_NEIGHBOR_EXCHANGE);
		else
			/* Communication step: ghost information */
			;//update_ghost_pos_comm.go();//ghost_communicator(&cell_structure.update_ghost_pos_comm);
	}
	else 
		cells_resort_particles(CELL_NEIGHBOR_EXCHANGE);
}

/*************************************************/

void Particles::cells_resort_particles(int global_flag)
{
	//CELL_TRACE(fprintf(stderr, "%d: entering cells_resort_particles %d\n", this_node, global_flag));

	//invalidate_ghosts();

	//particle_invalidate_part_node();
	//n_verlet_updates++;
	//cout << " 203 "<<local_cells[203]->p_objects.size()<<endl;
	dd_exchange_and_sort_particles(global_flag);
	//cout << " 203 "<<local_cells[203]->p_objects.size()<<endl;
	//update_ghosts();

	//cout << " 203 "<<local_cells[203]->p_objects.size()<<endl;
	//ghost_cells_comm.go();//ghost_communicator(&cell_structure.ghost_cells_comm);

	//exchange_ghosts_comm.go();//ghost_communicator(&cell_structure.exchange_ghosts_comm);

	//on_resort_particles();

	rebuild_verletlist = true;

	//recalc_forces = 1;

	//CELL_TRACE(fprintf(stderr, "%d: leaving cells_resort_particles\n", this_node));
}

void Particles::invalidate_ghosts()
{
	sPart *part;
	int c, np, p;

	vector<sPart*>::iterator first;
	vector<Cell*>::iterator last;
	int gcellsize = ghost_cells.size();

	/* remove ghosts, but keep Real Particles */
	for(c=gcellsize; c>0; c--) {
		//part = ghost_cells.cell[c]->part;

		np   = ghost_cells[c-1]->p_objects.size();
		for(p=0 ; p<np; p++) {
			//cout << " np " << np << " c "<< c <<endl; 
			/* Particle is stored as ghost in the local_particles array,
			if the pointer stored there belongs to a ghost celll
			particle array. */
			part   = ghost_cells[c-1]->p_objects[p];
			if( part == local_particles[part->identity] ) 
			{
				first = local_particles.begin();
				first += part->identity;
				local_particles.erase(first);
			}
		}
		ghost_cells.pop_back();
	}
}

/************************************************************/
void  Particles::dd_exchange_and_sort_particles(int global_flag)
{
	int dir, c, p, finished=0;
	Cell *sort_cell;
	Cell* cell;
	sPart *part;
	int np, p_index;
	//CELL_TRACE(fprintf(stderr,"%d: dd_exchange_and_sort_particles(%d):\n",this_node,global_flag));

	int lcellsize = local_cells.size();

	//cout << " lcellsize "<< lcellsize<<endl;

	while(finished == 0 ) {
		finished=1;
		/* direction loop: x, y, z */  
		for(dir=0; dir<3; dir++) { 
			/* Single node direction case (no communication) */
			/* Fold particles that have left the box */
			/* particle loop */
			for(c=0; c<lcellsize; c++) {
				cell = local_cells[c];
				np = cell->p_objects.size();
				p_index = 0;
				for (p = 0; p < np; p++) {

					part = cell->p_objects[p-p_index];
					//cout << " i "<<part->identity<<" p0 "<< part->p[0]<< " p1 "<< part->p[1]<< " p2 "<< part->p[2]<< " i0 "<< part->i[0]<< " i1 "<< part->i[1]<< " i2 "<< part->i[2]<<endl;
					fold_coordinate(part->p, part->i, dir);

					if (dir==2) {
						//cout << " i ";
						sort_cell = position_to_cell(part->p);

						if(sort_cell != cell) {
							//cout << " part "<<part->identity<<" c "<<c<<endl;

							if(sort_cell==NULL) {
								cout << " !!!part "<<part->identity<<endl;
								//CELL_TRACE(fprintf(stderr, "%d: dd_exchange_and_sort_particles: CP2 Particle %d (%f,%f,%f) not inside node domain.\n",
								//this_node,part->p.identity,part->r.p[0],part->r.p[1],part->r.p[2]));
								finished=0;
								sort_cell = local_cells[0];
								if(sort_cell != cell) {

									move_indexed_particle(sort_cell, cell, p-p_index);
									if(p < np) p_index++;
								}      
							}
							else {
								//CELL_TRACE(fprintf(stderr, "%d: dd_exchange_and_sort_particles: move particle id %d\n", this_node,part->p.identity));

								move_indexed_particle(sort_cell, cell, p-p_index);
								if(p < np) p_index++;
							}
						}
					}
				}
			}
		}
		//CELL_TRACE(fprintf(stderr,"%d: dd_exchange_and_sort_particles: finished value: %d\n",this_node,finished));
	}

	//CELL_TRACE(fprintf(stderr,"%d: dd_exchange_and_sort_particles finished\n",this_node));
}

/** fold a coordinate to primary simulation box.
\param pos         the position...
\param image_box   and the box
\param dir         the coordinate to fold: dir = 0,1,2 for x, y and z coordinate.

Both pos and image_box are I/O,
i. e. a previously folded position will be folded correctly.
*/
void Particles::fold_coordinate(double pos[3], int image_box[3], int dir)
{
	int tmp;
	image_box[dir] += (tmp = (int)floor(pos[dir]/local_box_l[dir]));
	pos[dir]        = pos[dir] - tmp*local_box_l[dir];
	if(pos[dir] < 0 || pos[dir] >= local_box_l[dir]) {
		/* slow but safe */
		if (fabs(pos[dir]/local_box_l[dir]) >= INT_MAX/2) {
			//char *errtext = runtime_error(128 + TCL_INTEGER_SPACE);
			//ERROR_SPRINTF(errtext,"{086 particle coordinate out of range, pos = %f, image box = %d} ", pos[dir], image_box[dir]);
			throw ComException("086 particle coordinate out of range. ", __FILE__, __LINE__ );
			image_box[dir] = 0;
			pos[dir] = 0;
			return;
		}
	}
}

void Particles::move_indexed_particle(Cell *dl, Cell *sl, int i)
{
	sPart *dst;
	sPart *src;
	//cout <<" sl0 "<<sl->p_objects.size()<<" dl0 "<<dl->p_objects.size();
	//cout << " 203 "<<local_cells[203]->p_objects.size()<<endl;
	src = sl->p_objects[i];
	dl->p_objects.push_back(src);
	//vector<sPart*>::iterator first;
	//first = sl->p_objects.begin();
	//first += i;

	//sl->p_objects.erase(first);
	int cSize = sl->p_objects.size();
	if (i<(cSize-1)) sl->p_objects[i]=sl->p_objects[cSize-1];
	sl->p_objects.pop_back();
	//cout <<" sl1 "<<sl->p_objects.size()<<" dl1 "<<dl->p_objects.size();
	//cout << " 203 "<<local_cells[203]->p_objects.size()<<endl;

}

void Particles::cells_pre_init()
{
	// CellPList tmp_local;
	// CELL_TRACE(fprintf(stderr, "%d: cells_pre_init\n",this_node));
	// /* her local_cells has to be a NULL pointer */
	// if(local_cells.cell != NULL) {
	//   fprintf(stderr,"INTERNAL ERROR: wrong usage of cells_pre_init!\n");
	//   errexit();
	//throw ComException("INTERNAL ERROR: wrong usage of cells_pre_init!. ", __FILE__, __LINE__ );
	// }
	// memcpy(&tmp_local,&local_cells,sizeof(CellPList));
	// dd_topology_init(&tmp_local);
}

void Particles::update_ghosts()
{
	int m, n, o, ind1, ind2;

	/*xy loop all local cells */
	for(n=1; n<(cell_grid[1]+1); n++)
		for(m=1; m<(cell_grid[0]+1); m++){
			o = 0;
			ind1 = get_linear_index(m,n,o,ghost_cell_grid);
			cells[ind1]->gPMove[0] = 0.0;
			cells[ind1]->gPMove[1] = 0.0;
			cells[ind1]->gPMove[2] =  (- local_box_l[2]);
			o = cell_grid[2];
			ind2 = get_linear_index(m,n,o,ghost_cell_grid);
			cells[ind1]->noGCell = ind2;

			o = cell_grid[2]+1;
			ind1 = get_linear_index(m,n,o,ghost_cell_grid);
			cells[ind1]->gPMove[0] = 0.0;
			cells[ind1]->gPMove[1] = 0.0;
			cells[ind1]->gPMove[2] =  (local_box_l[2]);
			o = 1;
			ind2 = get_linear_index(m,n,o,ghost_cell_grid);
			cells[ind1]->noGCell = ind2;
		}

		/* xz loop all local cells */
		for(o=1; o<(cell_grid[2]+1); o++) 
			for(m=1; m<(cell_grid[0]+1); m++){
				n = 0;
				ind1 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->gPMove[0] = 0.0;
				cells[ind1]->gPMove[1] = (- local_box_l[1]);
				cells[ind1]->gPMove[2] = 0.0;
				n = cell_grid[1];
				ind2 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->noGCell = ind2;

				n = cell_grid[1]+1;
				ind1 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->gPMove[0] = 0.0;
				cells[ind1]->gPMove[1] = (local_box_l[1]);
				cells[ind1]->gPMove[2] = 0.0;
				n = 1;
				ind2 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->noGCell = ind2;

			}

			/* yz loop all local cells */
			for(o=1; o<(cell_grid[2]+1); o++) 
				for(n=1; n<(cell_grid[1]+1); n++){
					m = 0;
					ind1 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->gPMove[0] = (- local_box_l[0]);
					cells[ind1]->gPMove[1] = 0.0;
					cells[ind1]->gPMove[2] = 0.0;
					m = cell_grid[0];
					ind2 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->noGCell = ind2;


					m = cell_grid[0]+1;
					ind1 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->gPMove[0] = (local_box_l[0]);
					cells[ind1]->gPMove[1] = 0.0;
					cells[ind1]->gPMove[2] = 0.0;
					m = 1;
					ind2 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->noGCell = ind2;
				}

				/* x loop all local cells */
				for(m=1; m<(cell_grid[0]+1); m++){
					o = 0;
					n = 0;
					ind1 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->gPMove[0] = 0.0;
					cells[ind1]->gPMove[1] = (- local_box_l[1]);
					cells[ind1]->gPMove[2] = (- local_box_l[2]);
					o = cell_grid[2];
					n = cell_grid[1];
					ind2 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->noGCell = ind2;

					o = cell_grid[2]+1;
					n = cell_grid[1]+1;
					ind1 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->gPMove[0] = 0.0;
					cells[ind1]->gPMove[1] = (local_box_l[1]);
					cells[ind1]->gPMove[2] = (local_box_l[2]);
					o =1;
					n =1;;
					ind2 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->noGCell = ind2;

					o = 0;
					n = cell_grid[1]+1;
					ind1 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->gPMove[0] = 0.0;
					cells[ind1]->gPMove[1] = (local_box_l[1]);
					cells[ind1]->gPMove[2] = (-local_box_l[2]);
					o = cell_grid[2];
					n = 1;
					ind2 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->noGCell = ind2;

					o = cell_grid[2]+1;
					n = 0;
					ind1 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->gPMove[0] = 0.0;
					cells[ind1]->gPMove[1] = (-local_box_l[1]);
					cells[ind1]->gPMove[2] = (local_box_l[2]);
					o = 1;
					n = cell_grid[1];
					ind2 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->noGCell = ind2;
				}

				/* y loop all local cells */
				for(n=1; n<(cell_grid[1]+1); n++){
					o = 0;
					m = 0;
					ind1 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->gPMove[0] = (-local_box_l[0]);
					cells[ind1]->gPMove[1] = 0.0;
					cells[ind1]->gPMove[2] = (-local_box_l[2]);
					o = cell_grid[2];
					m = cell_grid[0];
					ind2 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->noGCell = ind2;

					//cout << " 20312 "<<local_cells[203]->p_objects.size()<<endl;
					o = cell_grid[2]+1;
					m = cell_grid[0]+1;
					ind1 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->gPMove[0] = (local_box_l[0]);
					cells[ind1]->gPMove[1] = 0.0;
					cells[ind1]->gPMove[2] = (local_box_l[2]);
					o =1;
					m =1;;
					ind2 = get_linear_index(m,n,o,ghost_cell_grid);
					//cout << " 20315 "<<local_cells[203]->p_objects.size()<<endl;
					cells[ind1]->noGCell = ind2;

					//cout << " 20313 "<<local_cells[203]->p_objects.size()<<endl;
					o = 0;
					m = cell_grid[0]+1;
					ind1 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->gPMove[0] = (local_box_l[0]);
					cells[ind1]->gPMove[1] = 0.0;
					cells[ind1]->gPMove[2] = (-local_box_l[2]);
					o = cell_grid[2];
					m = 1;
					ind2 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->noGCell = ind2;

					o = cell_grid[2]+1;
					m = 0;
					ind1 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->gPMove[0] = (-local_box_l[0]);
					cells[ind1]->gPMove[1] = 0.0;
					cells[ind1]->gPMove[2] = (local_box_l[2]);
					o = 1;
					m = cell_grid[1];
					ind2 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->noGCell = ind2;
				}

				/* z loop all local cells */
				for(o=1; o<(cell_grid[2]+1); o++){
					m = 0;
					n = 0;
					ind1 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->gPMove[0] = (-local_box_l[0]);
					cells[ind1]->gPMove[1] = (-local_box_l[1]);
					cells[ind1]->gPMove[2] = 0.0;
					m = cell_grid[0];
					n = cell_grid[1];
					ind2 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->noGCell = ind2;

					m = cell_grid[0]+1;
					n = cell_grid[1]+1;
					ind1 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->gPMove[0] = (local_box_l[0]);
					cells[ind1]->gPMove[1] = (local_box_l[1]);
					cells[ind1]->gPMove[2] = 0.0;
					m =1;
					n =1;;
					ind2 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->noGCell = ind2;

					m = 0;
					n = cell_grid[1]+1;
					ind1 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->gPMove[0] = (-local_box_l[0]);
					cells[ind1]->gPMove[1] = (local_box_l[1]);
					cells[ind1]->gPMove[2] = 0.0;
					m = cell_grid[0];
					n = 1;
					ind2 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->noGCell = ind2;

					m = cell_grid[0]+1;
					n = 0;
					ind1 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->gPMove[0] = (local_box_l[0]);
					cells[ind1]->gPMove[1] = (-local_box_l[1]);
					cells[ind1]->gPMove[2] = 0.0;
					m = 1;
					n = cell_grid[1];
					ind2 = get_linear_index(m,n,o,ghost_cell_grid);
					cells[ind1]->noGCell = ind2;
				}

				m = 0;
				n = 0;
				o = 0;
				ind1 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->gPMove[0] = (-local_box_l[0]);
				cells[ind1]->gPMove[1] = (-local_box_l[1]);
				cells[ind1]->gPMove[2] = (-local_box_l[2]);
				m = cell_grid[0];
				n = cell_grid[1];
				o = cell_grid[2];
				ind2 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->noGCell = ind2;


				m = cell_grid[0]+1;
				n = 0;
				o = 0;
				ind1 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->gPMove[0] = (local_box_l[0]);
				cells[ind1]->gPMove[1] = (-local_box_l[1]);
				cells[ind1]->gPMove[2] = (-local_box_l[2]);
				m = 1;
				n = cell_grid[1];
				o = cell_grid[2];
				ind2 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->noGCell = ind2;


				m = 0;
				n = cell_grid[1]+1;
				o = 0;
				ind1 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->gPMove[0] = (-local_box_l[0]);
				cells[ind1]->gPMove[1] = (local_box_l[1]);
				cells[ind1]->gPMove[2] = (-local_box_l[2]);
				m = cell_grid[0];
				n = 1;
				o = cell_grid[2];
				ind2 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->noGCell = ind2;


				m = cell_grid[0]+1;
				n = cell_grid[1]+1;
				o = 0;
				ind1 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->gPMove[0] = (local_box_l[0]);
				cells[ind1]->gPMove[1] = (local_box_l[1]);
				cells[ind1]->gPMove[2] = (-local_box_l[2]);
				m = 1;
				n = 1;
				o = cell_grid[2];
				ind2 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->noGCell = ind2;


				m = 0;
				n = 0;
				o = cell_grid[2]+1;
				ind1 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->gPMove[0] = (-local_box_l[0]);
				cells[ind1]->gPMove[1] = (-local_box_l[1]);
				cells[ind1]->gPMove[2] = (local_box_l[2]);
				m = cell_grid[0];
				n = cell_grid[1];
				o = 1;
				ind2 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->noGCell = ind2;


				m = cell_grid[0]+1;
				n = 0;
				o = cell_grid[2]+1;
				ind1 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->gPMove[0] = (local_box_l[0]);
				cells[ind1]->gPMove[1] = (-local_box_l[1]);
				cells[ind1]->gPMove[2] = (local_box_l[2]);
				m = 1;
				n = cell_grid[1];
				o = 1;
				ind2 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->noGCell = ind2;


				m = 0;
				n = cell_grid[1]+1;
				o = cell_grid[2]+1;
				ind1 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->gPMove[0] = (-local_box_l[0]);
				cells[ind1]->gPMove[1] = (local_box_l[1]);
				cells[ind1]->gPMove[2] = (local_box_l[2]);
				m = cell_grid[0];
				n = 1;
				o = 1;
				ind2 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->noGCell = ind2;


				m = cell_grid[0]+1;
				n = cell_grid[1]+1;
				o = cell_grid[2]+1;
				ind1 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->gPMove[0] = (local_box_l[0]);
				cells[ind1]->gPMove[1] = (local_box_l[1]);
				cells[ind1]->gPMove[2] = (local_box_l[2]);
				m = 1;
				n = 1;
				o = 1;
				ind2 = get_linear_index(m,n,o,ghost_cell_grid);
				cells[ind1]->noGCell = ind2;

}

/** Init cell interactions for cell system domain decomposition.
* initializes the interacting neighbor cell list of a cell The
* created list of interacting neighbor cells is used by the verlet
* algorithm (see verlet.c) to build the verlet lists.
*/
void Particles::dd_init_cell_interactions()
{
	int m,n,o,p,q,r,ind1,ind2,c_cnt=0,n_cnt;
	Cell* temp_cell;

	/* initialize cell neighbor structures */
	//dd.cell_inter = (IA_Neighbor_List *) realloc(dd.cell_inter,local_cells.n*sizeof(IA_Neighbor_List));
	//for(m=0; m<local_cells.n; m++) { 
	//  dd.cell_inter[m].nList = NULL; 
	//  dd.cell_inter[m].n_neighbors=0; 
	//}

	/* loop all local cells */
	for(o=1; o<cell_grid[2]+1; o++) 
		for(n=1; n<cell_grid[1]+1; n++) 
			for(m=1; m<cell_grid[0]+1; m++){
				//dd.cell_inter[c_cnt].nList = (IA_Neighbor *) realloc(dd.cell_inter[c_cnt].nList, CELLS_MAX_NEIGHBORS*sizeof(IA_Neighbor));
				//dd.cell_inter[c_cnt].n_neighbors = CELLS_MAX_NEIGHBORS;

				n_cnt=0;
				ind1 = get_linear_index(m,n,o,ghost_cell_grid);

				temp_cell = local_cells[c_cnt];

				/* loop all neighbor cells */
				for(p=o-1; p<=o+1; p++)
					for(q=n-1; q<=n+1; q++)
						for(r=m-1; r<=m+1; r++) {
							ind2 = get_linear_index(r,q,p,ghost_cell_grid);
							if(ind2 >= ind1) {
								temp_cell->cell_inter[n_cnt].neighborCell = ind2;
								n_cnt++;
							}
						}
						c_cnt++;
			}
}

void Particles::addVelocity(double* vel)
{
	int i;
	double* tempv;

	int l_size = local_particles.size();

	for (i=0; i<l_size; i++ )
	{
		tempv = &local_particles[i]->v[0];
		*tempv += vel[0];
		tempv = &local_particles[i]->v[1];
		*tempv += vel[1];
		tempv = &local_particles[i]->v[2];
		*tempv += vel[2];
	}
}

void Particles::predict_momentum_particles(double* momentum)
{
	int i;
	sPart* pobjects;
	int l_size = local_particles.size();

	momentum[0] = 0.0;momentum[1] = 0.0;momentum[2] = 0.0;

	for (i=0; i<l_size; i++ )
	{
		pobjects = local_particles[i];
		momentum[0] +=  pobjects->v[0] + pobjects->f[0];
		momentum[1] +=  pobjects->v[1] + pobjects->f[1];
		momentum[2] +=  pobjects->v[2] + pobjects->f[2];
		//cout << " i "<<local_particles[i]->identity<< " v0 "<<local_particles[i]->v[0]<<" v1 "<<local_particles[i]->v[1]<<" v2 "<<local_particles[i]->v[2]<<" f0 "<<local_particles[i]->f[0]<<" f1 "<<local_particles[i]->f[1]<<" f2 "<<local_particles[i]->f[2]<<endl;
	}
	momentum[0] /= p_timestep;
	momentum[1] /= p_timestep;
	momentum[2] /= p_timestep;
}

void Particles::rescale_velocities(double scale) 
{
	int c;
	int l_size = local_particles.size();

	for (c = 0; c < l_size; c++) {
		local_particles[c]->v[0] *= scale;
		local_particles[c]->v[1] *= scale;
		local_particles[c]->v[2] *= scale;
	}
}

double Particles::calc_kinetic_energy()
{
	int i;
	sPart* pobjects;
	int l_size = local_particles.size();
	double k_energy;

	k_energy = 0.0;
	for (i=0; i<l_size; i++ )
	{
		pobjects = local_particles[i];
		k_energy += (SQR(pobjects->v[0]) + SQR(pobjects->v[1]) + SQR(pobjects->v[2]))*pobjects->mass;
		//cout << " i "<<local_particles[i]->identity<< " v0 "<<local_particles[i]->v[0]<<" v1 "<<local_particles[i]->v[1]<<" v2 "<<local_particles[i]->v[2]<<" f0 "<<local_particles[i]->f[0]<<" f1 "<<local_particles[i]->f[1]<<" f2 "<<local_particles[i]->f[2]<<endl;
	}
	k_energy /= (p_timestep*p_timestep*2.0);
	return k_energy;
}

/*  Particles: output  properties*/
void Particles::outProp(int t, FILE	*file_ptr)
{
	int n_sph;
	void   *io_ptr;

	io_ptr = (void *) &t;
	fwrite (io_ptr, sizeof(int), 1, file_ptr);
	for (n_sph = 0; n_sph < n_total_particles; n_sph++)
	{
		sPart* data = local_particles[n_sph];

		if ( data!= NULL) {
			unfold_position(data->p, data->i);
			io_ptr = (void *)     &data->p;
			fwrite (io_ptr, sizeof(data->p), 1, file_ptr);
			io_ptr = (void *)     &data->v;
			fwrite (io_ptr, sizeof(data->v), 1, file_ptr);
			io_ptr = (void *)     &data->f;
			fwrite (io_ptr, sizeof(data->f), 1, file_ptr);
			io_ptr = (void *)     &data->mass;
			fwrite (io_ptr, sizeof(data->mass), 1, file_ptr);
			io_ptr = (void *)     &data->type;
			fwrite (io_ptr, sizeof(data->type), 1, file_ptr);
			io_ptr = (void *)     &data->identity;
			fwrite (io_ptr, sizeof(data->identity), 1, file_ptr);
			//cout<<" data->p "<<data->p[0] << " data->i "<<data->identity<<" data->v "<<data->v[0]<<" data->f "<< data->f[0]<< " data->TYPE "<< data->type <<" data->mass " <<data->mass<<endl;
		}
	}
}

/*  Read properties file  */
void Particles::inProp(int t, FILE	*file_ptr)
{
	void   *io_ptr;
	int n_sph;
	int dt;

	sPart data;

	if (file_ptr != 0)
	{
		//cout << "open " << filename << " to read particle details." << endl;

		io_ptr = (void *) &dt;
		fread (io_ptr, sizeof(int), 1, file_ptr);
		if (t!=dt)throw ComException("time in the particle properties file does not match.", __FILE__, __LINE__ );

		for (n_sph = 0; n_sph < n_total_particles; n_sph++)
		{
			io_ptr = (void *)    &data.p;
			fread (io_ptr, sizeof(data.p), 1, file_ptr);
			io_ptr = (void *)    &data.v;
			fread (io_ptr, sizeof(data.v), 1, file_ptr);
			io_ptr = (void *)    &data.f;
			fread (io_ptr, sizeof(data.f), 1, file_ptr);
			io_ptr = (void *)    &data.mass;
			fread (io_ptr, sizeof(data.mass), 1, file_ptr);
			io_ptr = (void *)    &data.type;
			fread (io_ptr, sizeof(data.type), 1, file_ptr);
			io_ptr = (void *)     &data.identity;
			fread (io_ptr, sizeof(data.identity), 1, file_ptr);
			//cout << " ind "<<data.identity<<endl;

			//place_particle(data.identity, data.p, 1);
			sPart* temp_particle = local_particles[n_sph];
			temp_particle->identity = data.identity;
			temp_particle->mass = data.mass;
			temp_particle->type = data.type;
			memcpy(temp_particle->p, data.p, 3*sizeof(double));
			memcpy(temp_particle->v, data.v, 3*sizeof(double));
			memcpy(temp_particle->f, data.f, 3*sizeof(double));
		}
	}
}