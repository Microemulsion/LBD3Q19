////////////////////////////////////////////////////////////////////////////////
// 04/07/2008 Song Kaixu
// Mesh.h
// Ver. 0.1 04/07/2008
// Class of mesh.
////////////////////////////////////////////////////////////////////////////////

#if !defined(_MESH_H_)
#define _MESH_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <math.h>
#include<assert.h>

#include "CommonC.h"

//Struct of mesh.
class BaseMesh{
public:
		BaseMesh( int ptype, int pdim, int pboundary, int pnx  )
            : d_type( ptype )
            , d_dim( pdim )
            , d_boundary( pboundary )
            {
				assert( pdim == 1 );
				switch (d_type) {
					case 0:
                        nx = pnx;
						ny = 0;
						nz = 0;
						break;
                    case 1:
                        nx = pnx;
						ny = 0;
						nz = 0;
						break;
                    case 2:
                        nr = pnx;
						break;
					case 4:
                        nr = pnx;
						nh = 0;
						break;
                    default:
                        throw ComException("wrong mesh type", __FILE__, __LINE__ );
                }
                setSize(); }
		BaseMesh( int ptype, int pdim, int pboundary, int pnx, int pny )
			: d_type( ptype )
            , d_dim( pdim )
            , d_boundary( pboundary )
            {
				assert( pdim == 2 );
				switch (d_type) {
					case 0:
                        nx = pnx;
						ny = pny;
						nz = 0;
						break;
                    case 1:
                        nx = pnx;
						ny = pny;
						nz = 0;
						break;
					case 2:
                        nr = pnx;
						break;
					case 3:
                        nr = pnx;
						nh = pny;
						break;
                    default:
                        throw ComException("wrong mesh type", __FILE__, __LINE__ );
                }
				setSize(); }
		BaseMesh( int ptype, int pdim, int pboundary, int pnx, int pny, int pnz )
			: d_type( ptype )
            , d_dim( pdim )
            , d_boundary( pboundary )
			, nx( pnx )
			, ny( pny )
			, nz( pnz )
            {
				assert( pdim == 3 );
				switch (d_type) {
					case 0:
                        nx = pnx;
						ny = pny;
						nz = pnz;
						break;
                    case 1:
                        nx = pnx;
						ny = pny;
						nz = pnz;
						break;
                    case 2:
                        nr = pnx;
						nh = pny;
						break;
					case 3:
                        nr = pnx;
						nh = pny;
						break;
                    default:
                        throw ComException("wrong mesh type", __FILE__, __LINE__ );
                }
				setSize(); }

		~BaseMesh(void){};

        enum BoundaryCondition {
            PERIODIC, // default
            RESTRICTION
            };
		enum MeshType{
            REGULARMESH,
            RECTANGULARMESH,
            SPHERICALMESH,
            CYLINDRICALMESH
            };
        int getNx() const { return nx; }
        int getNy() const { return ny; }
        int getNz() const { return nz; }
        int getNr() const { return nr; }
        int getNh() const { return nh; }
        int getdim() const { return d_dim; }
        int gettype() const { return d_type; }
        int getBoundaryC() const { return d_boundary; }
		int getSize() const { return d_size; }

protected:
		int    d_size;
        int    d_dim;
		int    nx, ny, nz, nr, nh;	//Axis name is x, y, z, r, h.
		int    d_type;	//Meshtype such as REGULARMESH 0, RECTANGULARMESH 1, SPHERICALMESH 2, CYLINDRICALMESH 3.
		int    d_boundary; //BoundaryCondition such as PERIODIC..
		
		int setSize() {
            if (d_type == 0) {
                switch (d_dim) {
                    case 1:
                        d_size = nx;
						break;
                    case 2:
                        d_size = (nx * ny);
						break;
                    case 3:
                        d_size = (nx * ny * nz);
						break;
                    default:
                        throw ComException("wrong mesh dim", __FILE__, __LINE__ );
                }
            } else {
                throw ComException("only work for scalar field now", __FILE__, __LINE__ );
                }
			return 0;
		}
};

//Add regular mesh.
//------------------------------------------------------------------------
class NormMesh:public BaseMesh{
public:
    NormMesh(	const BaseMesh& pmesh,
				double xMin,     double xMax,     int nx,
                double yMin=0.0, double yMax=0.0, int ny = 0,
                double zMin=0.0, double zMax=0.0, int nz = 0);

    enum BoundaryID {
        BoundaryXMin, BoundaryXMax,
        BoundaryYMin, BoundaryYMax,
        BoundaryZMin, BoundaryZMax
    };

    double   x(int id) const;    // x coordinate of the mesh point id
    double   y(int id) const;    // y coordinate of the mesh point id
    double   z(int id) const;    // z coordinate of the mesh point id
    double   xMin() const { return d_xMin; }
    double   xMax() const { return d_xMax; }
    double   yMin() const { return d_yMin; }
    double   yMax() const { return d_yMax; }
    double   zMin() const { return d_zMin; }
    double   zMax() const { return d_zMax; }
	double   lengthx() const { return d_xMax - d_xMin; }
	double   lengthy() const { return d_yMax - d_yMin; }
	double   lengthz() const { return d_zMax - d_zMin; }
    double   dx() const { return d_dx; }
    double   dy() const { return d_dy; }
    double   dz() const { return d_dz; }
    double   dx( int id ) const;
    double   dy( int id ) const;
    double   dz( int id ) const;
    int      nx() const { return d_Nx; }
    int      ny() const { return d_Ny; }
    int      nz() const { return d_Nz; }
    int      geometricalNx() const { return d_geometricalNx; }
    int      geometricalNy() const { return d_geometricalNy; }
    int      geometricalNz() const { return d_geometricalNz; }

protected:

    double d_xMin, d_xMax, d_dx;    
    double d_yMin, d_yMax, d_dy;    
    double d_zMin, d_zMax, d_dz;    
    int    d_Nx,   d_Ny,   d_Nz;

	int    d_xyNumMesh;          // number of Mesh in XY plane

    // the d_N considered geometrical boundaryconditon
    int    d_geometricalNx, d_geometricalNy, d_geometricalNz;
};

#endif // !defined(_MESH_H_)