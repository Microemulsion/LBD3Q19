// Mesh.cpp: Class of mesh.
//
//////////////////////////////////////////////////////////////////////

#include "Mesh.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

NormMesh::
    NormMesh(	const BaseMesh& pmesh,
				double xMin, double xMax, int nx,
                double yMin, double yMax, int ny,
                double zMin, double zMax, int nz) 
				:   BaseMesh( pmesh.gettype(), pmesh.getdim(), pmesh.getBoundaryC(), pmesh.getNx(), pmesh.getNy(), pmesh.getNz() ),
                    d_xMin( xMin ), d_xMax( xMax ), d_Nx( nx ), d_geometricalNx( nx ),
                    d_yMin( yMin ), d_yMax( yMax ), d_Ny( ny ), d_geometricalNy( ny ),
                    d_zMin( zMin ), d_zMax( zMax ), d_Nz( nz ), d_geometricalNz( nz ),
                    d_dx( 0. ), d_dy( 0. ), d_dz( 0. )
{
    if( d_xMin > d_xMax || d_yMin > d_yMax || d_zMin > d_zMax ) {
        throw ComException("no good data", __FILE__, __LINE__ );
    }

    if( d_Nx <= 0 ) {
        throw ComException("nx <= 0", __FILE__, __LINE__ );
    }

    d_dx = ( d_xMax - d_xMin ) / d_Nx;
    if(  d_boundary == 1 ) {
        d_geometricalNx--;
    }

    if( d_Ny > 0 ) {
        d_dy = ( d_yMax - d_yMin ) / d_Ny;
        if( d_boundary == 1 ) {
            d_geometricalNy--;
        }
    }

    if( d_Nz > 0 ) {
        d_dz = ( d_zMax - d_zMin ) / d_Nz;
        if(  d_boundary == 1 ) {
            d_geometricalNz--;
        }
    }

    d_xyNumMesh = ( d_geometricalNx + 1 ) * ( d_geometricalNy + 1 );

    d_size = d_xyNumMesh * ( d_geometricalNz + 1 );
}