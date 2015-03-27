#ifndef _GPROBEFALL_H
#define	_GPROBEFALL_H

#include <fstream>
#include "Problem.h"
#include "Point.h"
#include "Cube.h"
#include "Sphere.h"
#include "Cone.h"
#include "Cylinder.h"

#include "ode/ode.h"


class GprobeFall: public Problem {
	private:
		Cube					experiment_box;
		Cube 					fluid;
		PointVect				parts;
		PointVect				boundary_parts;
		PointVect				boundary_elems;
		PointVect				vertex_parts;
		VertexVect				vertex_indexes;
		double					H;				// still water level
		double					lx, ly, lz;		// dimension of experiment box
		dQuaternion 			rcube;
		bool					m_usePlanes; // use planes or boundaries

		// ODE stuff
		Cube					cube;
		Cylinder				cylinder;
		Cone 					cone;
		dGeomID					planes[5];
		dJointID				joint;
		float3 					ODEGravity;

		// ode output writing
		ofstream 				outputData;
		int 					intTime1, intTime2;
	
	public:
		GprobeFall(GlobalData *);
		virtual ~GprobeFall(void);

		int fill_parts(void);
		void copy_planes(float4*, float*);
		uint fill_planes(void);
		void copy_to_array(BufferList &);

		void ODE_near_callback(void *, dGeomID, dGeomID);

		float3 g_callback(const double);
		void release_memory(void);
};
#endif	/* _GPROBEFALL_H */

