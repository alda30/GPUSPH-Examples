#ifndef _TRIAL_H
#define	_TRIAL_H

#include "Problem.h"
#include "Point.h"
#include "Cube.h"
#include "Sphere.h"
#include "Cone.h"
#include "Torus.h"
#include "Cylinder.h"

#include "ode/ode.h"

class Trial: public Problem {
	private:
		Cube		experiment_box;
		Cube		obstacle;
		PointVect	parts;
		PointVect	boundary_parts;
		PointVect	boundary_elems;
		PointVect	vertex_parts;
		VertexVect	vertex_indexes;
		PointVect	obstacle_parts;
		double		H;				// still water level
		double		lx, ly, lz;		// dimension of experiment box
		bool		wet;			// set wet to true have a wet bed experiment
		bool		m_usePlanes; // use planes or boundaries

		// ODE stuff
		Sphere		sphere;
		Cube		cube;
		Cylinder	cylinder;
		dGeomID		planes[5];
		dJointID	joint;
		float3 		ODEGravity;


	public:
		Trial(const GlobalData *);
		virtual ~Trial(void);

		int fill_parts(void);
		void copy_planes(float4*, float*);
		uint fill_planes(void);
		void copy_to_array(BufferList &);

		void ODE_near_callback(void *, dGeomID, dGeomID);

		float3 g_callback(const float);
		void release_memory(void);
};
#endif	/* _TRIAL_H */

