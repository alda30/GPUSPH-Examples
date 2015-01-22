#ifndef _GPROBE_H
#define	_GPROBE_H

#include "Problem.h"
#include "Point.h"
#include "Rect.h"
#include "Cube.h"
#include "Cylinder.h"

#include "ode/ode.h"

class Gprobe: public Problem {
	private:
		Cube		experiment_box;
		PointVect	parts;
		PointVect	boundary_parts;
		PointVect	boundary_elems;
		PointVect	vertex_parts;
		VertexVect	vertex_indexes;

		float		h, w, l;
		float		H; // still water level
		bool		m_usePlanes; // use planes or boundaries

		// ODE stuff
		dGeomID		planes[5];
		dJointID	joint;
		Cylinder	cylinder;
		bool 		wet;	// set wet to true have a wet bed experiment

	public:
		Gprobe(const GlobalData *);
		~Gprobe(void);

		int fill_parts(void);
		uint fill_planes(void);
		void copy_to_array(BufferList &);
		void copy_planes(float4*, float*);

		void ODE_near_callback(void *, dGeomID, dGeomID);

		void release_memory(void);
};


#endif	/* _GPROBE_H */