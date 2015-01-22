#include <math.h>
#include <iostream>

#include "Gprobe.h"
#include "GlobalData.h"

#define CENTER_DOMAIN 1
// set to coords (x,y,z) if more accuracy is needed in such point
// (waiting for relative coordinates)
#if CENTER_DOMAIN
#define OFFSET_X (-l/2)
#define OFFSET_Y (-w/2)
#define OFFSET_Z (-h/2)
#else
#define OFFSET_X 0
#define OFFSET_Y 0
#define OFFSET_Z 0
#endif

Gprobe::Gprobe(const GlobalData *_gdata) : Problem(_gdata)
{
	H = 0.6f;
	wet = false;

	set_deltap(0.02f);

	l = 1.0f; w = l; h = 1;
	m_usePlanes = true;

	// SPH parameters
	m_simparams.dt = 0.0001f;			// initial timestep
	m_simparams.xsph = false;
	m_simparams.dtadapt = true;			// adaptive time step is selected
	m_simparams.dtadaptfactor = 0.3;	// safety factor in the adaptive time step formula
	m_simparams.buildneibsfreq = 10;	// frequency (in iterations) of neib list rebuilding
	m_simparams.shepardfreq = 0;		// frequency (in iterations) of Shepard density filter
	m_simparams.mlsfreq = 0;			// frequency (in iterations) of MLS density filter

	// Ferrari correction parameter should be (L/deltap)/1000, with L characteristic
	// length of the problem
	m_simparams.ferrari = H/(m_deltap*1000);	// coefficient  of Ferrari correction
	//m_simparams.visctype = KINEMATICVISC;
	m_simparams.visctype = DYNAMICVISC;
	//m_simparams.visctype = ARTVISC;
	m_simparams.mbcallback = false;				// we have no moving boundary
	m_simparams.boundarytype = SA_BOUNDARY;
	//m_simparams.boundarytype = LJ_BOUNDARY;

	// Size and origin of the simulation domain
	m_size = make_double3(l, w ,h);
	m_origin = make_double3(OFFSET_X, OFFSET_Y, OFFSET_Z);

	m_simparams.tend = 5.0;
	if (m_simparams.boundarytype == SA_BOUNDARY) {
		m_simparams.maxneibsnum = 256; // needed during gamma initialization phase // maximum number of neibs (should be a multiple of NEIBS_INTERLEAVE)

	};


	// Physical parameters
	m_physparams.gravity = make_float3(0.0, 0.0, -9.81f);
	const float g = length(m_physparams.gravity);
	const float maxvel = sqrt(g*H);
	// purely for cosmetic reason, let's round the soundspeed to the next
	// integer
	const float c0 = ceil(10*maxvel);
	m_physparams.set_density(0, 1800.0, 7.0f, c0); // the second parameter is rho (total (or initial) density); third one is gamma

	m_physparams.dcoeff = 5.0f*g*H;					// what is this?	

	m_physparams.r0 = m_deltap;						// influence radius of boundary repulsive force (for LJ-boundary)
	//m_physparams.visccoeff = 0.05f;
	m_physparams.kinematicvisc = 0.125f;			// if the value of dynamic viscosity is chosen at "m_physparams.artvisccoeff" why we need to choose kinematic viscosity?
	//m_physparams.kinematicvisc = 1.0e-6f;
	m_physparams.artvisccoeff = 200.0f;
	m_physparams.epsartvisc = 0.01*m_simparams.slength*m_simparams.slength;	// what is this? 
																			// we have not chosen smoothing length and by default its set to zero. Why? We dont need that?
	m_physparams.epsxsph = 0.5f;	// XSPH correction coefficient

	// Set ODE parameters:
	// Allocate data for floating bodies
	allocate_ODE_bodies(1);	// originally it was 2
	dInitODE();				// Initialize ODE
	m_ODEWorld = dWorldCreate();	// Create a dynamic world
	m_ODESpace = dHashSpaceCreate(0);
	m_ODEJointGroup = dJointGroupCreate(0);
	dWorldSetGravity(m_ODEWorld, m_physparams.gravity.x, m_physparams.gravity.y, m_physparams.gravity.z);	// Set gravity（x, y, z)


	// Drawing and saving times
	set_timer_tick(1.0e-4);
	add_writer(VTKWRITER, 1000);

	// Name of problem used for directory creation
	m_name = "Gprobe";
}


Gprobe::~Gprobe(void)
{
	release_memory();
	dWorldDestroy(m_ODEWorld);
	dCloseODE();
}


void Gprobe::release_memory(void)		// which objects should be 'cleared' in 'release_memory' function? for example why boundary_elems is not cleared?
{
	parts.clear();
	boundary_parts.clear();
}


int Gprobe::fill_parts()
{
	// distance between fluid box and wall
	float wd = m_physparams.r0;

	parts.reserve(14000); // what is this? or for example in 'DamBreak3D' what is 'boundary_parts.reserve(2000);'?

	experiment_box = Cube(Point(m_origin), Vector(l, 0, 0), Vector(0, w, 0), Vector(0, 0, h));

	experiment_box.SetPartMass(wd, m_physparams.rho0[0]);	// where is this function defined? what are the parameters?

	if(!m_usePlanes){
		if(m_simparams.boundarytype == SA_BOUNDARY) {
			experiment_box.FillBorder(boundary_parts, boundary_elems, vertex_parts, vertex_indexes, wd, false); // the last parameters is a boolean one to determine if the top face to be filled or not. (false = open)
		}
		else {
			experiment_box.FillBorder(boundary_parts, wd, false);
		}
	}


	//==========================================
	// Do we need to define IDE planes for box? if so, create it here
	planes[0] = dCreatePlane(m_ODESpace, 0.0, 0.0, 1.0, -m_origin.z);
	planes[1] = dCreatePlane(m_ODESpace, 1.0, 0.0, 0.0, -m_origin.x);
	planes[2] = dCreatePlane(m_ODESpace, -1.0, 0.0, 0.0, -m_origin.x + l);
	planes[3] = dCreatePlane(m_ODESpace, 0.0, 1.0, 0.0, -m_origin.y);
	planes[4] = dCreatePlane(m_ODESpace, 0.0, -1.0, 0.0, -m_origin.y + w);
	//==========================================


	Cube fluid = Cube(m_origin + Point(wd, wd, wd), Vector(l-2*wd, 0, 0), Vector(0, w-2*wd, 0), Vector(0, 0, H-2*wd));	// we've set the margin between the fluid cube and experiment_box to be "wd".
	
	fluid.SetPartMass(m_deltap, m_physparams.rho0[0]);
	// InnerFill puts particle in the center of boxes of step m_deltap, hence at
	// m_deltap/2 from the sides, so the total distance between particles and walls
	// is m_deltap = r0 = wd
	fluid.Fill(parts, m_deltap);

	//DEBUG: set only one fluid particle
//	parts.clear();
//	parts.push_back(Point(0.0, w/2.f, 0.0));
//	for(int i=0; i < vertex_parts.size(); i++)
//		if(	vertex_parts[i](2) == 0 &&
//			vertex_parts[i](0) > 0.5*w && vertex_parts[i](0) < 0.5*w+2*m_deltap &&
//			vertex_parts[i](1) > 0.5*w && vertex_parts[i](1) < 0.5*w+2*m_deltap)
//			parts.push_back(Point(vertex_parts[i](0) + 0.5*m_deltap, vertex_parts[i](1) + 0.5*m_deltap, 0.0));


	//===============Rigidbody: Graviprobe====================
	float fallingHeight = 3.0f;
	float GprobeRadius = 0.025f;
	float GrpobeHeight = 1.0f;
	float GprobeMass = 8.3f;
	float r0 = m_physparams.r0;
	cylinder = Cylinder(Point(l/2, w/2, fallingHeight), Vector(GrpobeHeight, 0.0, 0.0), Vector(0.0, 0.0, GrpobeHeight));
	cylinder.SetPartMass(r0, m_physparams.rho0[0]*0.3);
	cylinder.SetMass(r0, GprobeMass);
	cylinder.Unfill(parts, r0);					// what does this function do?
	cylinder.FillBorder(cylinder.GetParts(), r0);
	cylinder.ODEBodyCreate(m_ODEWorld, m_deltap);
	cylinder.ODEGeomCreate(m_ODESpace, m_deltap);
	add_ODE_body(&cylinder);
	//========================================================

	return parts.size() + boundary_parts.size() + vertex_parts.size() + get_ODE_bodies_numparts();;
}

uint Gprobe::fill_planes() // where is the source?
{
	return (m_usePlanes ? 5 : 0);
}

void Gprobe::copy_planes(float4 *planes, float *planediv)
{
	if (!m_usePlanes) return;
	
	planes[0] = make_float4(0, 0, 1.0, -m_origin.z);		// bottom plane
	planediv[0] = 1.0;
	planes[1] = make_float4(0, 1.0, 0, -m_origin.y);		// side plane (y near)
	planediv[1] = 1.0;
	planes[2] = make_float4(0, -1.0, 0, -m_origin.y + w);	// side plane (y far)
	planediv[2] = 1.0;
	planes[3] = make_float4(1.0, 0, 0, -m_origin.x);		// side plane (x near)
	planediv[3] = 1.0;
	planes[4] = make_float4(-1.0, 0, 0, -m_origin.x + l);	// side plane (x far)
	planediv[4] = 1.0;
}

//======================= For ODE Object ===========================
void Gprobe::ODE_near_callback(void *data, dGeomID o1, dGeomID o2)
{
	const int N = 10;
	dContact contact[N];

	int n = dCollide(o1, o2, N, &contact[0].geom, sizeof(dContact));
	// if ((o1 == cylinder.m_ODEGeom && o2 == sphere.m_ODEGeom) || (o2 == cylinder.m_ODEGeom && o1 == sphere.m_ODEGeom)) {
	// 	cout << "Collision between cube and obstacle " << n << "contact points\n";
	// }
	for (int i = 0; i < n; i++) {
		contact[i].surface.mode = dContactBounce;
		contact[i].surface.mu   = dInfinity;
		contact[i].surface.bounce     = 0.0; // (0.0~1.0) restitution parameter
		contact[i].surface.bounce_vel = 0.0; // minimum incoming velocity for bounce
		dJointID c = dJointCreateContact(m_ODEWorld, m_ODEJointGroup, &contact[i]);
		dJointAttach (c, dGeomGetBody(contact[i].geom.g1), dGeomGetBody(contact[i].geom.g2));
	}
}
//==================================================================

void Gprobe::copy_to_array(BufferList &buffers)
{
	float4 *pos = buffers.getData<BUFFER_POS>();
	hashKey *hash = buffers.getData<BUFFER_HASH>();
	float4 *vel = buffers.getData<BUFFER_VEL>();
	particleinfo *info = buffers.getData<BUFFER_INFO>();
	vertexinfo *vertices = buffers.getData<BUFFER_VERTICES>();
	float4 *boundelm = buffers.getData<BUFFER_BOUNDELEMENTS>();

	std::cout << "Boundary parts: " << boundary_parts.size() << "\n";
	for (uint i = 0; i < boundary_parts.size(); i++) {
		vel[i] = make_float4(0, 0, 0, m_physparams.rho0[0]);
		info[i] = make_particleinfo(BOUNDPART, 0, i);
		calc_localpos_and_hash(boundary_parts[i], info[i], pos[i], hash[i]);
	}
	int j = boundary_parts.size();
	std::cout << "Boundary part mass: " << pos[j-1].w << "\n";

	std::cout << "Fluid parts: " << parts.size() << "\n";
	for (uint i = j; i < j + parts.size(); i++) {
		float rho = density(H - parts[i-j](2), 0);
		vel[i] = make_float4(0, 0, 0, rho);
		info[i] = make_particleinfo(FLUIDPART, 0, i);
		calc_localpos_and_hash(parts[i-j], info[i], pos[i], hash[i]);
	}
	j += parts.size();
	std::cout << "Fluid part mass: " << pos[j-1].w << "\n";

	if (m_simparams.boundarytype == SA_BOUNDARY) {
			uint j = parts.size() + boundary_parts.size();

			std::cout << "Vertex parts: " << vertex_parts.size() << "\n";
		for (uint i = j; i < j + vertex_parts.size(); i++) {
			float rho = density(H - vertex_parts[i-j](2), 0);
			vel[i] = make_float4(0, 0, 0, rho);
			info[i] = make_particleinfo(VERTEXPART, 0, i);
			calc_localpos_and_hash(vertex_parts[i-j], info[i], pos[i], hash[i]);
		}
		j += vertex_parts.size();
		std::cout << "Vertex part mass: " << pos[j-1].w << "\n";

		if(vertex_indexes.size() != boundary_parts.size()) {
			std::cout << "ERROR! Incorrect connectivity array!\n";
			exit(1);
		}
		if(boundary_elems.size() != boundary_parts.size()) {
			std::cout << "ERROR! Incorrect boundary elements array!\n";
			exit(1);
		}

		uint offset = parts.size() + boundary_parts.size();
		for (uint i = 0; i < boundary_parts.size(); i++) {
			vertex_indexes[i].x += offset;
			vertex_indexes[i].y += offset;
			vertex_indexes[i].z += offset;

			vertices[i] = vertex_indexes[i];

			boundelm[i] = make_float4(boundary_elems[i]);
		}
	}
}
