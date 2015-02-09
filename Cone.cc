/*  Copyright 2011-2013 Alexis Herault, Giuseppe Bilotta, Robert A. Dalrymple, Eugenio Rustico, Ciro Del Negro

    Istituto Nazionale di Geofisica e Vulcanologia
        Sezione di Catania, Catania, Italy

    Università di Catania, Catania, Italy

    Johns Hopkins University, Baltimore, MD

    This file is part of GPUSPH.

    GPUSPH is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GPUSPH is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GPUSPH.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "Cone.h"


Cone::Cone(void)
{
	m_origin = Point(0, 0, 0);
	m_rt = 0.0;
	m_rb = 0.0;
	m_h = 0.0;
	m_hg = 0.0;
	m_halfaperture = 0;
}


Cone::Cone(const Point& center, const double radiusbottom, const double radiustop, const Vector& height)
{
	m_origin = center;
	m_rb = radiusbottom;
	m_rt = radiustop;
	m_h = height.norm();

	m_halfaperture = atan((m_rb - m_rt)/m_h);

	Vector v(0, 0, 1);
	const double angle = acos(height*v/m_h);
	Vector rotdir = -height.cross(v);
	if (rotdir.norm() == 0)
		rotdir = Vector(0, 1, 0);
	m_ep = EulerParameters(rotdir, angle);
	m_ep.ComputeRot();

	m_hg = m_h*(m_rb*m_rb + 2.0*m_rb*m_rt + 3.0*m_rt*m_rt)/
				(4.0*(m_rb *m_rb + m_rb*m_rt + m_rt*m_rt));

	m_center = m_origin + m_ep.Rot(m_hg*v);
}


Cone::Cone(const Point& center, const double radiusbottom, const double radiustop, const double height, const EulerParameters&  ep)
{
	m_origin = center;
	m_rb = radiusbottom;
	m_rt = radiustop;
	m_h = height;

	m_halfaperture = atan((m_rb - m_rt)/m_h);

	m_ep = ep;
	m_ep.ComputeRot();

	m_hg = m_h*(m_rb*m_rb + 2.0*m_rb*m_rt + 3.0*m_rt*m_rt)/
				(4.0*(m_rb *m_rb + m_rb*m_rt + m_rt*m_rt));

	m_center = m_origin + m_hg*m_ep.Rot(Vector(0, 0, 1));
}


Cone::Cone(const Point& center, const Vector& radiusbottom, const Vector& radiustop, const Vector& height)
{
	if (abs(radiusbottom*height) > 1e-8*radiusbottom.norm()*height.norm()
		|| abs(radiustop*height) > 1e-8*radiustop.norm()*height.norm()) {
		std::cout << "Trying to construct a cone with non perpendicular radius and axis\n";
		std::exit(1);
	}

	m_origin = center;
	m_rb = radiusbottom.norm();
	m_rt = radiustop.norm();
	m_h = height.norm();

	Vector radiusdir = height.Normal();
	Vector generatrix(m_origin + m_rb*radiusdir, m_origin + height + m_rt*radiusdir);
	m_halfaperture = acos(height*generatrix/(height.norm()*generatrix.norm()));

	Vector v(0, 0, 1);
	const double angle = acos(height*v/m_h);
	Vector rotdir = height.cross(v);
	if (rotdir.norm() == 0)
		rotdir = Vector(0, 1, 0);
	m_ep = EulerParameters(rotdir, angle);
	m_ep.ComputeRot();

	m_hg = m_h*(m_rb*m_rb + 2.0*m_rb*m_rt + 3.0*m_rt*m_rt)/
				(4.0*(m_rb *m_rb + m_rb*m_rt + m_rt*m_rt));

	m_center = m_origin + m_hg*m_ep.Rot(Vector(0, 0, 1));

}


double
Cone::Volume(const double dx) const
{
	const double h = m_h + dx;
	const double rb = m_rb + dx/2.0;
	const double rt = m_rt + dx/2.0;

	const double volume = M_PI*h/3.0*(rb*rb + rb*rt + rt*rt);
	return volume;
}


void
Cone::SetInertia(const double dx)
{
	const double h = m_h + dx;
	const double rb = m_rb + dx/2.0;
	const double rt = m_rt + dx/2.0;

	const double d = 20.0*M_PI*(rb*rt + rb*rb + rt*rt);
	const double n1 = 2.0*h*h*(rb*rb + 3.0*rb*rt + 6.0*rt*rt);
	const double n2 = 3.0*(rb*rb*rb*rt + rb*rb*rt*rt + rb*rt*rt*rt + rb*rb*rb*rb + rt*rt*rt*rt);

	m_inertia[0] = m_mass*(n1 + n2)/d;
	m_inertia[1] = m_inertia[0];
	m_inertia[2] = 2.0*n2*m_mass/d;

	std::cout << "Inertia: " << m_inertia[0] << " " << m_inertia[1] << " " << m_inertia[2] << "\n";

}


void
Cone::FillBorder(PointVect& points, const double dx, const bool bottom, const bool top)
{
	m_origin(3) = m_center(3);
	const int nz = (int) ceil(m_h/dx);
	const double dz = m_h/nz;
	for (int i = 0; i <= nz; i++)
		FillDiskBorder(points, m_ep, m_origin, m_rb - i*dz*sin(m_halfaperture), i*dz, dx, 2.0*M_PI*rand()/RAND_MAX);
	if (bottom)
		FillDisk(points, m_ep, m_origin, m_rb - dx, 0.0, dx);
	if (top)
		FillDisk(points, m_ep, m_origin, m_rt - dx, nz*dz, dx);
}


int
Cone::Fill(PointVect& points, const double dx, const bool fill)
{
	m_origin(3) = m_center(3);
	int nparts = 0;
	const int nz = (int) ceil(m_h/dx);
	const double dz = m_h/nz;
	for (int i = 0; i <= nz; i++)
		nparts += FillDisk(points, m_ep, m_origin, m_rb - i*dz*sin(m_halfaperture), i*dz, dx, fill);

	return nparts;
}


bool
Cone::IsInside(const Point& p, const double dx) const
{
	Point lp = m_ep.TransposeRot(p - m_origin);
	const double h = m_h + dx;
	bool inside = false;
	const double z = lp(2);
	if (z > -dx && z < h) {
		const double r = m_rb - z*sin(m_halfaperture) + dx;
		if (lp(0)*lp(0) + lp(1)*lp(1) < r*r)
			inside = true;
	}

	return inside;
}

// Added to make it possible for Cone shape to be used in Rigid body part of GPUSPH

//######## GPUSPH FUNCTIONS FOR CONE SHAPE ########

void
Cone::ODEBodyCreate(dWorldID ODEWorld, const double dx, dSpaceID ODESpace)
{
	m_ODEBody = dBodyCreate(ODEWorld);
	dMassSetZero(&m_ODEMass);
	//dMassSetCylinderTotal(&m_ODEMass, m_mass, 3, m_r +dx/2.0, m_h + dx);
	//dMassSetSphereTotal(&m_ODEMass, m_mass, m_r + dx/2.0);
	//Presumably we have to use dMassSetTrimesh()
	dMassSetTrimeshTotal(&m_ODEMass, m_mass, m_ODEGeom);

	dBodySetMass(m_ODEBody, &m_ODEMass);
	dBodySetPosition(m_ODEBody, m_center(0), m_center(1), m_center(2));
	if (ODESpace)
		ODEGeomCreate(ODESpace, dx);
}

void
Cone::ODEGeomCreate(dSpaceID ODESpace, const double dx)
{
	// defining vertices
	// #0
	Vertices[0][0] = 0;
	Vertices[0][1] = m_rt;
	Vertices[0][2] = 0.25 * m_h;
	// #1
	Vertices[1][0] = 0.6 * m_rt;
	Vertices[1][1] = 0.8 * m_rt;
	Vertices[1][2] = 0.25 * m_h;
	// #2
	Vertices[2][0] = 0.96 * m_rt;
	Vertices[2][1] = 0.308 * m_rt;
	Vertices[2][2] = 0.25 * m_h;
	// #3
	Vertices[3][0] = 0.96 * m_rt;
	Vertices[3][1] = -0.308 * m_rt;
	Vertices[3][2] = 0.25 * m_h;
	// #4
	Vertices[4][0] = 0.6 * m_rt;
	Vertices[4][1] = -0.02;
	Vertices[4][2] = 0.25 * m_h;
	// #5
	Vertices[5][0] = 0;
	Vertices[5][1] = -m_rt;
	Vertices[5][2] = 0.25 * m_h;
	// #6
	Vertices[6][0] = -0.6 * m_rt;
	Vertices[6][1] = -0.8 * m_rt;
	Vertices[6][2] = 0.25 * m_h;
	// #7
	Vertices[7][0] = -0.96 * m_rt;
	Vertices[7][1] = -0.308 * m_rt;
	Vertices[7][2] = 0.25 * m_h;
	// #8
	Vertices[8][0] = -0.96 * m_rt;
	Vertices[8][1] = 0.308 * m_rt;
	Vertices[8][2] = 0.25 * m_h;
	// #9
	Vertices[9][0] = -0.6 * m_rt;
	Vertices[9][1] = 0.8 * m_rt;
	Vertices[9][2] = 0.25 * m_h;
	// #10
	Vertices[10][0] = 0;
	Vertices[10][1] = 0;
	Vertices[10][2] = -0.75 * m_h;

	// defining facets
	// #0
	Indices[0] = 0;
  	Indices[1] = 10;
  	Indices[2] = 1;
  	// #1
	Indices[3] = 1;
  	Indices[4] = 10;
  	Indices[5] = 2;
	// #2
	Indices[6] = 2;
  	Indices[7] = 10;
  	Indices[8] = 3;
  	// #3
	Indices[9] = 3;
  	Indices[10] = 10;
  	Indices[11] = 4;
  	// #4
	Indices[12] = 4;
  	Indices[13] = 10;
  	Indices[14] = 5;
  	// #5
	Indices[15] = 5;
  	Indices[16] = 10;
  	Indices[17] = 6;
  	// #6
	Indices[18] = 6;
  	Indices[19] = 10;
  	Indices[20] = 7;
  	// #7
	Indices[21] = 7;
  	Indices[22] = 10;
  	Indices[23] = 8;
  	// #8
	Indices[24] = 8;
  	Indices[25] = 10;
  	Indices[26] = 9;
  	// #9
	Indices[27] = 9;
  	Indices[28] = 10;
  	Indices[29] = 0;

	Data = dGeomTriMeshDataCreate();
	dGeomTriMeshDataBuildSingle(Data, Vertices[0], 3 * sizeof(float), VertexCount, &Indices[0], IndexCount, 3 * sizeof(dTriIndex));

	m_ODEGeom = dCreateTriMesh(ODESpace, Data, 0, 0, 0);

	if (m_ODEBody)
		dGeomSetBody(m_ODEGeom, m_ODEBody);
	else {
		dGeomSetPosition(m_ODEGeom, m_center(0), m_center(1), m_center(2));
	}
}