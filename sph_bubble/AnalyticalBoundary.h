//////////////////////////////////////////////////////////////////////////
// Analytical boundary for FluidSPHBubble
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __AnalyticalBoundary_h__
#define __AnalyticalBoundary_h__
#include "Common.h"
#include "GeometryPrimitives.h"

template<int d>
class AnalyticalBoundary {
	Typedef_VectorDii(d);
public:
	Array<std::shared_ptr<ImplicitGeometry<d> > > obstacles;//fluid should be outside of ALL of these
	Array<std::shared_ptr<ImplicitGeometry<d> > > boundaries;//fluid should be inside of ALL of these 
	void Add_Obstacle(std::shared_ptr<ImplicitGeometry<d> > obj) { obstacles.push_back(obj); }//outside of this
	void Add_Bounding(std::shared_ptr<ImplicitGeometry<d> > obj) { boundaries.push_back(obj); }//inside of this
	bool Available() { return (obstacles.size() + boundaries.size() > 0.); }
	bool Get_Nearest_Boundary(const VectorD& pos, real& dis, VectorD& normal) {
		bool detected = false;//for further extension. may add a list of "inside" geometries
		for (int i = 0; i < obstacles.size(); i++) {
			real phi = obstacles[i]->Phi(pos);
			if (!detected || phi < dis) {
				detected = true;
				dis = phi;
				normal = obstacles[i]->Normal(pos);
			}
		}
		for (int i = 0; i < boundaries.size(); i++) {
			real phi = boundaries[i]->Phi(pos);
			if (!detected || (-phi) < dis) {
				detected = true;
				dis = -phi;
				normal = -obstacles[i]->Normal(pos);
			}
		}
		if (!detected) { std::cerr << "AnalyticalBoundary<d>::Get_Nearest_Boundary error: no analytical boundaries\n"; exit(0); }
		return detected;
	}

	bool Get_All_Boundaries(const VectorD& pos, Array<real>& dis, VectorD& normal) {
		dis.clear();
		bool detected = false;//for further extension. may add a list of "inside" geometries
		for (int i = 0; i < obstacles.size(); i++) {
			real phi = obstacles[i]->Phi(pos);
			detected = true;
			dis.push_back(phi);
		}
		for (int i = 0; i < boundaries.size(); i++) {
			real phi = boundaries[i]->Phi(pos);
			detected = true;
			dis.push_back(phi);
		}
		return detected;
	}
};


#endif