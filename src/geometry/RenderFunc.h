//////////////////////////////////////////////////////////////////////////
// Assemble data for opengl_viewer. Especially for Points.
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "Common.h"
#include "GeometryParticles.h"

namespace RenderFunc {
	//each segment originates from xs[i], with direction normals[i] and length len_arr[i]*scale
	template<int d> void Write_Customized_Segments(std::string file_name, const Array<Vector<real, d> >& xs, const Array<Vector<real, d> >& normals, const Array<real>& len_arr, const real scale = 1.0);
	//write every scalar as a point. it's on the normal direction of particle i, and the distance between it and particle i is arr[i]*scale.
	template<int d> void Write_Scalars_As_Points(std::string file_name, const GeometryParticles<d>& particles, const Array<real>& arr, const real scale = 1.0);
}