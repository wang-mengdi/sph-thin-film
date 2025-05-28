//////////////////////////////////////////////////////////////////////////
// Assemble data for opengl_viewer. Especially for Points.
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "RenderFunc.h"
#include "AuxFunc.h"

namespace RenderFunc {
	template<int d>
	void Write_Customized_Segments(std::string file_name, const Array<Vector<real, d>>& xs, const Array<Vector<real, d>>& normals, const Array<real>& len_arr, const real scale)
	{
		Array<Vector<real, d> > to_write;
		int pn = (int)len_arr.size();
		to_write.resize(pn);
#pragma omp parallel for
		for (int i = 0; i < pn; i++) {
			to_write[i] = normals[i].normalized() * len_arr[i] * scale;
		}
		Write_Segments_To_File_3d_Fast<d, real>(xs, to_write, file_name);//in Particles.h
	}
	template void Write_Customized_Segments<2>(std::string file_name, const Array<Vector2>& xs, const Array<Vector2>& normals, const Array<real>& arr, const real scale);
	template void Write_Customized_Segments<3>(std::string file_name, const Array<Vector3>& xs, const Array<Vector3>& normals, const Array<real>& arr, const real scale);

	template<int d>
	void Write_Scalars_As_Points(std::string file_name, const GeometryParticles<d>& points, const Array<real>& arr, const real scale)
	{
		int pn = points.Size();
		if (arr.size() != pn) AuxFunc::Crash_With_Info("RenderFunc::Write_Scalars_As_Points error: size not match");
		Array<Vector<real, d>> to_write(pn);
#pragma omp parallel for
		for (int i = 0; i < pn; i++) {
			to_write[i] = points.X(i) + points.Normal(i).normalized() * arr[i] * scale;
		}
		Write_To_File_3d_Fast<d, real>(to_write, file_name);
	}
	template void Write_Scalars_As_Points<2>(std::string file_name, const GeometryParticles<2>& points, const Array<real>& arr, const real scale);
	template void Write_Scalars_As_Points<3>(std::string file_name, const GeometryParticles<3>& points, const Array<real>& arr, const real scale);
}
