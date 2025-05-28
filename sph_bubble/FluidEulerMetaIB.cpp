#include "FluidEulerMetaIB.h"


template<int d>
void FluidEulerMetaIB<d>::Initialize(const VectorDi& cell_counts, const real dx, const VectorD& domain_min)
{
	Base::Initialize(cell_counts, dx, domain_min);
	epsilon = dx*1;
	pressure.Resize(cell_counts, 0);
	flow_velocity.Resize(cell_counts, 0);
	body_velocity.Resize(cell_counts, 0);
	projection.use_alpha_for_correction = true;
	projection.use_multigrid_solver = false;
	projection.verbose = false;
}

template<int d>
void FluidEulerMetaIB<d>::Retrieve_Flow_Velocity(void)
{
	//meta_velocity=alpha*velocity+(1-alpha)*body_velocity
	for (int axis = 0; axis < d; axis++) {
		int face_num = mac_grid.face_grids[axis].node_counts.prod();
#pragma omp parallel for
		for (int idx = 0; idx < face_num; idx++) {
			const VectorDi& face = mac_grid.Face_Coord(axis, idx);
			real w = alpha(axis, face);
			if (w < 1e-15) flow_velocity(axis, face) = 0;
			else flow_velocity(axis, face) = (velocity(axis, face) - (1 - w) * body_velocity(axis, face)) / w;
		}
	}
}

template<int d>
void FluidEulerMetaIB<d>::Compose_Meta_Velocity(void)
{
	//meta_velocity=alpha*velocity+(1-alpha)*body_velocity
	for (int axis = 0; axis < d; axis++) {
		int face_num = mac_grid.face_grids[axis].node_counts.prod();
#pragma omp parallel for
		for (int idx = 0; idx < face_num; idx++) {
			const VectorDi& face = mac_grid.Face_Coord(axis, idx);
			real w = alpha(axis, face);
			velocity(axis, face) = w * flow_velocity(axis, face) + (1 - w) * body_velocity(axis, face);
		}
	}
}

template<int d>
void FluidEulerMetaIB<d>::Meta_Advect_Velocity(real dt, FaceField<real, d>& advected_vel)
{
	FaceField<real, d> old_vel = advected_vel;
	Advection::Semi_Lagrangian(dt, velocity, mac_grid, old_vel, mac_grid, advected_vel);
	Update_Cell_Types();
}

template<int d>
void FluidEulerMetaIB<d>::Apply_Body_Forces(const real dt)
{
	if (!use_body_force)return;
	for (int axis = 0; axis < d; axis++) {
		int face_num = mac_grid.Number_Of_Faces(axis);
#pragma omp parallel for
		for (int i = 0; i < face_num; i++) {
			VectorDi face = mac_grid.Face_Coord(axis, i);
			velocity.face_fields[axis](face) += g[axis] * dt * alpha(axis, face);
		}
	}
}

template<int d>
void FluidEulerMetaIB<d>::Update_Max_Velocity(void)
{
	max_vel = 0.0;
	iterate_face(axis, iter, mac_grid) {
		max_vel = std::max(max_vel, (real)fabs((double)velocity(axis, iter.Coord())));
	}
}

template<int d>
real FluidEulerMetaIB<d>::Pressure_At(const VectorD& pos)
{
	Interpolation<d> intp(mac_grid);
	return intp.Interpolate_Centers(pressure, pos);
}

template<int d>
void FluidEulerMetaIB<d>::Save_Snapshot(std::string snapshot_dir)
{
	velocity.Write_Binary(snapshot_dir + "/velocity.ib");
	flow_velocity.Write_Binary(snapshot_dir + "/flow_velocity.ib");
}

template<int d>
void FluidEulerMetaIB<d>::Load_Snapshot(std::string snapshot_dir)
{
	velocity.Read_Binary(snapshot_dir + "/velocity.ib");
	flow_velocity.Read_Binary(snapshot_dir + "/flow_velocity.ib");
}

template<int d>
void FluidEulerMetaIB<d>::Advance(const real dt)
{
	Meta_Advect_Velocity(dt, flow_velocity);
	Compose_Meta_Velocity();
	Apply_Body_Forces(dt);
	Base::Enforce_Incompressibility();
	projection.Pressure(pressure);
	Retrieve_Flow_Velocity();
	Update_Max_Velocity();
}

template class FluidEulerMetaIB<2>;
template class FluidEulerMetaIB<3>;