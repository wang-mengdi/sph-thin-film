//////////////////////////////////////////////////////////////////////////
// MIT-style immersed boundary fluid driver
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

//NOTE: now this solver can only handle the situation that almost everywhere is fluid.
//The non-fluid material (like some solid object) must be very thin.

#ifndef __FluidEulerMetaIB_h__
#define __FluidEulerMetaIB_h__
#include "FluidEuler.h"

template<int d>
class FluidEulerMetaIB :public FluidEuler<d> {
	Typedef_VectorDii(d);
public:
	using Base = FluidEuler<d>;
	using Base::velocity;
	using Base::alpha;
	using Base::Update_Cell_Types;
	using Base::use_body_force;
	using Base::g;
	using Base::mac_grid;
	using Base::projection;
public:
	real epsilon;
	real kernel_coeff = 0.5;
	real max_vel = 0;
	//Base::velocity is the meta velocity
	Field<real, d> pressure;
	FaceField<real, d> flow_velocity;
	FaceField<real, d> body_velocity;

	virtual void Initialize(const VectorDi& cell_counts, const real dx, const VectorD& domain_min = VectorD::Zero());
	virtual void Advance(const real dt);//before advance, alpha and body_velocity must be filled
	void Retrieve_Flow_Velocity(void);//calculate flow_velocity from Base::velocity and body_velocity
	void Compose_Meta_Velocity(void);//calculate Base::velocity from flow_velocity and body_velocity
	void Meta_Advect_Velocity(real dt, FaceField<real, d>& advected_vel);
	virtual void Apply_Body_Forces(const real dt);//add on Base::velocity
	void Update_Max_Velocity(void);
	real Pressure_At(const VectorD& pos);
	
	real IB_Kernel_Fluid(real phi) const
	{
		//See: Accurate Cartesian-grid simulations of near-body flows at intermediate Reynolds numbers, Maertens&Weymouth
		//phi is the signed distance from nearest non-fluid material.
		//like, if there is a solid body, phi is negative inside the body, and positive outside
		if (phi < -epsilon)return 0;
		else if (phi > epsilon)return (real)1;
		else {
			real wb = (1 + cos(pi * phi / epsilon)) / 2;
			return 1 - wb * kernel_coeff;
		}
		//else return (real).5 * ((real)1 + phi / epsilon + sin(phi / epsilon * pi) / pi);
	}
	void Save_Snapshot(std::string snapshot_dir);
	void Load_Snapshot(std::string snapshot_dir);
};

#endif

