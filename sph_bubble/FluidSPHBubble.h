//////////////////////////////////////////////////////////////////////////
// SPH bubble fluid driver
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __FluidSPHBubble_h__
#define __FluidSPHBubble_h__
#include "PointSet.h"
#include "SPHBubbleParticles.h"
#include "Kernels.h"
#include "RandomNumber.h"
//#include "FluidEulerMetaIB.h"
#include "EulerInitializer.h"
#include "AnalyticalBoundary.h"
#include "ArrayIO.h"
#include "BubbleParams.h"
#include "Fluid3DSPH.h"
#include "Timer.h"

#include <numeric>      // std::iota
#include<iomanip>
#include <fstream>

template<int d>
Vector<real, d> Project_Vector_Along(const Vector<real, d>& v, const Vector<real, d>& dir) {
	return v.dot(dir) / (dir.dot(dir)) * dir;//project vector a to b == (a dot b)/(b dot b) * b
}

template<class T>
std::function<T(const int)> Index_Function(const Array<T>& arr) {
	std::function<T(const int)> func = [&](const int idx)->T {return arr[idx]; };
	return func;
}

template<class T, int d>
void Print_Vector_Array(const Array<Vector<T, d> >& arr) {
	for (int i = 0; i < arr.size(); i++) {
		std::cout << "( ";
		for (int axis = 0; axis < d; axis++) {
			std::cout << arr[i][axis] << " ";
		}
		std::cout << "),";
		//std::cout << "i: " << arr[i].transpose() << "\n";
	}
	std::cout << "\n";
}

template<class T>
void Print_Scalar_Array(const Array<T>& arr) {
	for (int i = 0; i < arr.size(); i++) {
		std::cout << arr[i] << ",";
	}
	std::cout << "\n";
}



class DefaultSimParams {
public:
	real gravity_coeff;
	real viscosity_coeff;
	real normal_viscosity_coeff;
	real marangoni_coeff = 0.;
	OperatorParams grad_force_params;
	OperatorParams height_laplacian_params;
	OperatorParams geometry_params;
	DefaultSimParams() {}
	DefaultSimParams(real _dx, real _rho, real _viscosity_water, real _gamma_water, real _thickness, real _simulation_scale) {
		real alpha_vis = _dx * _dx * _rho / (_viscosity_water) * 50 * 2;
		viscosity_coeff = 0.1 * 2 * alpha_vis;
		normal_viscosity_coeff = 0.;
		gravity_coeff = 0.1 * 0.1 * 0.1;
		grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
		height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);
		geometry_params = OperatorParams("only_fluid", KernelType::SPIKY);
	}
};

template<int d> class FluidSPHBubble;

template<int d>
class MirrorCompensator {
	Typedef_VectorDii(d);
	using VectorT = Vector<real, d - 1>;
public:
	FluidSPHBubble<d>* fluid;
	int bnd_nb = -1;
	VectorT bnd_vec; real bnd_sqrval;
	bool self_bnd = false;
	real weight;
	MirrorCompensator() {}
	void Initialize(FluidSPHBubble<d>* _fluid, int idx, real _weight);
	real Coeff_Offset(VectorT lr_ij, int order);
};

template<int d> class FluidSPHBubble
{
	Typedef_VectorDii(d);
	Typedef_MatrixD(d);
public:
	using VectorT = Vector<real, d - 1>;								////tangential vector type
	SPHBubbleParticles<d> particles;
	SPHBubbleParticles<d> added_particles;
	PointSet<d> surface;
	std::shared_ptr<KernelSPH<d - 1>> kernel;						////kernel in local frame
	std::shared_ptr<KernelSPH<d>> kernel_v;

	OperatorParams geometry_params, viscosity_params, height_laplacian_params, grad_force_params, divergence_params;

	//// basic geometry
	int np_on_h = 6;							////number of particles on the kernel radius, given the input length scale, this parameter controls the radius size
	real length_scale=0.0;						////typical distance between two neighboring particles

	real max_vel = 0;
	real simulation_scale = 1;
	real total_volume = 0;

	//// physics
	real thickness = 1.0;//characteristic thickness
	real gamma_water = 1; //water surface tension
	real gamma_soap = 1; //soap surface tension --actually not soap surface tension, but how intense soap cancels out water surface tension
	real viscosity_water = 1.005e-3;//viscosity of water at 20 centigrade
	VectorD g = VectorD::Unit(1) * (real)-1.;	////gravity

	//// world forces
	real gravity_coeff = 1;
	real friction_coeff = 0;
	std::function<VectorD(const VectorD&)> air_velocity_func = nullptr;
	std::function<VectorD(const int)> external_force_func = nullptr;

	//// tangential force
	VorticityParams vorticity_params;//nothing by default structor
	TangentialPressureParams t_pressure_params;
	DefaultSimParams default_sim_params;//nothing by default structor
	real viscosity_coeff = 150;
	real marangoni_coeff = 4; ////marangoni force scale
	bool diffuse_soap = true;
	//// normal force
	NormalPressureParams n_pressure_params;
	real normal_viscosity_coeff = 150;

	//disabled for open source
	//FluidEulerMetaIB<d> air_solver;

	//// rh
	RenderHeightParams rh_params;
	//// boundary
	BoundaryParams boundary_params;
	AnalyticalBoundary<d> analytical_boundary;
	bool delete_solitary = true;
	bool delete_idle = false;

	//// control
	bool verbose = false;						////print out dbg information
	bool diagnosis = true, timing = false;

	//// temp or for debug use
	int iter = 0, current_frame = 0, last_frame = 0, first_frame = 0;
	real default_mirror_weight = 1.;

	real avg_height = 0.;
	///diagnosis arrays
	Array<real> total_vis_force;
	Array<real> total_ext_force;
	Array<real> total_pressure_force;
	Array<real> total_height_pressure_force;
	Array<real> total_height_laplacian_pressure_force;
	Array<real> total_div_pressure_force;
	Array<real> total_bnd_pressure_force;
	Array<real> total_marangoni_force;

	Array<real> phis;
	Array<VectorD> bnd_normals;

	bool clip_velocity = false;
	real vel_threshold = 1.; //allows max velocity of this for each force applied

	bool seeding_droplet = false;
	int last_droplet_frame = -1;
	real seed_droplet_rate = 0.3;
	real droplet_rh = thickness;
	real droplet_m = 1;
	int max_droplets = 28;
	int num_droplets = 0;
	Array<int> droplet_nozzles; //an array of boundary particles to yield droplet for a number of frames
	Array<int> drops;

	std::string exp_mode = "generic";
	real catenoid_speed = 0.03;
	bool delete_speedy_particles = false;
	bool quit_when_too_fast = true;

	real mark_solitary_interval = 1. / 200.;
	real last_mark_solitary_time = 0.;
	Array<int> to_explode;
	real explode_time = 0;
	bool has_exploded = 0;
	bool has_leaked = 0;
	bool merge_3d = false;
	int solitary_num = 0;
	real init_explosion_time = 0.;
	bool has_init_explosion = false;

	bool use_multiphase = false;

	VectorD g_init = 9.8 * VectorD::Unit(0);
	VectorD g_final = 9.8 * VectorD::Unit(0);
	bool interp_g = false;

	bool dynamic_seed_vortex = false;
	real dynamic_seed_rate = 0.05;
	int last_seed_frame = -1;

	Fluid3DSPH<d>* fluid_3d = nullptr;

	int catenoid_stop_frame = 100;
	real catenoid_stop_rate = 0.95;
	std::ofstream output_file;

	real replenish_V_rate = 0.5;

	void Initialize(real _length_scale, EulerInitializer<d>* perimeter = nullptr, Fluid3DSPH<d>* solver3d = nullptr)
	{
		if (perimeter != nullptr) {
			assert(false && "FluidSPHBubble::Initialize: perimeter is not supported in open source version");

			// std::cout << "initialize with grid solver\n";
			// n_pressure_params.air_pressure_mode = "ib";
			// air_solver.Initialize(perimeter->cell_counts, perimeter->dx, perimeter->domain_min);
			// perimeter->Fill_Boundary_Condition(air_solver.bc);
			// air_solver.use_body_force = true;
			// air_solver.g = g * gravity_coeff;
		}

		fluid_3d = solver3d;

		//lagrangian initialization
		length_scale = _length_scale;

		surface.Initialize(length_scale, np_on_h, &particles);
		surface.dif_type = PointSet<d>::DifType::SPH;
		//surface.v_r = 5 * _length_scale; //set surface neighbor searcher radius
		kernel = std::make_shared<KernelSPH<d - 1> >(surface.t_r);	////initialize kernel for the tangential plane
		kernel_v = std::make_shared<KernelSPH<d>>(surface.v_r);	////initialize kernel for the volumetric sphere	

		geometry_params = OperatorParams("only_fluid", KernelType::QUINTIC);
		
		viscosity_params= OperatorParams("only_fluid", KernelType::SPIKY);
		height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);
		grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
		divergence_params = OperatorParams("all", KernelType::SPIKY);

		Print_Params();
		Recompute_Geometries();

		avg_height = 0; real num = 0;
		for (int i = 0; i < particles.Size();i++) {
			if (particles.Is_Boundary(i)) continue;
			avg_height += particles.H(i); 
			num += 1.0;
		}
		avg_height /= num;

		std::cout << "avg_height at init: " << avg_height << std::endl;

		Update_Surface_Tension();
		Update_Fluid_Pressure();
		total_volume = 0; for (int i = 0; i < particles.Size(); i++) total_volume += particles.SA(i) * particles.RH(i);
		//total_volume = std::accumulate(particles.VolRef().begin(), particles.VolRef().end(), 0.0);
		//Fix_Render_Height();

		std::cout << "Total Number of Points: " << particles.Size() << std::endl;
	}

	//IB functions
	void Prepare_Meta_IB_Advance(void);
	real IB_Pressure_Gradient(int idx, real dt);//along normal

	//Particle geometry functions
	//calculate H(for thick force), SA(for operators)
	void Update_Particle_Heights(void);
	void Update_Max_Velocity(void);
	bool In_Calc_Field(int idx, const OperatorParams& params)const;

	MatrixD Build_Metric(const VectorD& t1_unprojected, const VectorD& unit_norm);
	VectorD Nominal_Vector(const VectorD& center, const VectorD& center_norm, const VectorD& pos, const VectorD& pos_norm, const VectorD& vec);//vec at pos -> nominal vec at center
	real Surface_Kernel_Weight(const real& norm, const KernelType kernel_type)const;
	VectorT Surface_Kernel_Grad(const VectorT& surface_vec, const KernelType kernel_type)const;
	VectorD Surface_Vector_To_World(int idx, const VectorT& v_surface);
	VectorT World_Vector_To_Surface(int idx, const VectorD& v_world)const;
	int Nearest_Boundary_Neighbor(int idx);//-1 for none
	int Nearest_Fluid_Neighbor(int idx);
	VectorD Relative_Vector_World(int p, int q) const { return particles.X(q) - particles.X(p); }//p->q
	//p->q, local frame of p
	VectorT Relative_Vector_Surface(int p, int q) const { return World_Vector_To_Surface(p, Relative_Vector_World(p, q)); }
	real Distance_World(int p, int q) { return Relative_Vector_World(p, q).norm(); }
	int Tangential_Neighbor_Size(int idx) { return (int)surface.Tangential_Neighbor_Of(idx).size(); }
	int Opposite_Side_Nb_Num(int idx, int ref_idx);

	real Compute_Unsampled_Area(int i, real tube_r);

	real Compute_Enclosed_Volume();

	void Seed_Droplet(void);
	void Update_Catenoid_Boundary(real dt);

	int Merge_3d_Particles(void);	////return how many particles are merged
	bool Is_Solitary(int idx, const std::string &mode)const;

	void Mark_Idle_Particles(Array<int>& is_solitary);//mark the particles that are too close to the boundary 
	void Mark_Solitary_Particles(Array<int>& is_solitary);//does not clear is_solitary
	void Mark_Fast_Particles(Array<int>& is_fast, real vel_threshold);//does not clear is_fast
	void Transform_3D_Particles(const Array<int>& to_transform);

	void Apply_3D_Surface_Tension(void);//apply surface tension from 3d particles

	real Cohesion_Kernel(const real& r, const real& h) const {
		real alpha = 32.0 / (pi * pow(h, 9));
		if (0 <= r && 2 * r <= h) {//0~0.5h
			return alpha * 2 * pow((h - r) * r, 3) - pow(h, 6) / 64.0;
		}
		else if (2 * r > h && r <= h) {
			return alpha * pow((h - r) * r, 3);
		}
		else return 0;
	}

	void Apply_Multiphase_Cohesion(void);//apply surface tension from 3d particles

	real Rand_Number(void) const{
		real number;
#pragma omp critical
		{
			number = (real)rand() / ((real)RAND_MAX + 1);
		}
		return number;
	}

	real Surface_Tension_Coefficient(real conc) {
		//return water_st - film_elasticity * particles.Gamma(idx) * particles.ND(idx); --- this is not necessary
		//return particles.Gamma(idx);
		return gamma_water - gamma_soap * conc;
	}

	void Update_Surface_Tension(void) {
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) particles.Gamma(i) = Surface_Tension_Coefficient(particles.Conc(i));
	}

	void Scheduled_Explosion(Array<int>& is_solitary, real current_time) {
		if (current_time >= init_explosion_time && !has_init_explosion) {
#pragma omp parallel for
			for (int i = 0; i < to_explode.size(); i++) {
				int idx = to_explode[i];
				//particles.V(idx) = -0. * VectorD::Unit(0);
				particles.V(idx) = -.2 * VectorD::Unit(0);
			}
			has_init_explosion = true;
			std::cout << "Begin Explotion Procedur, set velocity" << std::endl;
			return;
		}
		if (has_exploded) {
			if (has_leaked) return;
			if (solitary_num > 0) {
				//n_pressure_params.closed = true;
				//n_pressure_params.pressure_decay = 0.1;
				//n_pressure_params.capillary_coeff *= 0.015;
				has_leaked = true;
				for (int i = 0; i < particles.Size(); i++) {
					particles.V(i) *= 1./10.; //slow down time 50 times
				}
				std::cout << "Propagation Begins!" << std::endl;
			}
			return;
		}
		if (current_time < explode_time) return;
#pragma omp parallel for
		for (int i = 0; i < to_explode.size(); i++) {
			int idx = to_explode[i];
			//particles.V(idx) = -0. * VectorD::Unit(0);
			particles.V(idx) *= -0.2;
			is_solitary[idx] = true;
		}
		has_exploded = true;
		std::cout << "Explosion Begins" << std::endl;
		n_pressure_params.closed = false;
		n_pressure_params.capillary_coeff = n_pressure_params.init_capillary_coeff * 0.002;
		//n_pressure_params.closed = false;
		//n_pressure_params.closed = true;
	}


	//use poly6 kernel
	//for f==1, pass a nullptr
	real Surface_Weighted_Sum(int i, std::function<real(const int)> f, const OperatorParams& params, MirrorCompensator<d>& mcp);
	//use spiky kernel
	VectorT Surface_Gradient_Symmetric(int i, std::function<real(const int)> f, const OperatorParams& params, bool force_symmetric = false);
	VectorT Surface_Gradient_Difference(int i, std::function<real(const int)> f, const OperatorParams& params);
	real Surface_Divergence(int i, std::function<VectorD(const int)> f, const OperatorParams& params);
	template<class T> T Surface_Laplacian(int i, std::function<T(const int)> f, const OperatorParams& params, bool print_things = false)const;
	VectorD Gravity_Adhesion(int i, const OperatorParams& params);

	//// World forces
	VectorD Boundary_Force(int i);

	void Diffuse_Scalar(Array<real>& f, const real& dt, const real& coeff, const OperatorParams& params);

	//// Update of forces
	void Update_Fluid_Pressure(void);
	//world forces
	void Update_Gravity_Forces(void);
	void Update_Friction_Forces(void);
	void Update_External_Forces(void);
	void Update_Boundary_Forces(void);
	void Update_Viscosity_Forces(void);
	//tangential forces
	void Update_Vorticity_Confinement_Forces(void);
	void Update_Tangential_Pressure_Forces(void);
	void Update_Marangoni_Forces(void);
	//normal forces
	void Update_Normal_Pressure_Forces(real dt);
	//recompute
	void Update_Divergence(void);//it's nominal divergence
	void Update_Render_Height(real dt);
	void Fix_Render_Height(void);
	void Recompute_Velocities(real dt);
	void Recompute_Geometries(void);


	void Correct_Velocity_With_Analytical_Boundary(int i, real dt);
	void Enforce_Boundary(int i, real dt);
	void Replenish_Boundary_Points_To(SPHBubbleParticles<d>& added_pts, real farthest_dist = 0.1, real replenish_rate = 0.5);//do not clear added_pts
	void Update_Metric_Tensor(void);

	void Print_Iter_Info(real dt) { std::cout << "\niter = " << iter << " dt= " << dt << " length scale= " << length_scale << " max vel= " << max_vel << " cfl= " << max_vel * dt / length_scale << "\n"; }
	void Print_Param_Diagnosis(void) {
		std::cout << "=======       Begin      ======" << std::endl;
		std::cout << "=======   Force Details  ======" << std::endl;
		std::cout << "Avg Tang Pressure Force: " << AuxFunc::Mean(total_pressure_force) << std::endl;
		std::cout << "Avg Viscosity Force: " << AuxFunc::Mean(total_vis_force) << std::endl;
		std::cout << "======= Pressure Fprce Details ======" << std::endl;		
		std::cout << "------     Breakdown     ------" << std::endl;
		std::cout << "Height Pressure Force: " << AuxFunc::Mean(total_height_pressure_force) << std::endl;
		std::cout << "Div Pressure Force: " << AuxFunc::Mean(total_div_pressure_force) << std::endl;
		std::cout << "Height Laplacian Pre. Force: " << AuxFunc::Mean(total_height_laplacian_pressure_force) << std::endl;
		std::cout << "Boundary Pressure Force: " << AuxFunc::Mean(total_bnd_pressure_force) << std::endl;
		std::cout << "=======        End       ======" << std::endl;
		std::cout << "\n" << std::endl;
	};
	void Print_Particle_Dynamics(int i);
	void Check_Neighbor_Num(void);
	void Numerical_Check_SA(void);
	void Numerical_Check_Force(void);
	void Numerical_Check(void);
	void Resize_Debug_Structure(void);
	void Update_Phis(void);

	virtual void Advance(const real dt, const real current_time)
	{
		//if (current_frame > 145) g *= 0.;
		if (diagnosis) Print_Iter_Info(dt);
		if (diagnosis) Resize_Debug_Structure();

		Timer<real> timer;
		timer.Reset();

		if (n_pressure_params.air_pressure_mode == "ib") {
			assert(false && "FluidSPHBubble::Advance: n_pressure_params.air_pressure_mode == ib is not supported in open source version");
			// Prepare_Meta_IB_Advance();
			// air_solver.Advance(dt);
		}

		if (n_pressure_params.pressure_decay > 0.) {
			n_pressure_params.decayed_pressure = n_pressure_params.pressure_constant * pow(1 - n_pressure_params.pressure_decay, 50. * (current_time - explode_time));
			std::cout << "how much left " << pow(1 - n_pressure_params.pressure_decay, 50. * (current_time - explode_time)) << std::endl;
			//params.pressure_constant *= (1 - params.pressure_decay);
		}

		if (interp_g) {
			real progress = (real)(current_frame - first_frame) / (real)(last_frame - first_frame);
			g = progress * g_final + (1 - progress) * g_init;
		}

		if (fluid_3d != nullptr) {
			fluid_3d->Advance(dt, current_time);
			if (merge_3d){
				int merged_num = Merge_3d_Particles();
				if (merged_num > 0) {
					std::cout << "merge " << merged_num << " 3d particles to 2d\n";
					Recompute_Geometries();
				}
			}
		}

		const int pn = particles.Size(); // pn : number of particles

		if (seeding_droplet) Seed_Droplet();
		Update_Surface_Tension();
		Update_Divergence();
		Update_Fluid_Pressure();
		Update_Render_Height(dt);

		if (timing) timer.Elapse_And_Output("surface tension & pressure update\n");

		//world forces
		Update_Gravity_Forces();
		Update_Friction_Forces();
		Update_External_Forces();

		if (timing) timer.Elapse_And_Output("world force update\n");

		//tangential forces
		Update_Vorticity_Confinement_Forces();
		Update_Tangential_Pressure_Forces();
		Update_Marangoni_Forces();

		if (timing) timer.Elapse_And_Output("tangential force update\n");

		//normal forces
		Update_Normal_Pressure_Forces(dt);
		if (diffuse_soap) Diffuse_Scalar(particles.ConcRef(), dt, 0.002, viscosity_params);
		if (vorticity_params.enable_diffuse) Diffuse_Scalar(particles.VrtRef(), dt, vorticity_params.vorticity_diffuse_coeff, viscosity_params);
		
		if (timing) timer.Elapse_And_Output("normal force update\n");

		if (exp_mode == "bursting_bubble") Apply_3D_Surface_Tension();
		if (use_multiphase) Apply_Multiphase_Cohesion();

		Recompute_Velocities(dt);

		if (timing) timer.Elapse_And_Output("recompute velocity-1st\n");

		Update_Boundary_Forces();
		Update_Viscosity_Forces();

		if (timing) timer.Elapse_And_Output("boundary and viscosity forces update\n");

		Recompute_Velocities(dt);

		if (timing) timer.Elapse_And_Output("recompute velocity-2nd\n");

#pragma omp parallel for
		for (int i = 0; i < pn; i++) {
			Enforce_Boundary(i, dt);
			if (particles.Is_Boundary(i))continue;
			particles.X(i) += particles.V(i) * dt;
		}
		if (exp_mode == "catenoid") Update_Catenoid_Boundary(dt);

		//delete and add particles. do not insert into these lines
		static Array<int> to_delete; to_delete.resize(particles.Size()); AuxFunc::Fill(to_delete, 0);
		if (delete_solitary) Mark_Solitary_Particles(to_delete);
		if (delete_idle) Mark_Idle_Particles(to_delete);
		if (exp_mode == "bursting_bubble") Scheduled_Explosion(to_delete, current_time);
		if (fluid_3d != nullptr) Transform_3D_Particles(to_delete);
		if (delete_speedy_particles) Mark_Fast_Particles(to_delete, 100.0);
		added_particles.Resize(0);
		static real last_replenish_time = 0.0;
		if (boundary_params.replenish && current_time >= last_replenish_time + boundary_params.replenish_interval) {
			Replenish_Boundary_Points_To(added_particles, boundary_params.replenish_dx_num * length_scale, boundary_params.replenish_proportion);
			last_replenish_time = current_time;
		}
		particles.Delete_Elements(to_delete);
		particles.Join(added_particles);

		Recompute_Geometries();

		if (timing) timer.Elapse_And_Output("recompute geometries\n");

		if (dynamic_seed_vortex) {
			if (current_frame != last_seed_frame) {
				if (Rand_Number() < dynamic_seed_rate) {
					std::cout << "dynamically seeding a vortex" << std::endl;
					Seed_Vortex_Rand(1, 0.000001, 0.2 * simulation_scale);
				}
				last_seed_frame = current_frame;
			}
		}

		if (diagnosis)Print_Statistics();
		if (diagnosis) timer.Elapse_And_Output_And_Reset("iter");
		Update_Max_Velocity();
		Numerical_Check_SA();
		Numerical_Check();
		iter++;

		if (timing) timer.Elapse_And_Output("post-process\n");
		if (diagnosis) { std::cout << diagnosis << std::endl; Print_Param_Diagnosis(); };

		/*int half_idx = particles.Size() / 2;
		int half_sgn = (particles.V(half_idx).dot(particles.Normal(half_idx)) > 0) ? 1 : -1;
		int knot_idx = -1;
		for (int i = half_idx; i < particles.Size(); i++) {
			int i_sgn = (particles.V(i).dot(particles.Normal(i)) > 0) ? 1 : -1;
			if (i_sgn != half_sgn) {
				knot_idx = i;
				break;
			}
		}
		real pass_omega = pi / 2 - atan2(particles.X(knot_idx)[1], particles.X(knot_idx)[0]);
		real pass_len = pass_omega * simulation_scale;
		output_file << "frame " << current_frame << " iter " << iter << " time " << current_time << " idx " << knot_idx << " len " << pass_len << "\n";*/
	}

	real Center_Smooth_Kernel(real r, real h, real mtp = 2) {//apply in [-h,h], it's an even function
		return std::min(cos(r / h * pi / 2.0) * mtp, 1.0);
	}

	void Seed_Vortex_Rand(int seed_num, real max_strength, real sample_r)
	{
		static std::shared_ptr<RandomInt> rand_id = std::make_shared<RandomInt>(0, particles.Size() - 1);
		static std::shared_ptr<RandomNumber> rand_w = std::make_shared<RandomNumber>(-1, 1);
#pragma omp parallel for
		for (int sid = 0; sid < seed_num; sid++) {
			real h = (0.8 * Rand_Number() + 0.2) * sample_r;
			int id = rand_id->Value();
			const real vorticity = max_strength * rand_w->Value();
			const VectorD& pos = particles.X(id);
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				real r = (particles.X(i) - pos).norm();
				real w = Center_Smooth_Kernel(r, h, sqrt(2.0) + 0.3 * rand_w->Value());
				particles.Vrt(i) += w * vorticity;
			}
		}
	}

	void Update_Sphere_Local_Frame()
	{
		using namespace AuxFunc;
		int pn=particles.Size();
		for(int i=0;i<particles.Size();i++){
			VectorD normal=particles.X(i);
			normal=normal.normalized();
			VectorD t1=-Orthogonal_Vector(normal).normalized();
			
			particles.E(i).col(0)=t1;
			if constexpr (d==2){
				particles.E(i).col(1)=normal;}
			else if constexpr (d==3){
				VectorD t2=(t1.cross(normal)).normalized();
				particles.E(i).col(1)=t2;
				particles.E(i).col(2)=normal;}
		}
	}

	////Verbose prints
	void Print_Params()
	{
		std::cout<<"length_scale: "<<length_scale<<std::endl;
	}

	void Print_Statistics()
	{
		const int pn=particles.Size();

		int avg_nb_num=0;
		int max_nb_num=-1;
		int min_nb_num=std::numeric_limits<int>::max();

		for (int i = 0; i < pn; i++) {
			const auto& now_nbs = surface.Tangential_Neighbor_Of(i);
			int i_nb_num = (int)now_nbs.size();
			avg_nb_num += i_nb_num;
			if (i_nb_num > max_nb_num)max_nb_num = i_nb_num;
			if (i_nb_num < min_nb_num)min_nb_num = i_nb_num;
		}

		std::cout<<"nbs_num avg: "<<avg_nb_num/pn<<", max: "<<max_nb_num<<", min: "<<min_nb_num<<std::endl;
	}
};

template<int d>
template<class T>
inline T FluidSPHBubble<d>::Surface_Laplacian(int i, std::function<T(const int)> f, const OperatorParams& params, bool verbose)const
{
	const auto& nbs = surface.Tangential_Neighbor_Of(i);
	T lap = f(i) - f(i);
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		if (i == j) continue;
		if (!In_Calc_Field(j, params)) continue;
		real S_j = particles.SA(j);
		if (particles.Is_Boundary(j)) { S_j = particles.SA(i); }
		VectorT lr_ji = -Relative_Vector_Surface(i, j);
		real norm_ji = std::max(lr_ji.norm(), length_scale * 1e-1);
		VectorT grad_W = Surface_Kernel_Grad(lr_ji, params.kernel);
		lap += S_j * (f(j) - f(i)) * 2 * grad_W.norm() / norm_ji;
		//if (verbose&&i == 3) { std::cout << "i=" << i << ",j=" << j << ",term=" << S_j << "*" << (f(j) - f(i)) << "*" << 2 * grad_W.norm() / lr_ji.norm() << "\n"; }
	}
	return lap;
}


#endif



