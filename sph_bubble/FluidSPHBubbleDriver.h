//////////////////////////////////////////////////////////////////////////
// SPH bubble fluid driver
// Copyright (c) (2018-), Xiangxin Kong
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __FluidSPHBubbleDriver_h__
#define __FluidSPHBubbleDriver_h__
#include "Driver.h"
#include "FluidSPHBubble.h"
#include "PointSetFunc.h"
#include "TinyObjLoader.h"
#include "AuxFunc.h"
#include "PerlinNoise.hpp"
#include "Fluid3DSPH.h"

#include "omp.h"

template<int d> class FluidSPHBubbleDriver : public Driver
{
	Typedef_VectorDii(d);
	Typedef_MatrixD(d);
	using Base=Driver;
public:
	using VectorT=Vector<real,d-1>;								////tangential vector type
	using VectorTi=Vector<int,d-1>;								////tangential vector int

	FluidSPHBubble<d> fluid;
	Fluid3DSPH<d> fluid_3d;

	std::string scene_file_name = "";
	std::string init_snapshot_name = "";
	bool save_all = false;
	real default_mass = 1.0;
	real total_mass = 1000.;
	

	virtual void Advance_To_Target_Time(const real target_time)
	{
		bool done = false; 
		for (int substep = 1; !done; substep++) {
			real vel = std::max(CFL() * fluid.length_scale * frame_rate, fluid.max_vel);
			vel = std::max(vel, fluid_3d.max_vel);
			real dt = CFL() * fluid.length_scale / vel;
			//std::cout << "CFL? " << CFL() << std::endl;
			//std::cout << "fluid.length_scale " << fluid.length_scale << std::endl;
			//std::cout << "frame_rate " << frame_rate << std::endl;;
			//std::cout << "fluid.max_vel " << fluid.max_vel << std::endl;
			//std::cout << "Vel " << vel << std::endl;
			if (max_iter_per_frame > 0) {
				dt = std::max(dt, 1.0 / (frame_rate * max_iter_per_frame));
			}
			if (time + dt >= target_time) { dt = target_time - time; done = true; }
			else if (time + 2 * dt >= target_time) { dt = (real).5 * (target_time - time); }
			Advance_One_Time_Step(dt, time);
			time += dt;
		}
	}

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		//std::cout << "this is called" << std::endl;
		double begin_time = omp_get_wtime();
		//std::cout << "Frame rate: " << frame_rate << std::endl;
		//std::cout << "dt: " << dt << std::endl;
		if (verbose) std::cout << "[" << std::setw(6) << (int)(1.0 / frame_rate / dt + 0.5) << " Iterations Per Frame] ";
		std::cout << std::defaultfloat;

		fluid.current_frame = current_frame;
		fluid.Advance(dt,time);

		double end_time = omp_get_wtime();
		if (verbose) std::cout << "    ... " << std::setw(7) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << end_time - begin_time << "s used" << std::endl;
		std::cout << std::defaultfloat;
	}

	void Save_Snapshot(const int frame);
	void Load_Snapshot(const int frame);
	void Write_Scalar_Field(std::string file_name, const Array<real>& arr, const real scale);
	void Write_Air_Solver(const int frame);
	virtual void Write_Output_Files(const int frame);
	
	
	virtual void Initialize()
	{
		cfl=.001;
		snapshot_stride = 5;
		switch (test) {
		case 1:	Case_1(); break;//Perlin Noise Test 3D sphere
		case 2: Case_2(); break;//3D sphere
		case 3: Case_3(); break;//obj test
		case 4: Case_4(); break;
		case 5: Case_5(); break;
		case 6: Case_6(); break;
		case 7: Case_7(); break;
		case 8: Case_8(); break;//Capillary wave on a segment
		case 9: Case_9(); break;//Capillary Wave on an enclosed circle
		case 10: Case_10(); break;//Newton ring with dynamic gravity
		case 11: Case_11(); break;//test if a 3d "rod" will merge to a few droplets
		case 12: Case_12(); break;//"blow bubble"
		case 13: Case_13(); break;//a half sphere falling down with gravity
		case 14: Case_14(); break;//a line, for debugging SPH operators
		case 15: Case_15(); break;//circle points with initial tangential velocity toward origin
		case 16: Case_16(); break;//test different kernels
		case 17: Case_17(); break;//random thin film oscillate
		case 18: Case_18(); break;//half sphere
		case 19: Case_19(); break;//test marangoni effect
		case 20: Case_20(); break;//rotating boundary in a plane
		case 21: Case_21(); break;//catenoid
		case 22: Case_22(); break;//hand designated height field, for shader calibration
		case 23: Case_23(); break;//surface water flow around obstacle
		case 24: Case_24(); break;//a rotating film fixed on a ring --2D
		case 25: Case_25(); break;//a rotating film fixed on a ring --3D
		case 26: Case_26(); break;//surface flow in a circular ring
		case 27: Case_27(); break;//a copy of case 26, for numerical tests
		case 28: Case_28(); break;//a copy of case 26, scale up size instead of fineness
		case 29: Case_29(); break;//a copy of case 28, with horizontal gravity
		case 30: Case_30(); break;//Test vorticity on a plane
		case 31: Case_31(); break;//a copy of cas 28, with vertical gravity
		case 32: Case_32(); break;//initial velocity towards center to test tangential response (upgraded case 15)
		case 33: Case_33(); break;//a copy of case 29, used to test scalability of cappilary actions
		case 34: Case_34(); break;//a copy of case 29, but with Default Sim Params
		case 35: Case_35(); break;//test marangoni random motion
		case 36: Case_36(); break;//test RT instability
		case 37: Case_37(); break;//test RT instability 2 (with height pressure)
		case 38: Case_38(); break;//Newton Ring plain
		case 39: Case_39(); break;//used for point seeding (generating initial snapshot)
		case 40: Case_40(); break;//Newton-Ring with Perlin Noise
		case 41: Case_41(); break;//Newton Ring with random initialization
		case 42: Case_42(); break;//RT finger with thick heavy top and thin light bottom
		case 43: Case_43(); break;//Newton Ring with init Vortex
		case 44: Case_44(); break;//case 44 but with thick heavy drips on top
		case 45: Case_45(); break;//test Marangoni
		case 46: Case_46(); break;//Catenoid
		case 47: Case_47(); break;//Incorporate 3d sph
		case 48: Case_48(); break;//bursting bubble
		case 49: Case_49(); break;//test 3d surface tension
		case 50: Case_50(); break;//merge3d into 2d marangoni
		case 51: Case_51(); break;//merge3d into 2d marangoni
		case 52: Case_52(); break;//case 38 newton ring BUT WITH SOME HACKS AND RANDOMIZATIONS
		case 53: Case_53(); break;//test difference in laplacian pressure and height pressure
		case 54: Case_54(); break;
		case 55: Case_55(); break;
		}
		if (fluid.particles.Size() == 0) {
			std::cerr << "FluidSPHBubbleDriver::Initialize error: points not initialized\n";
			exit(0);
		}
		fluid.verbose = verbose;

		if (first_frame > 0) {
			int last_saved = int(first_frame / snapshot_stride) * snapshot_stride;
			std::cout << "Run from snapshot frame " << last_saved << "\n";
			current_frame = last_saved;
			time = Time_At_Frame(current_frame);
			Load_Snapshot(last_saved);
			fluid.Numerical_Check_Force();
			std::cout << "loaded\n";
		}
		if (fluid.analytical_boundary.Available()) {
			//std::cout << "analytical boundary is available!" << std::endl;
			fluid.Update_Phis();
		}
		fluid.surface.Update();

		
		//for (int i = 0;i < fluid.particles.Size();i++) {
		//	auto z_func_i = [&](const int idx)->real {return fluid.particles.X(idx).dot(fluid.particles.Normal(i)); };
		//	fluid.particles.KH(i) = fluid.Surface_Laplacian<real>(i, z_func_i, fluid.n_pressure_params.opt_mode, true);
		//}
	}

	virtual void Initialize_Particle(int i, real mass, real render_thickness)	////Initialize particle with index i
	{
		//total_mass/num_fluid_particles
		fluid.particles.M(i) = mass;
		fluid.particles.F(i)=VectorD::Zero();
		fluid.particles.V(i)=VectorD::Zero();
		//fluid.particles.B(i)=0;
		fluid.particles.Gamma(i) = (real)1;
		//fluid.particles.H0(i) = (real)1;
		fluid.particles.Conc(i) = (real)0.0;
		//fluid.particles.Conc_F(i) = (real)0.0;
		fluid.particles.Vrt(i) = 0.0;

		fluid.particles.Vol(i) = fluid.particles.M(i) / 1e3;
		fluid.particles.RH(i) = render_thickness;
		fluid.particles.Phase(i) = 0;
	}

	std::function<VectorD(const VectorD&)> Corridor_Flow_Func(real R, real line_vel);
	void Seed_Vortex(int seed_num, real max_strength, real sample_r);
	void Seed_Vortex_Rand(int seed_num, real max_strength, real sample_r);

	void Initialize_Center_Oil(void);
	void Set_Fixed_Gamma(real gamma);
	void Set_Physical_Parameters(void);
	real Perlin_Noise(VectorD pos, std::uint32_t seed, real perlin_freq, real perlin_scale, std::uint32_t octaves=4);
	void Case_1(void);
	void Case_2(void);
	void Case_3(void);
	void Case_4(void);
	void Case_5(void);
	void Case_6(void);
	void Case_7(void);
	void Case_8(void);
	void Case_9(void);
	void Case_10(void);
	void Case_11(void);
	void Case_12(void);
	void Case_13(void);
	void Case_14(void);
	void Case_15(void);
	void Case_16(void);
	void Case_17(void);
	void Case_18(void);
	void Case_19(void);
	void Case_20(void);
	void Case_21(void);
	void Case_22(void);
	void Case_23(void);
	void Case_24(void);
	void Case_25(void);
	void Case_26(void);
	void Case_27(void);
	void Case_28(void);
	void Case_29(void);
	void Case_30(void);
	void Case_31(void);
	void Case_32(void);
	void Case_33(void);
	void Case_34(void);
	void Case_35(void);
	void Case_36(void);
	void Case_37(void);
	void Case_38(void);
	void Case_39(void);
	void Case_40(void);
	void Case_41(void);
	void Case_42(void);
	void Case_43(void);
	void Case_44(void);
	void Case_45(void);
	void Case_46(void);
	void Case_47(void);
	void Case_48(void);
	void Case_49(void);
	void Case_50(void);
	void Case_51(void);
	void Case_52(void);
	void Case_53(void);
	void Case_54(void);
	void Case_55(void);

	real Rand_Number(void) {
		real number = (real)rand() / ((real)RAND_MAX + 1);
		return number;
	}

};

	
#endif

