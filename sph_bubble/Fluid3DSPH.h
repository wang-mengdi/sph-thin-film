#ifndef __Fluid3DSPH_h__
#define __Fluid3DSPH_h__
#include "FluidSPH.h"
#include "AnalyticalBoundary.h"

template<int d>
class Fluid3DSPH : public FluidSPH<d, 1024> {
	Typedef_VectorDii(d); using Base = FluidSPH<d, 1024>;
public:
	//using Base::particles;
	using Base::nbs;
	using Base::Update_Neighbors;
	using Base::spatial_hashing;
	using Base::kernel;
	using Base::Collision;
	using Base::Print_Statistics;
	using Base::Print_Params;

	using Base::den_0;			////rest density	
	using Base::nden_0;		////rest number density
	using Base::V_0;			////rest particle volume
	using Base::mass_0;		////rest particle mass
	using Base::avg_nb_num;		////averaged #particles
	using Base::h;				////supporting radius			
	using Base::kp;			////density-pressure coefficient
	using Base::vis;			////viscosity
	using Base::st;			////surface tension

	using Base::use_body_force;
	using Base::use_fixed_dt;
	using Base::use_central_gravity;
	using Base::verbose;
public:
	real max_vel=0;
	real length_scale = 1.;
	real cohesion_coeff = 0.1;
	real curvature_coeff = 0.1;
	real sn_coeff = 1.;
	VectorD g = VectorD::Zero();
	std::shared_ptr<NeighborSearcher<d> > nbs_searcher;
	bool use_kd_tree_nb_search = true;
	Array<Array<int> > vol_nbs;
	real v_r;
	real simulation_scale;
	SPHBubbleParticles<d> particles;
	bool use_local_gravity = false;
	bool prune_far_away = false;
	bool prune_low = false;
	real low_bound = 0.0;
	bool use_surface_tension = false;
	real speed_decay = 0.999;
	bool clip_velocity = true;
	real vel_threshold = 0.1;
	Array<real> cs;
	Array<real> curvatures;
	real pressure_multiplier = 1.;
	AnalyticalBoundary<d> analytical_boundary;
	real curvature_multiplier = 1.;
	real viscosity_multiplier = 1.;
	real cohesion_multiplier = 1.;
	real adhesion_multiplier = 0.00002;


public:

	void Initialize(const real _length_scale, const real _rho0 = (real)1000., const real _simulation_scale = 0.1, const VectorD gravity = VectorD::Zero())
	{
		max_vel = 0.1;
		length_scale = _length_scale;
		den_0 = _rho0;
		simulation_scale = _simulation_scale;
		g = gravity;

		kp = (real)1e3;
	//	vis = den_0 * (real)1e-2;
		vis = 1e1;
		st = (real)1;
		
		v_r = 9 * length_scale;
		if (use_kd_tree_nb_search) nbs_searcher = std::make_shared<NeighborKDTree<d> >();
		else nbs_searcher = std::make_shared<NeighborHashing<d, 1024> >(v_r);

		kernel = std::make_shared<KernelSPH<d> >(v_r);
		spatial_hashing.Initialize(v_r);

		////update nden_0
		nbs.resize(particles.Size());
		vol_nbs.resize(particles.Size());
		//Update_Neighbors();
		Update_Neighbors_2();
		int pn = particles.Size();
		nden_0 = (real)0;
		for (int i = 0;i < pn;i++) {
			real nden = (real)0;
			int nb_n = vol_nbs[i].size();

			for (int k = 0; k < nb_n; k++) {
				int j = vol_nbs[i][k];
				real dis_ij = (particles.X(i) - particles.X(j)).norm();
				nden += kernel->W_Poly6(dis_ij);
			}
			nden_0 += nden;
		}
		nden_0 /= (real)pn;
		std::cout << "Real nden_0: " << nden_0 << " BUT WE SET IT TO 1 HERE\n";
		nden_0 = 1;

		if (verbose)Print_Params();
	}

	real Cohesion_Kernel(const real &r, const real &h) const{
		real alpha = 32.0 / (pi * pow(h, 9));
		if (0 <= r && 2 * r <= h) {//0~0.5h
			return alpha * 2 * pow((h - r) * r, 3) - pow(h, 6) / 64.0;
		}
		else if (2 * r > h && r <= h) {
			return alpha * pow((h - r) * r, 3);
		}
		else return 0;
	}

	real Adhesion_Kernel(const real& r, const real& h) const {
		real alpha = 0;
		if (r < h && 2 * r > h) {//0~0.5h
			alpha = 0.007 / pow(h, 3.25) * pow(-4*pow(r,2)/h + 6*r -2*h, 1. / 4.);
		}
		return alpha;
	}

	void Update_Neighbors_2()
	{
		//update data structure nbs_searcher
		std::function<bool(const int)> valid = [&](const int idx)->bool {return true; };
		nbs_searcher->Update_Points(particles.XRef(), valid);
		//use updated nbs_searcher to update tang_nbs
		int pn = particles.Size();
		vol_nbs.resize(pn);
#pragma omp parallel for
		for (int i = 0; i < pn; i++) {
			const VectorD& pos = particles.X(i);
			nbs_searcher->Find_Neighbors(pos, v_r, vol_nbs[i]);
		}
	}

	void Add_Particle(const VectorD& X, const VectorD& V, const real& m, const real& rh) {
		int k = particles.Add_Element();
		particles.X(k) = X;
		particles.V(k) = V;
		particles.F(k) = VectorD::Zero();
		particles.M(k) = m;
		particles.RH(k) = rh;
	}


	void Correct_Velocity_With_Analytical_Boundary(int i, real dt)
	{
		if (!analytical_boundary.Available())return;
		real radius = length_scale * 0.5;
		real phi; VectorD normal;
		bool flg = analytical_boundary.Get_Nearest_Boundary(particles.X(i), phi, normal);
		if (phi < radius) {
			
			VectorD tangential_velocity = particles.V(i) - particles.V(i).dot(normal) * normal;
			if (phi > 0) {//in "normal" fluid region
				particles.V(i) = tangential_velocity;
			}
			else {
				real normal_rate = (-phi + length_scale * 0.5) / dt;
				particles.V(i) = tangential_velocity + normal_rate * normal;
			}
		}
	}

	virtual void Advance(const real dt, const real time = 0)
	{	
		//std::cout << "advance 3d with particle size: " << particles.Size() << std::endl;
		//////update nbs
		nbs.resize(particles.Size());
		vol_nbs.resize(particles.Size());
		cs.resize(particles.Size());
		curvatures.resize(particles.Size());
		//Update_Neighbors();
		Update_Neighbors_2();
		const int pn = particles.Size();

		////update number density and pressure
#pragma omp parallel for
		for (int i = 0; i < pn; i++) {
			real nden = (real)0;
			int nb_n = vol_nbs[i].size();
			for (int k = 0; k < nb_n; k++) {
				int j = vol_nbs[i][k];
				real dis_ij = (particles.X(i) - particles.X(j)).norm();
				nden += kernel->W_Poly6(dis_ij);
			}

			particles.Vol(i) = (real)1 / nden;
			//particles.P(i) = nden / nden_0 - 1.0;
			particles.P(i) = pow(nden / nden_0, 5) - 1.0;
		}

#pragma omp parallel for
		// update normal -->gradient of color field
		for (int i = 0; i < pn; i++)
		{
			int nb_n = vol_nbs[i].size();
			VectorD sn_i = VectorD::Zero();
			for (int k = 0; k < nb_n; k++) {
				int j = vol_nbs[i][k];
				VectorD r_ij = particles.X(i) - particles.X(j);
				//sn_i += cs[j] * particles.Vol(j) * kernel->Grad_Poly6(r_ij);
				sn_i += particles.Vol(j) * kernel->Grad_Quintic(r_ij);
			}
			particles.SN(i) = sn_i;
		}
#pragma omp parallel for
		// update curvature -->laplacian of color field
		for (int i = 0; i < pn; i++)
		{
			int nb_n = vol_nbs[i].size();
			real curvature_i = 0;
			for (int k = 0; k < nb_n; k++) {
				int j = vol_nbs[i][k];
				if (i == j) continue;
				VectorD r_ij = particles.X(i) - particles.X(j);
				real r_ij_norm = std::max(length_scale * 0.1, r_ij.norm());
				curvature_i += -(cs[i]-cs[j]) * particles.Vol(j) * 2*kernel->Grad_Poly6(r_ij).norm()/r_ij_norm;
			}
			curvatures[i] = curvature_i;
		}

		// compute geometric center of points
		VectorD center = AuxFunc::Mean(particles.XRef());

#pragma omp parallel for
		//////update forces
		for (int i = 0; i < pn; i++) {
			real one_over_m = 1./particles.M(i);
			VectorD aggregated_f_p = VectorD::Zero();
			VectorD aggregated_f_v = VectorD::Zero();
			VectorD aggregated_f_cohesion = VectorD::Zero();
			VectorD aggregated_f_curvature = VectorD::Zero();
			int nb_n = vol_nbs[i].size();
			for (int k = 0; k < nb_n; k++) {
				int j = vol_nbs[i][k];
				if (particles.Is_Boundary(j)) continue;
				VectorD r_ij = particles.X(i) - particles.X(j);
				real r2 = r_ij.squaredNorm();
				real one_over_r2 = (r2 == (real)0 ? (real)0 : (real)1 / r2);
				VectorD v_ij = particles.V(i) - particles.V(j);
				VectorD f_p = -(particles.P(i) * pow(particles.Vol(i), 2) + particles.P(j) * pow(particles.Vol(j), 2)) * kernel->Grad_Spiky(r_ij) * particles.M(i);
				VectorD f_v = vis * particles.Vol(i) * particles.Vol(j) * v_ij * one_over_r2 * r_ij.dot(kernel->Grad_Spiky(r_ij))*particles.M(i);
				//Surface Tension model
				VectorD f_cohesion;
				//f_cohesion = 1 / (particles.P(i) + particles.P(j)) * -1 * particles.M(i) * particles.M(j) * r_ij.normalized() * 0.15e2 * 0.7e7 * 5 * 4000;
				f_cohesion = -1 * particles.M(i) * r_ij.normalized() * Cohesion_Kernel(r_ij.norm(), v_r) * 1e-5;

				aggregated_f_p += f_p;
				aggregated_f_v += 1000 * f_v;
				aggregated_f_cohesion +=  f_cohesion;
				aggregated_f_curvature += -1 * particles.M(i) * (particles.SN(i) - particles.SN(j)) * 1e-4;
			}
			
			VectorD aggregated_f_adhesion = VectorD::Zero();
			for (int k = 0; k < nb_n; k++) {
				int j = vol_nbs[i][k];
				if (!particles.Is_Boundary(j)) continue;
				VectorD r_ij = particles.X(i) - particles.X(j);
				aggregated_f_adhesion += -adhesion_multiplier * particles.M(i) * nden_0 * particles.Vol(j) * Adhesion_Kernel(r_ij.norm(), v_r) * r_ij.normalized();
			}
			VectorD f_st;
			//muller way
			//if (particles.SN(i).norm()> 1e-8) f_st = .33e-11 * 1 * curvatures[i] * particles.SN(i)/ particles.SN(i).norm();

			//f_st = (1*2e8 * aggregated_f_cohesion + 1e3*aggregated_f_curvature);
			f_st = cohesion_multiplier * aggregated_f_cohesion + curvature_multiplier * aggregated_f_curvature;
			
			 
			VectorD a_g = VectorD::Zero();
			if (use_central_gravity) {
				a_g = 3e-5 * (center - particles.X(i));
				//a_g_2 = f_st;
				aggregated_f_p *= 7e-9;
				aggregated_f_v *= 6e-8;
			}
			else if (use_local_gravity) {
				int nbs_num = 1;
				VectorD local_center = particles.X(i);
				for (int k = 0; k < nb_n; k++) {
					int j = vol_nbs[i][k];
					if ((particles.X(i)-particles.X(j)).norm() < 2*length_scale){
						local_center += particles.X(j);
						nbs_num += 1;
					}
				}
				local_center *= (1. / (real)nbs_num);
				a_g = 6e-7 * (local_center - particles.X(i));
				aggregated_f_p *= 7e-11;
				aggregated_f_v *= 6e-10;
			}
			else {
				a_g = VectorD::Zero();
				//aggregated_f_p *= 7e-10; // this works for changing box to sphere
				/*aggregated_f_p *= 7e-10;
				aggregated_f_v *= 1e-8;*/
			}

			aggregated_f_p *= pressure_multiplier;

			particles.F(i) = one_over_m * (a_g + aggregated_f_p + viscosity_multiplier * aggregated_f_v + aggregated_f_adhesion);
			if (use_surface_tension) {
				particles.F(i) += one_over_m * f_st;
			}
			
			//particles.F(i) *= 0.1;

			//if (i == 955) { 
			//	std::cout << "x is:\n" << particles.X(i) << std::endl;
			//	std::cout << "a_g is: \n" << a_g << std::endl; 
			//	std::cout << "f_st is: \n" << a_g_2 << std::endl;
			//	std::cout << "cohesion is: \n" << aggregated_f_cohesion << std::endl;
			//	std::cout << "curvature is: \n" << aggregated_f_curvature << std::endl;
			//	std::cout << "f p is: \n" << aggregated_f_p << std::endl;
			//	std::cout << "f v is: \n" << aggregated_f_v << std::endl;
			//}

			particles.F(i) += g;
		}

#pragma omp parallel for
		////time integration
		for (int i = 0; i < pn; i++) {
			particles.V(i) *= speed_decay;
			VectorD increment;
			increment = particles.F(i) * dt;
			if (clip_velocity) {
				real vel_norm = increment.norm();
				if (vel_norm > vel_threshold) {
					increment = increment * vel_threshold / vel_norm;
				}
			}
			particles.V(i) += increment;
			Correct_Velocity_With_Analytical_Boundary(i, dt);
			if (particles.Is_Boundary(i)) particles.V(i) *= 0;
			particles.X(i) += particles.V(i) * dt;
		}

		max_vel = 0.;
		for (int i = 0; i < pn; i++) {
			//update max vel
			if (particles.V(i).norm() > max_vel) {
				max_vel = particles.V(i).norm();
			}
		}

		static Array<int> to_delete; to_delete.resize(particles.Size()); AuxFunc::Fill(to_delete, 0);
		static std::shared_ptr<RandomNumber> rand_w = std::make_shared<RandomNumber>(0, 1);
		if (prune_far_away) {
#pragma omp parallel for
			for (int i = 0; i < pn; i++) {
				if ((particles.X(i)).norm() > 0.201 * simulation_scale) {
					if(rand_w->Value() < 0.1) to_delete[i] = true;
				}
			}
		}
		if (prune_low) {
#pragma omp parallel for
			for (int i = 0; i < pn; i++) {
				if (particles.X(i)[1] < low_bound) {
					to_delete[i] = true;
				}
			}
		}
		particles.Delete_Elements(to_delete);
		////collision
		//if (Collision != nullptr)Collision(dt);

		if (verbose)Print_Statistics();
	}
};


#endif
