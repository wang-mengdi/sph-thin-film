#include "FluidSPHBubble.h"

template<int d>
void MirrorCompensator<d>::Initialize(FluidSPHBubble<d>* _fluid, int idx, real _weight) {
	weight = _weight;
	fluid = _fluid;
	bnd_nb = fluid->Nearest_Boundary_Neighbor(idx);
	if (bnd_nb != -1) {
		bnd_vec = fluid->Relative_Vector_Surface(idx, bnd_nb);
		bnd_sqrval = bnd_vec.squaredNorm();
	}
	if (bnd_nb == idx) self_bnd = true;
}

template<int d>
real MirrorCompensator<d>::Coeff_Offset(VectorT lr_ij, int order) {
	if (bnd_nb == -1) return 0;
	real val = lr_ij.dot(bnd_vec);
	if (self_bnd) {
		return 0.;
	}
	if (val > bnd_sqrval) {//leave it as it is if it is to the same side as bnd_nb but further away
		return 0;
	}
	else if (val < -bnd_sqrval) {//if in range
		if (order % 2) return -1 * weight;
		else return 1 * weight;
	}
	else {
		return 0.;
	}
}


template<int d>
void FluidSPHBubble<d>::Prepare_Meta_IB_Advance(void)
{
	assert(false && "Prepare_Meta_IB_Advance is not implemented for open source version.");
	// const auto& mac_grid = air_solver.mac_grid;
	// iterate_face(axis, iter, mac_grid) {
	// 	const VectorDi face = iter.Coord();
	// 	VectorD pos = mac_grid.Face_Center(axis, face);
	// 	int idx = surface.nbs_searcher->Find_Nearest_Nb(pos);
	// 	if (idx != -1) {
	// 		real phi = (particles.X(idx) - pos).norm();
	// 		air_solver.alpha(axis, face) = air_solver.IB_Kernel_Fluid(phi);
	// 		air_solver.body_velocity(axis, face) = particles.V(idx)[axis];
	// 	}
	// 	else {
	// 		air_solver.alpha(axis, face) = 1.0;
	// 		air_solver.body_velocity(axis, face) = 0;
	// 	}
	// }
}

template<int d>
real FluidSPHBubble<d>::IB_Pressure_Gradient(int idx, real dt)
{
	assert(false && "IB_Pressure_Gradient is not implemented for open source version.");
	// real dx = air_solver.mac_grid.grid.dx * 0.5;
	// VectorD norm = particles.Normal(idx);
	// VectorD up_pt = particles.X(idx) + norm * dx;
	// VectorD down_pt = particles.X(idx) - norm * dx;
	// real grad_p = (air_solver.Pressure_At(up_pt) - air_solver.Pressure_At(down_pt)) / (2 * dx);
	// //note: pressure in Eularian fluid solver is -dt/rho*p actually
	// grad_p = -grad_p / dt * n_pressure_params.air_density;
	// return grad_p;
}

template<int d>
real FluidSPHBubble<d>::Compute_Unsampled_Area(int i, real tube_r) {
	real total_wasted_area = 0;
	const auto& nbs = surface.Tangential_Neighbor_Of(i);
	std::vector<int> bnd_nbs;
	std::vector<real> polar_angles;
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		if (!particles.Is_Boundary(j))continue;//only boundary neighbors
		if (i == j)continue; //not itself
		VectorT lr_ij = Relative_Vector_Surface(i, j); //i pointing to j
		if (lr_ij.norm() >= tube_r) continue; // only boundary neighbors within the TUBE
		//std::cout << "distance" << lr_ij.norm() << std::endl;
		bnd_nbs.push_back(j);
		polar_angles.push_back(atan2(lr_ij[1], lr_ij[0]));
	}
	//std::cout << "B" << std::endl;
	if (bnd_nbs.size() > 1) { //if has more than one boundary neighbors, calculate unsampled neighborhood area

		std::vector<int> V(polar_angles.size());
		std::iota(V.begin(), V.end(), 0); //Initializing
		sort(V.begin(), V.end(), [&](int u, int v) {return polar_angles[u] < polar_angles[v]; });

		std::vector<int> new_bnd_nbs;
		for (int k = 0; k < polar_angles.size(); k++) {
			new_bnd_nbs.push_back(bnd_nbs[V[k]]);
		}
		bnd_nbs = new_bnd_nbs;

		std::sort(polar_angles.begin(), polar_angles.end());

		polar_angles.push_back(polar_angles[0] + 2 * pi); //append the smallest polar angle + 2 pi to the list

		//compute biggest gap, start from after the gap
		real max_diff = 0;
		int bnd_nb_after_max_diff = 0;
		for (int k = 0; k < polar_angles.size() - 1; k++) {
			real curr_diff = polar_angles[k + 1] - polar_angles[k];
			if (curr_diff > max_diff) {
				max_diff = curr_diff;
				bnd_nb_after_max_diff = k + 1;
			}
		}
		if (bnd_nb_after_max_diff == polar_angles.size() - 1)bnd_nb_after_max_diff = 0; // remember we added in a "fake" point, which is just the same as the first bnd_nb
		std::vector<int> sorted_bnd_nbs = std::vector<int>(bnd_nbs.begin() + bnd_nb_after_max_diff, bnd_nbs.end());
		sorted_bnd_nbs.insert(sorted_bnd_nbs.end(), bnd_nbs.begin(), bnd_nbs.begin() + bnd_nb_after_max_diff);

		real aggregate_angle = 0;
		//loop through the boundary neighbors sorted properly
		for (int k = 0; k < sorted_bnd_nbs.size() - 1; k++) {
			int j1 = sorted_bnd_nbs[k];
			int j2 = sorted_bnd_nbs[k + 1];
			VectorT lr_ij1 = Relative_Vector_Surface(i, j1); //i pointing to j1
			VectorT lr_ij2 = Relative_Vector_Surface(i, j2); //i pointing to j2
			real angle = acos(lr_ij1.dot(lr_ij2) / (lr_ij1.norm() * lr_ij2.norm()));

			aggregate_angle += angle;

			real sector_area = (angle / (2 * pi)) * pi * pow(tube_r, 2);
			VectorT p0, p1, p2;
			p0 = VectorT::Zero();
			p1 = lr_ij1;
			p2 = lr_ij2;
			real triangle_area = abs(0.5 * ((p0[0] * (p1[1] - p2[1])) + (p1[0] * (p2[1] - p0[1])) + (p2[0] * (p0[1] - p1[1]))));

			if (sector_area < triangle_area) {
				std::cout << "this should never happen" << std::endl;
				std::cout << "ij1: " << lr_ij1 << std::endl;
				std::cout << "ij2: " << lr_ij2 << std::endl;
				std::cout << "angle: " << angle << std::endl;
				std::cout << "tri_area: " << triangle_area << std::endl;
				std::cout << "tube r: " << tube_r << std::endl;
				std::cout << "sec_area: " << sector_area << std::endl;
			}
			else {
				total_wasted_area += (sector_area - triangle_area);
			}
		}
	}
	return total_wasted_area;
}

template<int d>
void FluidSPHBubble<d>::Update_Particle_Heights(void)
{
	const BoundaryParams& params = boundary_params;
	if (params.particle_geometry_mode == "mirror_compensate") {
		real tube_r = surface.t_r;
		bool weight_mirror_compensation = false;
		//std::cout << "A" << std::endl;
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			real weight;//weight is used for mirro compensation
			if (!weight_mirror_compensation) weight = default_mirror_weight;
			else if (particles.Is_Boundary(i)) weight = 1;//don't bother for adjusting boundary particles, because they are not mirror compensated anyway
			else {
				//std::cout << "B" << std::endl;
				real total_wasted_area = Compute_Unsampled_Area(i, tube_r);
				if (total_wasted_area < 1e-8) {
					weight = 1;
				}
				else {
					//std::cout << "total wasted area" << total_wasted_area << std::endl;
					int bnd_nb = Nearest_Boundary_Neighbor(i);
					real nearest_nb_dist = Relative_Vector_Surface(i, bnd_nb).norm();
					real half_angle = acos(nearest_nb_dist / tube_r);
					real whole_angle = 2 * half_angle;
					real area1 = (whole_angle / (2 * pi)) * pi * pow(tube_r, 2);
					real area2 = nearest_nb_dist * (tube_r * sin(half_angle));
					real compensated_area = area1 - area2;
					if (compensated_area < 0) {
						std::cout << "compensated area < 0, that's wrong" << std::endl;
					}
					weight = total_wasted_area / compensated_area;
					if (weight < 1) {
						std::cout << i << std::endl;
						std::cout << "total_wasted_area/compensated_area < 1, that's wrong" << std::endl;
						//std::cout << "total_wasted_area" << total_wasted_area << std::endl;
						//std::cout << "aggregate angle:" << aggregate_angle << std::endl;
						//std::cout << "compensation angle:" << whole_angle << std::endl;
						//std::cout << "compensated_area" << compensated_area << std::endl;
					}
				}
			}
			MirrorCompensator<d> mcp; mcp.Initialize(this, i, weight);
			//std::cout << "weight: " << weight << std::endl;
			OperatorParams params = geometry_params; params.calculate_field = "only_fluid";
			real sph_vol_density = Surface_Weighted_Sum(i, Index_Function(particles.VolRef()), params, mcp); //volume per surface area
			if (sph_vol_density == 0) {
				std::cerr << "Update_Particle_Heights error: sph_vol_density==0 of " << i << "\n";
				exit(0);
			}
			particles.SA(i) = particles.Vol(i) / sph_vol_density; //surface area for the 2d sph
			particles.H(i) = particles.Vol(i) / particles.SA(i); //height for
		}
	}
	else if (params.particle_geometry_mode == "no_compensate") {
		for (int i = 0; i < particles.Size(); i++) {
			MirrorCompensator<d> mcp; mcp.Initialize(this, i, 0.0);
			OperatorParams params = geometry_params; params.calculate_field = "all";
			particles.H(i) = Surface_Weighted_Sum(i, Index_Function(particles.VolRef()), params, mcp); //volume per surface area
			particles.SA(i) = particles.Vol(i) / particles.H(i); //surface area for the 2d sph
		}
	}
	else {
		std::cerr << "Update_Particle_Geometries error: undefined geometry mode\n";
	}
}


template<int d>
void FluidSPHBubble<d>::Update_Max_Velocity(void)
{
	max_vel = 0;
	for (int i = 0; i < particles.Size(); i++) {
		if (particles.Is_Boundary(i)) continue;
		if (particles.V(i).norm() > 100) {
			std::cout << "[SPEED WARNING]\n";
			Print_Particle_Dynamics(i);
			if (quit_when_too_fast) exit(0);
		}
		max_vel = std::max(max_vel, particles.V(i).norm());
	}
	if (n_pressure_params.air_pressure_mode == "ib") {
		assert(false && "IB in Update_Max_Velocity is not implemented for open source version.");
		// real nominal_vel = air_solver.max_vel * length_scale / air_solver.mac_grid.grid.dx;
		// max_vel = std::max(max_vel, nominal_vel);
	}
}

template<int d>
bool FluidSPHBubble<d>::In_Calc_Field(int idx, const OperatorParams& params)const
{
	if (params.calculate_field == "only_fluid") {
		return !particles.Is_Boundary(idx);
	}
	else if (params.calculate_field == "only_boundary") {
		return particles.Is_Boundary(idx);
	}
	else if (params.calculate_field == "all") {
		return true;
	}
	else {
		std::cerr << "In_Calc_Field error: unknown calculation field\n";
		return false;
	}
}

template<int d>
Matrix<real, d> FluidSPHBubble<d>::Build_Metric(const VectorD& t1_unprojected, const VectorD& unit_norm)
{
	MatrixD E;
	E.col(0) = AuxFunc::Eliminate_Unit_Component<real,d>(t1_unprojected, unit_norm).normalized();
	E.col(1) = unit_norm;
	if constexpr (d > 2) {
		E.col(2) = E.col(0).cross(E.col(1));
	}
	return E;
}

template<int d>
Vector<real, d> FluidSPHBubble<d>::Nominal_Vector(const VectorD& center, const VectorD& center_norm, const VectorD& pos, const VectorD& pos_norm, const VectorD& vec)
{
	//NOTE: length of center_norm and pos_norm MUST be 1
	if (center == pos) return vec;
	VectorD r_ij = pos - center;
	MatrixD nominal_center = Build_Metric(r_ij, center_norm);
	MatrixD nominal_pos = Build_Metric(r_ij, pos_norm);
	return nominal_center * nominal_pos.transpose() * vec;
}

template<int d>
real FluidSPHBubble<d>::Surface_Kernel_Weight(const real& norm, const KernelType kernel_type)const
{
	if (kernel_type == KernelType::SPIKY) {
		return kernel->W_Spiky(norm);
	}
	else if (kernel_type == KernelType::POLY6) {
		return kernel->W_Poly6(norm);
	}
	else if (kernel_type == KernelType::CUBIC) {
		return kernel->W_Cubic(norm);
	}
	else if (kernel_type == KernelType::GAUSSIAN) {
		return kernel->W_Gaussian(norm);
	}
	else if (kernel_type == KernelType::QUINTIC) {
		return kernel->W_Quintic(norm);
	}
	else {
		std::cerr << "Surface_Kernel_Weight error: unknown type\n";
		assert(false);
		return 0.0;
	}
}

template<int d>
Vector<real, d - 1> FluidSPHBubble<d>::Surface_Kernel_Grad(const VectorT& surface_vec, const KernelType kernel_type)const
{
	if (kernel_type == KernelType::SPIKY) {
		return kernel->Grad_Spiky(surface_vec);
	}
	else if (kernel_type == KernelType::POLY6) {
		return kernel->Grad_Poly6(surface_vec);
	}
	else if (kernel_type == KernelType::CUBIC) {
		return kernel->Grad_Cubic(surface_vec);
	}
	else if (kernel_type == KernelType::GAUSSIAN) {
		return kernel->Grad_Gaussian(surface_vec);
	}
	else if (kernel_type == KernelType::QUINTIC) {
		return kernel->Grad_Quintic(surface_vec);
	}
	else {
		std::cerr << "Surface_Kernel_Grad error: unknown type\n";
		assert(false);
		return VectorT::Zero();
	}
}

template<int d>
Vector<real, d> FluidSPHBubble<d>::Surface_Vector_To_World(int idx, const VectorT& v_surface)
{
	VectorT scaled_v_surface = particles.G(idx).inverse() * v_surface;
	VectorD v_world; surface.Unproject_To_World(scaled_v_surface, particles.E(idx), v_world);
	return v_world;
}

template<int d>
Vector<real, d - 1> FluidSPHBubble<d>::World_Vector_To_Surface(int idx, const VectorD& v_world)const
{
	VectorT v_surface = surface.Project_To_TPlane(v_world, particles.E(idx));
	//surface.Project_To_TPlane(v_world, particles.E(idx), v_surface);
	return v_surface;
}

template<int d>
int FluidSPHBubble<d>::Nearest_Boundary_Neighbor(int idx)
{
	const auto& nbs = surface.Tangential_Neighbor_Of(idx);
	int bnd_nb = -1;
	real sqrlen = 0.0;
	for (int i = 0; i < nbs.size(); i++) {
		int q = nbs[i];
		if (particles.Is_Boundary(q)) {
			real now_sqr = (particles.X(idx) - particles.X(q)).squaredNorm();
			if (bnd_nb == -1 || now_sqr < sqrlen) {
				sqrlen = now_sqr;
				bnd_nb = q;
			}
		}
	}
	return bnd_nb;
}

template<int d>
int FluidSPHBubble<d>::Nearest_Fluid_Neighbor(int idx)
{
	const auto& nbs = surface.Tangential_Neighbor_Of(idx);
	int nearest_nb = -1;
	real sqrlen = 0.0;
	for (int i = 0; i < nbs.size(); i++) {
		int q = nbs[i];
		if (idx == q) continue;
		if (analytical_boundary.Available()) {
			if (phis[q] < 0.05 * length_scale) { //if a fluid neighbor is too close to the boundary, count it out
				continue;
			}
		}
		if (!particles.Is_Boundary(q)) {
			//std::cout << q << " is not boundary" << std::endl;
			real now_sqr = (particles.X(idx) - particles.X(q)).squaredNorm();
			if (nearest_nb == -1 || now_sqr < sqrlen) {
				sqrlen = now_sqr;
				nearest_nb = q;
			}
		}
	}
	//if (idx == 859)std::cout << "nearest_fluid_nb: " << nearest_nb << std::endl;
	return nearest_nb;
}

template<int d>
int FluidSPHBubble<d>::Opposite_Side_Nb_Num(int i, int ref)
{
	const auto& nbs = surface.Tangential_Neighbor_Of(i);
	VectorT lr_ri = surface.Project_To_TPlane(particles.X(ref) - particles.X(i), particles.E(i));
	int diff_nbs = 0;
	for (int k = 0; k < nbs.size(); k++) {//loop through
		int j = nbs[k];
		VectorT lr_ji = surface.Project_To_TPlane(particles.X(j) - particles.X(i), particles.E(i));
		if (lr_ji.dot(lr_ri) < 0) { // if there is one neighbor that is to the opposite of nearest boundary
			diff_nbs++; //then not out of bound, stop checking
		}
	}
	return diff_nbs;
}

template<int d>
real FluidSPHBubble<d>::Compute_Enclosed_Volume() {
	VectorD center;
	center = VectorD::Zero();
	for (int i = 0; i < particles.Size(); i++) {
		center += particles.X(i);
	}
	center *= (real)1 / particles.Size();
	real vol = 0;
	for (int i = 0; i < particles.Size(); i++) {
		real height = fabs(particles.Normal(i).dot(particles.X(i) - center));
		vol += real(1.0 / d) * particles.SA(i) * height;
		//vol += (real)(1. / d) * particles.SA(i) * (particles.X(i) - center).norm();
		//std::cout << (real)(1) / 3 * particles.SA(i) * (particles.X(i) - center).norm()
	}
	return vol;
}

template<int d>
void FluidSPHBubble<d>::Seed_Droplet(void)
{
	if (num_droplets >= max_droplets) return;
	//if (current_frame == last_droplet_frame) return;
	//if (Rand_Number() > seed_droplet_rate) return;
	//if (current_frame > 1) return;
	real tmp = Rand_Number() * pi / 2. + ((3. * pi) / 4.);
	VectorD drip_center;
	//drip_center << simulation_scale * -1, 0, 0;
	real drip_r = (2.5*length_scale) * (0.8 + 0.0 * Rand_Number());
	drip_center << (simulation_scale-0.5*drip_r) * cos(tmp), 0, (simulation_scale-0.5*drip_r)* sin(tmp);
#pragma omp parallel for
	for (int i = 0; i < particles.Size(); i++) {
		real dist = (particles.X(i) - drip_center).norm();
		if (dist < drip_r) {
			particles.M(i) = droplet_m;
			//particles.V(i) = VectorD::Unit(0);
			//particles.RH(i) = thickness + droplet_rh * std::min(cos(dist / drip_r * pi / 2.0), 1.0);
			particles.RH(i) = droplet_rh;
			particles.Phase(i) = 1;
		};
	}
	last_droplet_frame = current_frame;
	num_droplets++;
	std::cout << "current_droplets: " << num_droplets << std::endl;
	std::cout << "max_droplets: " << max_droplets << std::endl;
}

//template<int d>
//void FluidSPHBubble<d>::Seed_Droplet(void)
//{
//#pragma omp critical
//	{
//		for (int i = 0; i < droplet_nozzles.size();i++) {
//			if (drops[i] >= max_droplets) continue;
//			int k = particles.Add_Element();
//			particles.Copy_Element_From(k, particles, droplet_nozzles[i]);
//			particles.V(k) = VectorD::Zero();
//			particles.X(k) -= 0.6 * length_scale * particles.X(k).normalized();
//			particles.B(k) = 0;
//			particles.Phase(k) = 1;
//			particles.RH(k) = droplet_rh;
//			particles.RH_V(k) = 0;
//			particles.M(k) = droplet_m;
//			drops[i] += 1;
//		}
//	}
//}

template<int d>
void FluidSPHBubble<d>::Update_Catenoid_Boundary(real dt) {
	real curr_catenoid_speed = catenoid_speed;
	if (current_frame > catenoid_stop_frame) { curr_catenoid_speed = catenoid_speed * pow(catenoid_stop_rate, current_frame - catenoid_stop_frame); }
	VectorD trans_vec = VectorD::Unit(1) * dt * curr_catenoid_speed;
#pragma omp parallel for
	for (int i = 0; i < particles.Size(); i++) {
		if (particles.B(i) == 1) {
			//std::cout << particles.BG(i) << std::endl;
			if (particles.X(i)[1] > 0) {
				particles.Apply_Translation(i, trans_vec);
			}
			else if (particles.X(i)[1] < 0) {
				particles.Apply_Translation(i, -trans_vec);
			}
		}
	}
	std::shared_ptr<HalfBowl<d>> lower = std::dynamic_pointer_cast<HalfBowl<d>> (analytical_boundary.obstacles[0]);
	std::shared_ptr<HalfBowl<d>> upper = std::dynamic_pointer_cast<HalfBowl<d>> (analytical_boundary.obstacles[1]);
	lower->center += -VectorD::Unit(1) * dt * curr_catenoid_speed;
	upper->center += VectorD::Unit(1) * dt * curr_catenoid_speed;
}

template<int d>
int FluidSPHBubble<d>::Merge_3d_Particles(void) {
	Array<std::string> att_list{ "x","m","conc","rh","e","v"};
	static Array<int> to_delete; to_delete.resize(fluid_3d->particles.Size()); AuxFunc::Fill(to_delete, 0);
	for (int i = 0; i < fluid_3d->particles.Size(); i++) {
		real near_radius = 2 * fluid_3d->length_scale;
		Array<int> surface_nbs = surface.nbs_searcher->Find_Neighbors(fluid_3d->particles.X(i), near_radius);
		if (surface_nbs.size() < 1) continue;
		VectorD diff = particles.X(surface_nbs[0]) - fluid_3d->particles.X(i);
		if (abs(diff.dot(particles.Normal(surface_nbs[0]))) < 0.3 * length_scale) {
			int k = particles.Add_Element();
			particles.Copy_Element_From(k, fluid_3d->particles, i, att_list);
			particles.V(k) *= 0.1;
			particles.F(k) = VectorD::Zero();
			particles.E(k) = particles.E(surface_nbs[0]);
			particles.Vol(k) = particles.Vol(surface_nbs[0]);
			to_delete[i] = true;
		}
	}
	int before_size = fluid_3d->particles.Size();
	int after_size = fluid_3d->particles.Delete_Elements(to_delete);
	return before_size - after_size;
}

template<int d>
bool FluidSPHBubble<d>::Is_Solitary(int i, const std::string& mode)const
{
	const NormalPressureParams& params = n_pressure_params;
	if (mode == "catenoid") {
		auto normal_func = [&](const int idx)->VectorD {return particles.Normal(idx); };
		VectorD normal_laplacian = Surface_Laplacian<VectorD>(i, normal_func, params.opt_mode);
		Array<int> solitary_nbs;
		real nearby = 5 * length_scale;
		fluid_3d->nbs_searcher->Find_Neighbors(particles.X(i), surface.v_r, solitary_nbs);
		int nbs_3d_num = solitary_nbs.size(), nbs_2d_num = surface.Tangential_Neighbor_Of(i).size();
		//return particles.X(i).norm() < 0.2 * simulation_scale || normal_laplacian.norm() > 7000 || surface.Tangential_Neighbor_Of(i).size() < 2;
		return normal_laplacian.norm() > 10000 || surface.Tangential_Neighbor_Of(i).size() < 2||nbs_3d_num>3*nbs_2d_num;
	}
	else if (mode == "bursting_bubble") {
		assert(fluid_3d != nullptr);
		real nearby = 2 * length_scale;
		Array<int> solitary_nbs;
		fluid_3d->nbs_searcher->Find_Neighbors(particles.X(i), nearby, solitary_nbs);
		auto normal_func = [&](const int idx)->VectorD {return particles.Normal(idx); };
		VectorD normal_laplacian = Surface_Laplacian<VectorD>(i, normal_func, params.opt_mode);
		return surface.Tangential_Neighbor_Of(i).size() < 2 || normal_laplacian.norm() > 10000;
		//if (i == 959) std::cout << particles.H(i) << std::endl;
		//return normal_laplacian.norm() > 10000 || surface.Tangential_Neighbor_Of(i).size() < 2 || particles.H(i) < 0.6 * avg_height;
	}
	else if (mode == "dripping") {
		//if (current_frame <= 100) return false;
		real near_radius = surface.t_r;
		Array<int> surface_nbs = fluid_3d->nbs_searcher->Find_Neighbors(particles.X(i), near_radius);
		int nbs_3d_num = surface_nbs.size(), nbs_2d_num = surface.Tangential_Neighbor_Of(i).size();
		auto normal_func = [&](const int idx)->VectorD {return particles.Normal(idx); };
		VectorD normal_laplacian = Surface_Laplacian<VectorD>(i, normal_func, params.opt_mode);
		return normal_laplacian.norm() > 8000 || particles.H(i) > 2.5e-6 || nbs_3d_num >= 2 * nbs_2d_num;
	}
	else if (mode=="generic"){
		return false;
	}
	else {
		std::cerr << "FluidSPHBubble<d>::Is_Solitary error: unknown mode\n"; exit(0);
	}
	return false;
}

template<int d>
void FluidSPHBubble<d>::Mark_Solitary_Particles(Array<int>& is_solitary)
{	
	if (d == 2) return;
	assert(is_solitary.size() == particles.Size());
	//last_mark_solitary_frame = current_frame;

#pragma omp parallel for
	for (int i = 0; i < particles.Size(); i++) {
		if (particles.Is_Boundary(i))continue;
		if(Is_Solitary(i,exp_mode)){
			is_solitary[i] = true;
			if (fluid_3d != nullptr) {
#pragma omp critical
				{
					if (exp_mode == "catenoid"){ fluid_3d->Add_Particle(particles.X(i), 0.2 * particles.V(i), particles.M(i), particles.RH(i)); }
					//else { fluid_3d->Add_Particle(particles.X(i), particles.V(i), particles.M(i), particles.RH(i)); }
				}
			}
		}
	}

	//Perform a BFS for exp_mode=="dripping"
	/*if (exp_mode == "dripping") {
		real threshold = 2 * length_scale;
		int pn = particles.Size();
		Array<int> is_bubble(pn);
		std::queue<int> Q;
		for (int i = 0; i < pn; i++) {
			if (particles.X(i)[1] > simulation_scale * 2.5 && !is_solitary[i]) {
				is_bubble[i] = true;
				Q.push(i);
			}
			else is_bubble[i] = false;
		}
		while (!Q.empty()) {
			int x = Q.front(); Q.pop();
			const auto& nbs = surface.Tangential_Neighbor_Of(x);
			for (int i = 0; i < nbs.size(); i++) {
				int y = nbs[i];
				if ((particles.X(x) - particles.X(y)).norm() < threshold && !is_bubble[y] && !is_solitary[y]) {
					is_bubble[y] = true; Q.push(y);
				}
			}
		}
		for (int i = 0; i < pn; i++) is_solitary[i] = is_solitary[i] || (!is_bubble[i]);
	}*/

	int temp_solitary_num = 0;
#pragma omp parallel for reduction (+: temp_solitary_num)
	for (int i = 0; i < particles.Size(); i++) {
		temp_solitary_num += (int)is_solitary[i];
	}
	solitary_num = temp_solitary_num;
	if (temp_solitary_num > 0)std::cout << "Mark " << solitary_num << " solitary numbers\n";
}

template<int d>
void FluidSPHBubble<d>::Mark_Idle_Particles(Array<int>& is_idle)
{
	if (d == 2) return;
	//assert(is_solitary.size() == particles.Size());

#pragma omp parallel for
	for (int i = 0; i < particles.Size(); i++) {
		if (particles.Is_Boundary(i))continue;
		if (analytical_boundary.Available() && phis[i] < 0.51 * length_scale) {
			is_idle[i] = true;
		}
	}

	int temp_solitary_num = 0;
#pragma omp parallel for reduction (+: temp_solitary_num)
	for (int i = 0; i < particles.Size(); i++) {
		temp_solitary_num += (int)is_idle[i];
	}
	solitary_num = temp_solitary_num;
	if (temp_solitary_num > 0)std::cout << "Mark " << solitary_num << " idle numbers\n";
}

template<int d>
void FluidSPHBubble<d>::Apply_3D_Surface_Tension(void)
{
	if (d == 2) return;
	if (fluid_3d == nullptr) return;
	fluid_3d->nbs_searcher->Update_Points(fluid_3d->particles.XRef());
#pragma omp parallel for
	for (int i = 0; i < particles.Size(); i++) {
		//surface.nbs_searcher->Update_Points(particles.XRef());
		real nearby = 2 * length_scale;
		Array<int> solitary_nbs;
		//Array<int> surface_nbs;
		//Array<VectorD> tmp_array = fluid_3d->particles.XRef();
		//tmp_array.push_back(particles.X(i));
	    fluid_3d->nbs_searcher->Find_Neighbors(particles.X(i), nearby, solitary_nbs);
		//surface.nbs_searcher->Find_Neighbors(particles.X(i), nearby, surface_nbs);
		//if (solitary_nbs.size() > 0.2 * surface.Tangential_Neighbor_Of(i).size() || surface.Tangential_Neighbor_Of(i).size() < 2) {
		//if ((solitary_nbs.size() > 0.3 * surface_nbs.size() && Rand_Number() < 0.2) || surface.Tangential_Neighbor_Of(i).size() < 2){
		//if (Rand_Number() < ((real)solitary_nbs.size()/ (real)surface_nbs.size()) || surface.Tangential_Neighbor_Of(i).size() < 2) {
		for (int k = 0; k < solitary_nbs.size(); k++) {
			int j = solitary_nbs[k];
			VectorD r_ij; 
			r_ij = particles.X(i) - fluid_3d->particles.X(j);
			VectorD f_cohesion = VectorD::Zero(); 
			if (acos(particles.Normal(i).dot(r_ij)/(r_ij.norm()))> pi/4.) f_cohesion = -1 * -particles.M(i) * fluid_3d->particles.M(j) * kernel_v->W_Cubic(r_ij.norm()) * r_ij.normalized();
			//f_cohesion -= f_cohesion.dot(particles.Normal(i)) * particles.Normal(i);
			particles.F(i) += 3000. * f_cohesion;
		}
	}
}


template<int d>
void FluidSPHBubble<d>::Apply_Multiphase_Cohesion(void)
{
	for (int i = 0; i < particles.Size(); i++) {
		if (particles.Is_Boundary(i)) continue;
		if (particles.Phase(i) < 1) continue;
		const auto& nbs = surface.Tangential_Neighbor_Of(i);
		for (int k = 0; k < nbs.size(); k++) {
			int j = nbs[k];
			if (particles.Phase(j) != particles.Phase(i)) continue;
			VectorD r_ij = particles.X(i) - particles.X(j);
			particles.F(i) += -.3 * particles.M(i) * r_ij.normalized() * Cohesion_Kernel(r_ij.norm(), surface.v_r) * 1e-5;
		}
	}
}

template<int d>
void FluidSPHBubble<d>::Mark_Fast_Particles(Array<int>& is_fast, real vel_threshold)
{
	assert(is_fast.size() == particles.Size());
#pragma omp parallel for
	for (int i = 0; i < particles.Size(); i++) {
		if (particles.V(i).norm() > vel_threshold) {
			is_fast[i] = 1;
		}
	}
}

template<int d>
void FluidSPHBubble<d>::Transform_3D_Particles(const Array<int>& to_transform)
{
	assert(to_transform.size() == particles.Size());
	assert(fluid_3d != nullptr);
	for (int i = 0; i < particles.Size(); i++) {
		if (to_transform[i]) {
			//fluid_3d->Add_Particle(particles.X(i), particles.V(i), particles.M(i), particles.RH(i));
			if (exp_mode == "catenoid") { fluid_3d->Add_Particle(particles.X(i), 0.2 * particles.V(i), particles.M(i), particles.RH(i)); }
			else { fluid_3d->Add_Particle(particles.X(i), particles.V(i), particles.M(i), particles.RH(i)); }
		}
	}
}

template<int d>
real FluidSPHBubble<d>::Surface_Weighted_Sum(int i, std::function<real(const int)> f, const OperatorParams& params, MirrorCompensator<d>& mcp)
{
	const auto& nbs = surface.Tangential_Neighbor_Of(i);
	real sum = 0, coeff = 1;
	for (int k = 0; k < nbs.size(); k++) {
		coeff = 1;
		int j = nbs[k];
		if (!In_Calc_Field(j, params) && j != i) continue;
		real val = 1; if (f != nullptr) val = f(j);
		VectorT lr_ji = -Relative_Vector_Surface(i, j);
		coeff += mcp.Coeff_Offset(-lr_ji, 0);

		sum += coeff * val * Surface_Kernel_Weight(lr_ji.norm(), params.kernel);
	}
	return sum;
}

template<int d>
Vector<real, d - 1> FluidSPHBubble<d>::Surface_Gradient_Symmetric(int i, std::function<real(const int)> f, const OperatorParams& params, bool force_symmetric)
{
	//std::cout << "wtf:: " << params.calculate_field << std::endl;
	const auto& nbs = surface.Tangential_Neighbor_Of(i);
	VectorT grad_f = VectorT::Zero();
	const real& V_i = particles.Vol(i);
	const real& H_i = particles.H(i);
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		if (!In_Calc_Field(j, params) && j != i) continue;
		const real& V_j = particles.Vol(j);
		real H_j = particles.H(j);
		VectorT lr_ji = -Relative_Vector_Surface(i, j);
		VectorT grad_W = Surface_Kernel_Grad(lr_ji, params.kernel);
		if (!force_symmetric) { grad_f += V_j * (f(i) / (H_i * H_i) + f(j) / (H_j * H_j)) * grad_W; }
		else { grad_f += std::max(V_j, V_i) * (f(i) / (H_i * H_i) + f(j) / (H_j * H_j)) * grad_W; }
	}
	return grad_f * H_i;
}

template<int d>
Vector<real, d - 1> FluidSPHBubble<d>::Surface_Gradient_Difference(int i, std::function<real(const int)> f, const OperatorParams& params)
{
	//if (i == 10)std::cout << "calculate symmetric fi: " << f(i) << std::endl;
	const auto& nbs = surface.Tangential_Neighbor_Of(i);
	VectorT grad_f = VectorT::Zero();
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		if (!In_Calc_Field(j, params) && j != i) continue;
		real S_j = particles.SA(j);
		if (particles.Is_Boundary(j)) S_j = particles.SA(i);
		VectorT lr_ji = -Relative_Vector_Surface(i, j);
		VectorT grad_W = Surface_Kernel_Grad(lr_ji, params.kernel);
		grad_f += S_j * (f(j) - f(i)) * grad_W;
		//if (i == 0) { std::cout << "grad diff of " << i << " j= " << j << " term= " << S_j << " *( " << f(j)<<" - "<<f(i) << " )* " << grad_W << "\n"; }
	}
	return grad_f;
}

template<int d>
real FluidSPHBubble<d>::Surface_Divergence(int i, std::function<VectorD(const int)> f, const OperatorParams& params)
{
	const auto& nbs = surface.Tangential_Neighbor_Of(i);
	real div = 0;
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		if (!In_Calc_Field(j, params) && j != i) continue;
		real S_j = particles.SA(j);
		if (particles.Is_Boundary(j))S_j = particles.SA(i);

		VectorT lr_ji = -Relative_Vector_Surface(i, j);
		VectorT grad_W = Surface_Kernel_Grad(lr_ji, params.kernel);
		VectorT fr_ij = World_Vector_To_Surface(i, f(j) - f(i));
		div += S_j * fr_ij.dot(grad_W);

		//if (i == 2) {
		//	std::cout << "====" << std::endl;
		//	std::cout << "j is: " << j << std::endl;
		//	if (particles.Is_Boundary(j)) {
		//		std::cout << "fucked up " << std::endl;
		//		std::cout << "wtf? " << f(j) << std::endl;
		//	}
		//	std::cout << "contribution: " <
		//	std::cout << "====" << std::endl;
		//}
	}
	return div;
}

template<int d>
Vector<real, d> FluidSPHBubble<d>::Gravity_Adhesion(int i, const OperatorParams& params) {
	const auto& nbs = surface.Tangential_Neighbor_Of(i);
	VectorD grav_force = VectorD::Zero();
	int num_bnd_nbs = 0;
	for (int k = 0; k < nbs.size(); k++) {
		int j = nbs[k];
		if (i == j) continue;
		if (!In_Calc_Field(j, params)) continue;
		VectorD lr_ij = Relative_Vector_World(i, j);
		//VectorT lr_center; surface.Project_To_TPlane(VectorD::Zero() - particles.X(i), particles.E(i), lr_center); //project origin to the current surface to represent "inner". Assumes origin is always within the boundary
		//if (lr_ji.dot(lr_center)<0)
		VectorD lr_ij_scaled = lr_ij / simulation_scale;
		if (lr_ij_scaled.norm() > 0.1) continue; //condition to determine if is REALLY near boundary
		if (lr_ij_scaled.norm() < 0.001) continue; //if TOO close then stop apply this gravity
		grav_force += 50 * simulation_scale * (lr_ij / lr_ij.norm()) * (1 / (lr_ij_scaled.norm() * lr_ij_scaled.norm()));
		num_bnd_nbs++;
	}

	//if (num_bnd_nbs < 1) { return grav_force; }
	//else { return 1 / num_bnd_nbs * grav_force; }
	return  boundary_params.grav_boundary * (particles.M(i) * particles.M(i)) * grav_force;
}

template<int d>
Vector<real, d> FluidSPHBubble<d>::Boundary_Force(int i)
{
	const BoundaryParams& params = boundary_params;
	if (params.boundary_force_mode == "binary") {
		auto nominal_vel_i = [&](const int idx)->VectorD {return Nominal_Vector(particles.X(i), particles.Normal(i), particles.X(idx), particles.Normal(idx), particles.V(idx));  };
		OperatorParams vis_bnd_opt = viscosity_params, elas_bnd_opt("only_boundary", KernelType::SPIKY);
		vis_bnd_opt.calculate_field = "only_boundary";
		//VectorD vis_world = particles.Vol(i)* viscosity_water* Surface_Laplacian(i, Index_Function(particles.VRef()), vis_bnd_params, false)* vis_boundary;
		VectorD vis_world = particles.Vol(i) * viscosity_water * Surface_Laplacian<VectorD>(i, nominal_vel_i, vis_bnd_opt, false) * params.vis_boundary;
		//if (Gravity_Adhesion(i, elas_bnd_params).norm() > 0) std::cout << "what is gravity adhesion: " << Gravity_Adhesion(i, elas_bnd_params) << std::endl;
		return vis_world + Gravity_Adhesion(i, elas_bnd_opt);
	}
	if (params.boundary_force_mode == "vis") {
		auto nominal_vel_i = [&](const int idx)->VectorD {return Nominal_Vector(particles.X(i), particles.Normal(i), particles.X(idx), particles.Normal(idx), particles.V(idx));  };
		OperatorParams vis_bnd_opt = viscosity_params, elas_bnd_opt("only_boundary", KernelType::SPIKY);
		vis_bnd_opt.calculate_field = "only_boundary";
		VectorD vis_world = particles.Vol(i) * viscosity_water * Surface_Laplacian<VectorD>(i, nominal_vel_i, vis_bnd_opt, false) * params.vis_boundary;
		return vis_world;
	}
	else if (params.boundary_force_mode == "none") {
		return VectorD::Zero();
	}
	else if (params.boundary_force_mode == "grav") {
		OperatorParams elas_bnd_params("only_boundary", KernelType::SPIKY);
		if (particles.H(i) < avg_height) return Gravity_Adhesion(i, elas_bnd_params);
	}
	else {
		std::cerr << "Boundary_Force error: unknown boundary force mode\n";
		return VectorD::Zero();
	}
}

template<int d>
void FluidSPHBubble<d>::Diffuse_Scalar(Array<real>& f, const real& dt, const real& coeff, const OperatorParams& params)
{
	Array<real> diff_v;
	diff_v.resize(particles.Size());
#pragma omp parallel for
	for (int i = 0; i < particles.Size(); i++) {
		if (particles.Is_Boundary(i)) continue;
		diff_v[i] = coeff * Surface_Laplacian(i, Index_Function(f), params);
	}
#pragma omp parallel for
	for (int i = 0; i < particles.Size(); i++) {
		f[i] += diff_v[i] * dt;
	}
}

template<int d>
void FluidSPHBubble<d>::Update_Fluid_Pressure(void)
{
	const TangentialPressureParams& params = t_pressure_params;
#pragma omp parallel for //reduction (+:total_div_pressure_tmp,total_height_pressure_tmp,total_height_laplacian_pressure_tmp)
	for (int i = 0; i < particles.Size(); i++) {
		particles.KH(i) = Surface_Laplacian<real>(i, Index_Function(particles.HRef()), height_laplacian_params, false);
		real div = particles.Div(i);

		real divergence_pressure = -params.divergence_pressure_coeff * div;
		real height_pressure = params.height_pressure_coeff * particles.H(i) / thickness;
		real height_laplacian_pressure = params.laplacian_pressure_coeff * -2 * particles.Gamma(i) * particles.KH(i);
		real boundary_pressure = 0.;
		if (particles.Is_Boundary(i)) {
			boundary_pressure = 0.;
		}
		else {
			boundary_pressure = params.boundary_pressure_coeff * (avg_height - particles.H(i));
		}

		if (analytical_boundary.Available() && phis[i] < 2 * length_scale) {
			height_laplacian_pressure *= 0.;
			divergence_pressure *= 0.;
		}

		divergence_pressure *= params.tangential_pressure_coeff;
		height_laplacian_pressure *= params.tangential_pressure_coeff;
		height_pressure *= params.tangential_pressure_coeff;
		boundary_pressure *= params.tangential_pressure_coeff;
		particles.P0(i) = height_pressure;
		particles.P1(i) = divergence_pressure;
		particles.P2(i) = height_laplacian_pressure;
		particles.PB(i) = boundary_pressure;
		particles.P(i) = divergence_pressure;
		particles.P(i) += height_pressure;
		particles.P(i) += height_laplacian_pressure;
		particles.P(i) += boundary_pressure;
	}
}

template<int d>
void FluidSPHBubble<d>::Update_Gravity_Forces(void)
{
#pragma omp parallel for
	for (int i = 0; i < particles.Size(); i++) {
		particles.F(i) += g * particles.M(i) * gravity_coeff;
	}
}

template<int d>
void FluidSPHBubble<d>::Update_Friction_Forces(void)
{
	if (air_velocity_func == nullptr) return;//no friction
#pragma omp parallel for
	for (int i = 0; i < particles.Size(); i++) {
		VectorD air_vel = air_velocity_func(particles.X(i));
		VectorD vel_diff = air_vel - particles.V(i);
		particles.F(i) += vel_diff * particles.M(i) * friction_coeff;
	}
}

template<int d>
void FluidSPHBubble<d>::Update_External_Forces(void)
{
	if (external_force_func != nullptr) {
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			particles.F(i) += external_force_func(i);
		}
	}
}

template<int d>
void FluidSPHBubble<d>::Update_Boundary_Forces(void)
{
#pragma omp parallel for
	for (int i = 0; i < particles.Size(); i++) {
		particles.F(i) += Boundary_Force(i) * boundary_params.boundary_force_coeff;
	}
}

template<int d>
void FluidSPHBubble<d>::Update_Viscosity_Forces(void)
{
	//real total_vis_force_tmp = 0;
	VectorD curr_vis_force;
#pragma omp parallel for //reduction (+:total_vis_force_tmp)
	for (int i = 0; i < particles.Size(); i++) {
		VectorD vis_force = particles.Vol(i) * viscosity_water * Surface_Laplacian(i, Index_Function(particles.VRef()), viscosity_params);
		VectorD norm = particles.Normal(i);
		VectorD vis_normal = norm * vis_force.dot(norm) / norm.squaredNorm();//normal component
		VectorD vis_tangential = vis_force - vis_normal;//tangential component
		curr_vis_force = vis_normal * normal_viscosity_coeff + vis_tangential * viscosity_coeff;
		particles.F(i) += curr_vis_force;
		if (!diagnosis) continue;
		if (particles.Is_Boundary(i)) {
			total_vis_force[i] = 0.;
		}
		else {
			total_vis_force[i] = curr_vis_force.norm();
		}
	}
}

template<int d>
void FluidSPHBubble<d>::Update_Vorticity_Confinement_Forces(void)
{
	if constexpr (d == 3) {
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			VectorD norm = particles.Normal(i);
			VectorD force = VectorD::Zero();
			const auto& nbs = surface.Tangential_Neighbor_Of(i);
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				if (i == j) continue;
				VectorT lr_ij = Relative_Vector_Surface(i, j);
				VectorD r_ij = particles.X(j) - particles.X(i);
				real w = Surface_Kernel_Weight(lr_ij.norm(), vorticity_params.weight_kernel);
				VectorD omega = w * particles.Vrt(j) * norm;
				VectorD n_p = AuxFunc::Eliminate_Unit_Component(r_ij, norm).normalized();
				force += vorticity_params.force_coeff * (n_p.cross(omega));
			}
			particles.F(i) += force * particles.M(i);
		}
	}
}

template<int d>
void FluidSPHBubble<d>::Update_Tangential_Pressure_Forces(void)
{
#pragma omp parallel for
	for (int i = 0; i < particles.Size(); i++) {
		VectorT P0F, P1F, P2F;
		VectorT PBF;

		bool is_inside_particle = true;
		if (boundary_params.boundary_mode == "analytical" && analytical_boundary.Available()) {
			if (phis[i] <= np_on_h * length_scale) {
				is_inside_particle = false;
			}
		}
		if (is_inside_particle) {
			P0F = particles.Vol(i) * -1 * Surface_Gradient_Symmetric(i, Index_Function(particles.P0Ref()), grad_force_params);
			P1F = particles.Vol(i) * -1 * Surface_Gradient_Symmetric(i, Index_Function(particles.P1Ref()), grad_force_params);
			P2F = particles.Vol(i) * -1 * Surface_Gradient_Symmetric(i, Index_Function(particles.P2Ref()), grad_force_params);
		}
		else {
			P0F = particles.Vol(i) * -1 * Surface_Gradient_Difference(i, Index_Function(particles.P0Ref()), grad_force_params);
			P1F = particles.Vol(i) * -1 * Surface_Gradient_Difference(i, Index_Function(particles.P1Ref()), grad_force_params);
			P2F = particles.Vol(i) * -1 * Surface_Gradient_Difference(i, Index_Function(particles.P2Ref()), grad_force_params);
		}
		//P0F = particles.Vol(i) * -1 * Surface_Gradient_Difference(i, Index_Function(particles.P0Ref()), grad_force_params);
		//P1F = particles.Vol(i) * -1 * Surface_Gradient_Symmetric(i, Index_Function(particles.P1Ref()), divergence_params);
		//P2F = particles.Vol(i) * -1 * Surface_Gradient_Symmetric(i, Index_Function(particles.P2Ref()), grad_force_params);
		PBF = particles.Vol(i) * -1 * Surface_Gradient_Difference(i, Index_Function(particles.PBRef()), OperatorParams("only_boundary", KernelType::SPIKY));


		VectorT tang_pressure_force = P0F + P1F + P2F + PBF; //p0

		VectorD curr_pressure_force = Surface_Vector_To_World(i, tang_pressure_force);

		//int bnd_nb = Nearest_Boundary_Neighbor(i);
		//if (bnd_nb != -1) {
		//	VectorT bnd_pressure_force = particles.Vol(i) * -1 * Surface_Gradient_Difference(i, Index_Function(particles.PBRef()), boundary_params); //p4
		//	//VectorD diff = particles.X(bnd_nb) - particles.X(i);
		//	PBF = Surface_Vector_To_World(i, bnd_pressure_force);
		//	std::cout << "curr_boundary_pressure_force" << bnd_pressure_force.norm() << std::endl;
		//	curr_pressure_force += PBF;
		//}


		particles.F(i) += curr_pressure_force;
		if (!diagnosis) continue;
		if (particles.Is_Boundary(i)) {
			total_pressure_force[i] = 0.;
			total_div_pressure_force[i] = 0.;
			total_height_pressure_force[i] = 0.;
			total_height_laplacian_pressure_force[i] = 0;
			total_bnd_pressure_force[i] = 0;
		}
		else {
			total_pressure_force[i] = curr_pressure_force.norm();
			total_div_pressure_force[i] = P1F.norm();
			total_height_pressure_force[i] = P0F.norm();
			total_height_laplacian_pressure_force[i] = P2F.norm();
			total_bnd_pressure_force[i] = PBF.norm();
		}
	}
}

template<int d>
void FluidSPHBubble<d>::Update_Marangoni_Forces(void)
{
#pragma omp parallel for
	for (int i = 0; i < particles.Size(); i++) {
		VectorT tang_mar_force = particles.SA(i) * Surface_Gradient_Difference(i, Index_Function(particles.GammaRef()), grad_force_params) * marangoni_coeff;
		particles.F(i) += Surface_Vector_To_World(i, tang_mar_force);
	}
}

template<int d>
void FluidSPHBubble<d>::Update_Normal_Pressure_Forces(real dt)
{
	const NormalPressureParams& params = n_pressure_params;
//	std::cout << "current capillary coeff " << params.capillary_coeff << std::endl;
	real enclosed_volume = 0, enclosed_pressure = 0;
	if (params.closed) {
		enclosed_volume = Compute_Enclosed_Volume();
		if (enclosed_volume > 1e-8) { 
			if (params.pressure_decay > 0.) {
				enclosed_pressure = params.decayed_pressure / enclosed_volume;
				std::cout << "pressure_decay is added, current_enclosed pressure: " << enclosed_pressure << std::endl;
			}
			else{
				enclosed_pressure = params.pressure_constant / enclosed_volume;
			}
		}
	}
#pragma omp parallel for
	for (int i = 0; i < particles.Size(); i++) {
		auto z_func_i = [&](const int idx)->real {return particles.X(idx).dot(particles.Normal(i)); };
		real pressure_difference = 0;
		particles.KH(i) = Surface_Laplacian<real>(i, z_func_i, params.opt_mode, true);
		pressure_difference += 4 * particles.Gamma(i) * particles.KH(i);//capillary
		//std::cout << particles.KH(i) << std::endl;
		if (params.air_pressure_mode == "pv") {
			//along normal direction, if normal points outward then it's inner-outside
			real sgn = 1;
			if (particles.X(i).dot(particles.Normal(i)) < 0) sgn = -1;
			else sgn = 1;
			pressure_difference += sgn * (enclosed_pressure - params.atmo_pressure) * params.pv_air_force_coeff;
		}
		//real c02 = 4 * particles.Gamma(i) * particles.SA(i) * params.capillary_coeff / particles.M(i);
		//real c02 = 4 * particles.Gamma(i) * params.capillary_coeff / (1e3 * thickness);
		//real predicted_multiplier = 4 * n_pressure_params.capillary_coeff * gamma_soap / (1e3 * thickness);
		//std::cout << "gamma soap: " << gamma_soap << " gamma particle: " << particles.Gamma(i) << " capillary param: " << n_pressure_params.capillary_coeff << " capillary now: " << params.capillary_coeff << "\n";
		//std::cout << "M: " << particles.M(i) << " Vol: " << particles.Vol(i) << " SA: " << particles.SA(i) << " H: " << particles.H(i) << "\n";
		//std::cout << "c02: " << c02 << " predicted: " << predicted_multiplier << "\n";
		//std::cout << "phase velocity predicted: " << std::setiosflags(std::ios::fixed) << std::setprecision(8) << sqrt(c02) << "\n";
		particles.F(i) += pressure_difference * particles.SA(i) * params.capillary_coeff * particles.Normal(i);
		if (params.air_pressure_mode == "ib") {
			assert(false && "air pressure mode ib is not implemented for open source version.");
			// Interpolation<d> intp(air_solver.mac_grid);
			// VectorD ib_vel = intp.Interpolate_Face_Vectors(air_solver.velocity, particles.X(i));
			// VectorD ib_friction_force = (ib_vel - particles.V(i)) * particles.M(i) * params.ib_force_coeff;
			// particles.F(i) += ib_friction_force;
			// //particles.F(i) += ib_friction_force.dot(particles.Normal(i)) * particles.Normal(i);
		}
		//particles.F(i) += -particles.M(i) * 0.5 * particles.Normal(i);// this is just to see what happens, please delete
	}
	//IB_Pressure_Gradient(idx, dt) * params.ib_pressure_coeff * particles.H(idx);
}

template<int d>
void FluidSPHBubble<d>::Update_Divergence(void)
{
#pragma omp parallel for
	for (int i = 0; i < particles.Size(); i++) {
		auto nominal_vel_i = [&](const int idx)->VectorD {return Nominal_Vector(particles.X(i), particles.Normal(i), particles.X(idx), particles.Normal(idx), particles.V(idx));  };
		//real div = Surface_Divergence(i, Index_Function(particles.VRef()), divergence_params);
		real div = Surface_Divergence(i, nominal_vel_i, divergence_params);
		particles.Div(i) = div;
	}
}

template<int d>
void FluidSPHBubble<d>::Fix_Render_Height(void) {
	if (rh_params.fix_RH) {
		real total_render_vol = 0;
#pragma omp parallel for reduction(+:total_render_vol)
		for (int i = 0; i < particles.Size(); i++) {
			//std::cout << particles.SA(i) << std::endl;
			total_render_vol += particles.SA(i) * particles.RH(i);
		}
		real coeff = total_volume / total_render_vol;
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			particles.RH(i) *= coeff;
			//std::cout << "RH: " << particles.RH(i) << std::endl;
		}
	}
}

template<int d>
void FluidSPHBubble<d>::Update_Render_Height(real dt)
{
	if (rh_params.mode == "divergence") {
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			real rh = particles.RH(i);
			real div = particles.Div(i);
			real DhDt = -rh * div;
			particles.RH(i) += DhDt * dt;
		}
	}
	else if (rh_params.mode == "laplacian") {
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			real kh = Surface_Laplacian<real>(i, Index_Function(particles.RHRef()), height_laplacian_params, false);
			real gamma = Surface_Tension_Coefficient(particles.Conc(i));
			real rho = particles.M(i) / particles.Vol(i);
			particles.RH_V(i) += gamma / rho * kh * dt;
		}
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			particles.RH(i) += particles.RH_V(i) * dt;
		}
	}
	else {
		std::cerr << "FluidSPHBubble<d>::Update_Render_Height error: unrecognized render height mode\n"; exit(0);
	}
	//if blend with SPH h, blend
	if (rh_params.blend_h) {
		real alpha = rh_params.blend_constant;
#pragma omp parallel for
		for (int i = 0; i < particles.Size(); i++) {
			particles.RH(i) = (1. - alpha) * particles.RH(i) + alpha * particles.H(i);
		}
	}
	//now fix with total SA
	Fix_Render_Height();
}

template<int d>
void FluidSPHBubble<d>::Recompute_Velocities(real dt)
{
#pragma omp parallel for
	for (int i = 0; i < particles.Size(); i++) {
		if (particles.Is_Boundary(i)) continue;
		VectorD increment;
		increment = particles.F(i) / particles.M(i) * dt;
		if (clip_velocity) {
			real vel_norm = increment.norm();
			if (vel_norm > vel_threshold) {
				//std::cout << "clipped velocity! " << std::endl;
				//std::cout << "original: "  << particles.V(i) << std::endl;
				increment = increment * vel_threshold / vel_norm;
				//std::cout << "clipped: " << particles.V(i) << std::endl;
			}
		}
		particles.V(i) += increment;
	}
#pragma omp parallel for
	for (int i = 0; i < particles.Size(); i++) {
		for (int axis = 0; axis < d; axis++) particles.F(i)[axis] = 0;
	}
}

template<int d>
void FluidSPHBubble<d>::Recompute_Geometries(void)
{
	if (analytical_boundary.Available()) {
		//std::cout << "analytical boundary is available!" << std::endl;
		Update_Phis();
	}
	surface.Update();
	Update_Particle_Heights();

	avg_height = 0; real num = 0;
	for (int i = 0; i < particles.Size();i++) {
		if (particles.Is_Boundary(i)) continue;
		avg_height += particles.H(i);
		num += 1.0;
	}
	avg_height /= num;

	real total_SA = 0.;

	for (int i = 0; i < particles.Size(); i++) {
		//total_SA += particles.SA(i);
		if (!particles.Is_Boundary(i))total_SA += particles.SA(i);
	}
	if (diagnosis)std::cout << "Currently, total SA = " << total_SA << std::endl;
	Update_Metric_Tensor();
	//Compute_Color();
}

template<int d>
void FluidSPHBubble<d>::Correct_Velocity_With_Analytical_Boundary(int i, real dt)
{
	if (!analytical_boundary.Available()) {
		std::cerr << "FluidSPHBubble<d>::Correct_Velocity_With_Analytical_Boundary error: analytical boundary not available\n";
	}
	if (particles.Is_Boundary(i)) {
		return;
	}
	real radius = length_scale * 0.5;
	//if (i == 170 && iter >= 231) { std::cout << "analytical boundary particle " << i << " at " << particles.X(i).transpose() << " with phi=" << phi<<",radius="<<radius << "\n"; }
	real phi = phis[i];
	VectorD normal = bnd_normals[i];
	if (phi < radius) {
		VectorD tangential_velocity = particles.V(i) - particles.V(i).dot(normal) * normal;
		if (phi > 0) {//in "normal" fluid region
			particles.V(i) = tangential_velocity;
		}
		else {
			real normal_rate = (-phi + length_scale * 0.5) / dt;
			particles.V(i) = tangential_velocity + normal_rate * normal;
		}
		if constexpr (d == 3) {
			VectorD unit_vel = particles.Normal(i).cross(phi * normal);
			//particles.Vrt(i) = unit_vel.dot(tangential_velocity);
		}
	}
	if (boundary_params.keep_xz_plane && phi > 0 && phi < 1 * length_scale) particles.X(i)[1] = 0.;
}

template<int d>
void FluidSPHBubble<d>::Enforce_Boundary(int i, real dt)
{
	const BoundaryParams& params = boundary_params;
	if (params.boundary_mode == "slippery") {
		int bnd_nb = Nearest_Boundary_Neighbor(i);
		if (bnd_nb != -1) {
			VectorD bnd_vec_3d = Relative_Vector_World(i, bnd_nb);
			if (bnd_vec_3d.norm() < length_scale) { //first make sure that the point is very very close to the nearest boundary
				//find the velocity component that is INTO the boundary (to the same direction of bnd_vec)
				//if V_into is actually "in to", not "away from", cancel out this component
				VectorD V_into = Project_Vector_Along<d>(particles.V(i), bnd_vec_3d);
				if (V_into.dot(bnd_vec_3d) > 0) particles.V(i) -= V_into;
			}
		}
	}
	else if (params.boundary_mode == "partial") {
		int bnd_nb = Nearest_Boundary_Neighbor(i);
		if (bnd_nb != -1) {
			if (Opposite_Side_Nb_Num(i, bnd_nb) < 0.4 * Tangential_Neighbor_Size(i)) {
				VectorT bnd_vec = Relative_Vector_Surface(i, bnd_nb);
				if (bnd_vec.norm() < 0.05) {
					////then we need to eliminate the tangential normal component
					VectorT V_surface = surface.Project_To_TPlane(particles.V(i), particles.E(i));
					real v_dot = V_surface.dot(bnd_vec);
					if (v_dot < 0) { // penalize velocity only when it is moving farther away from the boundary  
						VectorT V_normal = Project_Vector_Along<d - 1>(V_surface, bnd_vec);
						VectorD V_penalty = Surface_Vector_To_World(i, V_normal);
						particles.V(i) -= V_penalty;
					}
				}
			}
		}
	}
	else if (params.boundary_mode == "analytical") {
		Correct_Velocity_With_Analytical_Boundary(i, dt);
	}
	else if (params.boundary_mode == "rotating_boundary") {
		Correct_Velocity_With_Analytical_Boundary(i, dt);
		if (particles.Is_Boundary(i)) {
			if constexpr (d == 3) {
				VectorD tmp_vel;
				tmp_vel << particles.X(i)[2], 0, -particles.X(i)[0];
				particles.V(i) = tmp_vel;
				particles.X(i) += tmp_vel * dt;
			}
		}
	}
	else if (params.boundary_mode == "shear_flow") {
		if (particles.X(i)[1] > 0.5 - length_scale * 0.5) {
			real target_vel = 1;
			//std::cout << "particle " << i << " " << " net force: " << (drag_force + particles.F(i)).transpose() << "\n";
			VectorD dv = VectorD::Unit(0) * (target_vel - particles.V(i)[0]);
			VectorD acc = dv / dt;
			VectorD ghost_force = acc * particles.M(i);
			particles.F(i) = ghost_force;
			particles.V(i) += dv;
		}
		Correct_Velocity_With_Analytical_Boundary(i, dt);
		//real source_rate = 1; VectorD vel_boundary = AuxFunc::V<d>(source_rate);
		/*if (particles.Is_Boundary(i)) {
			if (particles.X(i)[1] > 0) {
				particles.V(i) = vel_boundary;
				particles.X(i) += particles.V(i) * dt;
			}
		}*/
		if (particles.X(i)[0] > 1.0 + length_scale * 0.5) {
			particles.X(i)[0] -= 2.0 + length_scale;
		}
	}
	else if (params.boundary_mode == "none") {
		return;
	}
	else {
		std::cerr << "Enforce_Boundary error: undefined boundary mode\n";
	}
}

template<int d> 
void FluidSPHBubble<d>::Replenish_Boundary_Points_To(SPHBubbleParticles<d>& added_pts, real farthest_dist, real replenish_rate) {
	const int pn = particles.Size();
	static Array<int> to_replenish; to_replenish.resize(pn); AuxFunc::Fill(to_replenish, 0);
	static Array<int> replenish_fluid_nb; replenish_fluid_nb.resize(pn);

	bool check_with_point = !analytical_boundary.Available() || boundary_params.replenish_criteria == "particle";
#pragma omp parallel for
	for (int i = 0; i < pn; i++) {
		if (!particles.Is_Boundary(i)) continue;
		int fluid_nb = Nearest_Fluid_Neighbor(i);
		if (fluid_nb == -1) continue;
		real phi;
		if (check_with_point) {
			phi = (particles.X(i) - particles.X(fluid_nb)).norm();
		}
		else phi = phis[fluid_nb];
		if (phi > farthest_dist) {
			to_replenish[i] = true;
			replenish_fluid_nb[i] = fluid_nb;
		}
	}

	//rand() is thread-unsafe, so please leave this to serial!
	int replenish_num = 0;
	for (int i = 0; i < pn; i++) {
		if (to_replenish[i] && Rand_Number() <= replenish_rate) {
			replenish_num++;
			int new_idx = added_pts.Add_Element();
			int fluid_idx = replenish_fluid_nb[i];
			VectorD mean_pos = 0.5 * (particles.X(i) + particles.X(fluid_idx));
			VectorD rand_offset = 0.3 * Surface_Vector_To_World(fluid_idx, VectorT::Unit(0) * Rand_Number() * length_scale + VectorT::Unit(1) * Rand_Number() * length_scale);
			added_pts.Copy_Element_From(new_idx, particles, fluid_idx);
			added_pts.X(new_idx) = mean_pos + rand_offset;
			added_pts.V(new_idx) = replenish_V_rate * particles.V(fluid_idx);
			if (!boundary_params.inherit_rh) added_pts.RH(new_idx) = avg_height;
		}
	}
	if (replenish_num > 0) { std::cout << " Replenish " << replenish_num << " points\n"; }
}

template<int d>
void FluidSPHBubble<d>::Update_Metric_Tensor(void)
{
	////update metric tensor particles.G --- we need to know more about what this is for
#pragma omp parallel for
	for (int i = 0; i < particles.Size(); i++) {
		if (particles.I(i) == -1)continue;
		const auto& now_nbs = surface.Tangential_Neighbor_Of(i);
		int nb_n = (int)now_nbs.size();
		//std::cout << "for particle: " << i << " number of nbs: " << nb_n << std::endl;
		VectorT dhdx = VectorT::Zero();
		for (int k = 0; k < nb_n; k++) {
			int j = now_nbs[k];
			VectorT lr_ij = surface.Project_To_TPlane(particles.X(i) - particles.X(j), particles.E(i));
			real h_j = surface.Project_To_TPlane_H(particles.X(j) - particles.X(i), particles.E(i)); real h_i = (real)0;
			dhdx += (h_i * pow(particles.SA(i), 2) + h_j * pow(particles.SA(j), 2)) * kernel->Grad_Spiky(lr_ij);
		}
		particles.G(i) = surface.Metric_Tensor(dhdx);
	}
}


template<int d>
void FluidSPHBubble<d>::Update_Phis(void)
{
	phis.resize(particles.Size());
	bnd_normals.resize(particles.Size());
#pragma omp parallel for
	for (int i = 0; i < particles.Size(); i++) {
		real phi; VectorD normal;
		bool flg = analytical_boundary.Get_Nearest_Boundary(particles.X(i), phi, normal);
		phis[i] = phi;
		bnd_normals[i] = normal;
	}
}




template<int d>
void FluidSPHBubble<d>::Print_Particle_Dynamics(int i)
{
	std::cout << "particle " << i << ",X=(" << particles.X(i).transpose() << "),M=" << particles.M(i) << ",F=(" << particles.F(i).transpose() << "),V=(" << particles.V(i).transpose() << "),vel rate=" << particles.V(i).norm() << "\n";
}

template<int d>
void FluidSPHBubble<d>::Numerical_Check(void)
{
	bool success = true;
	int pn = particles.Size();
	for (int i = 0; i < pn; i++) {
		for (int axis = 0; axis < d; axis++) {
			real val = particles.X(i)[axis];
			if (std::isnan(val) || !std::isfinite(val)) {
				std::cerr << "Error: Numerical Check fails for X(" << i << ")[" << axis << "]=" << val << "\n";
				success = false;
			}
		}
	}
	if (!success) {
		exit(0);
	}
}



template<int d>
void FluidSPHBubble<d>::Numerical_Check_SA(void)
{
	bool success = true;
	int pn = particles.Size();
	for (int i = 0; i < pn; i++) {
		real val = particles.SA(i);
		if (std::isnan(val) || !std::isfinite(val)) {
			std::cerr << "Error: Numerical Check fails for SA(" << i << ")" << val << "\n";
			success = false;
		}
	}
	if (!success) {
		exit(0);
	}
}

template<int d>
void FluidSPHBubble<d>::Check_Neighbor_Num(void)
{
	bool success = true;
	for (int i = 0; i < particles.Size(); i++) {
		auto nbs = surface.Tangential_Neighbor_Of(i);
		if (nbs.size() == 0) {
			std::cerr << "Error: Neighbor Number of " << i << " is 0\n";
			success = false;
		}
	}
	if (!success) exit(0);
}

template<int d>
void FluidSPHBubble<d>::Numerical_Check_Force()
{
	bool success = true;
	int pn = particles.Size();
	for (int i = 0; i < pn; i++) {
		for (int axis = 0; axis < d; axis++) {
			real val = particles.F(i)[axis];
			if (std::isnan(val) || !std::isfinite(val)) {
				std::cerr << "Error: Numerical Check fails for F(" << i << ")[" << axis << "]=" << val << "\n";
				success = false;
			}
		}
	}
	if (!success) {
		exit(0);
	}
}

template<int d>
void FluidSPHBubble<d>::Resize_Debug_Structure() {
	std::cout << "diagnosis: " << diagnosis << std::endl;
	if (diagnosis) {
		total_vis_force.resize(particles.Size());
		total_ext_force.resize(particles.Size());
		total_pressure_force.resize(particles.Size());
		total_height_pressure_force.resize(particles.Size());
		total_height_laplacian_pressure_force.resize(particles.Size());
		total_div_pressure_force.resize(particles.Size());
		total_bnd_pressure_force.resize(particles.Size());
		total_marangoni_force.resize(particles.Size());
	}
}


template class MirrorCompensator<2>;
template class MirrorCompensator<3>;

template class FluidSPHBubble<2>;
template class FluidSPHBubble<3>;
