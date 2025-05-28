#include "FluidSPHBubbleDriver.h"
#include "MeshFunc.h"
#include "PointSetFunc.h"
#include "ArrayIO.h"

real Center_Smooth_Kernel(real r, real h, real mtp = 2) {//apply in [-h,h], it's an even function
	return std::min(cos(r / h * pi / 2.0) * mtp, 1.0);
}

template<int d>
void FluidSPHBubbleDriver<d>::Write_Scalar_Field(std::string file_name, const Array<real>& arr, const real scale)
{
	/*int pn = fluid.particles.Size(); Array<Vector<real, d> > normals(pn);
	for (int i = 0; i < pn; i++) {
		normals[i] = fluid.particles.Normal(i) * arr[i] * scale;
	}
	Write_Segments_To_File_3d_Fast<d, real>(fluid.particles.XRef(), normals, file_name);*/
	Array<VectorD> xs, normals;
	int pn = fluid.particles.Size();
	for (int i = 0; i < pn; i++) {
		xs.push_back(fluid.particles.X(i));
		//xs.push_back(fluid.particles.X(i));
		VectorD n_pos = fluid.particles.Normal(i) * arr[i] * scale * 0.5;
		
		//////TO SHOW no BOUNDARY PARTICLES
		if (fluid.particles.Is_Boundary(i)) {
			n_pos *= 0;
		}

		normals.push_back(n_pos);
		//normals.push_back(-n_pos);
	}
	Write_Segments_To_File_3d_Fast<d, real>(xs, normals, file_name);
}



template<int d>
void FluidSPHBubbleDriver<d>::Save_Snapshot(const int frame)
{
	std::cout << "save snapshot for frame " << frame << "\n"; fluid.Numerical_Check_Force();
	std::string snapshot_dir = output_dir + "/snapshot/" + std::to_string(frame);
	std::string snapshot_dir_3d = snapshot_dir + "/3d";
	if (!File::Directory_Exists(snapshot_dir.c_str())) File::Create_Directory(snapshot_dir);
	if (!File::Directory_Exists(snapshot_dir_3d.c_str())) File::Create_Directory(snapshot_dir_3d);
	fluid.particles.Save_Snapshot(snapshot_dir);
	if (fluid.fluid_3d != nullptr)fluid_3d.particles.Save_Snapshot(snapshot_dir + "/3d");
	if (fluid.n_pressure_params.air_pressure_mode == "ib") fluid.air_solver.Save_Snapshot(snapshot_dir);
}

template<int d>
void FluidSPHBubbleDriver<d>::Load_Snapshot(const int frame)
{
	std::string snapshot_dir = output_dir + "/snapshot/" + std::to_string(frame);
	fluid.particles.Load_Snapshot(snapshot_dir);
	if (fluid.fluid_3d != nullptr) fluid_3d.particles.Load_Snapshot(snapshot_dir + "/3d");
	if (fluid.n_pressure_params.air_pressure_mode == "ib") fluid.air_solver.Load_Snapshot(snapshot_dir);
}

template<int d>
void FluidSPHBubbleDriver<d>::Write_Air_Solver(const int frame)
{
	if (frame == 0) {
		fluid.air_solver.mac_grid.grid.Write_To_File_3d(frame_dir + "/grid");
		fluid.air_solver.bc.Write_Psi_D_To_File_3d(frame_dir + "/psi_D");
		fluid.air_solver.bc.Write_Psi_N_To_File_3d(frame_dir + "/psi_N");
	}
	fluid.air_solver.velocity.Write_To_File_3d(frame_dir + "/velocity");
}

template<int d>
void FluidSPHBubbleDriver<d>::Write_Output_Files(const int frame)
{
	static double begin_time = omp_get_wtime();
	static double last_time = begin_time;
	static double now_time = last_time;
	//std::cout << "suspended\n"; return;
	Base::Write_Output_Files(frame);

	std::cout << "#  Write Frame " << frame << " to: " << Base::frame_dir << std::endl;

	BinaryDataIO::Write_Vector_Array_3D<real, d>(frame_dir + "/x_bin", fluid.particles.XRef());
	Array<VectorD> normals; normals.resize(fluid.particles.Size());
	for (int i = 0; i < normals.size(); i++) normals[i] = fluid.particles.Normal(i);
	BinaryDataIO::Write_Vector_Array_3D<real, d>(frame_dir + "/norm_bin", normals);
	BinaryDataIO::Write_Scalar_Array(frame_dir + "/h_bin", fluid.particles.HRef());
	BinaryDataIO::Write_Scalar_Array(frame_dir + "/rh_bin", fluid.particles.RHRef());
	BinaryDataIO::Write_Scalar_Array(frame_dir + "/kh_bin", fluid.particles.KHRef());

	//output ib things
	if (fluid.n_pressure_params.air_pressure_mode == "ib") Write_Air_Solver(frame);

	fluid.particles.Write_To_File_3d_Fast(frame_dir + "/tracker_points"); 
	Write_Segments_To_File_3d_Fast<d, real>(fluid.particles.XRef(), fluid.particles.VRef(), frame_dir + "/point_velocity"); 
	PointSetFunc::Write_Tracker_Circles_To_File<d>(frame_dir + "/tracker_circles", fluid.particles); 
	Write_Scalar_Field(frame_dir + "/point_height", fluid.particles.HRef(), 0.2 / fluid.thickness * fluid.simulation_scale);
	//Write_Scalar_Field(frame_dir + "/point_gamma", fluid.particles.ConcRef(), 1);
	//Write_Scalar_Field(frame_dir + "/point_gamma", fluid.particles.RHRef(), 0.2 / fluid.thickness * fluid.simulation_scale);
	Write_Scalar_Field(frame_dir + "/point_gamma", fluid.particles.KHRef(), 1000000);
	
	if (save_all) {
		Write_Segments_To_File_3d_Fast<d, real>(fluid.particles.XRef(), fluid.particles.FRef(), frame_dir + "/point_force");
		Write_Scalar_Field(frame_dir + "/point_sa", fluid.particles.SARef(), 2);
		PointSetFunc::Write_Local_Frames_To_File(frame_dir + "/segment_mesh", fluid.particles, 0.002);
	}

	if (fluid.fluid_3d != nullptr) {
		{std::string file_name = frame_dir + "/tracker_points_3d";
		fluid_3d.particles.Write_To_File_3d_Fast(file_name); }
		BinaryDataIO::Write_Scalar_Array(frame_dir + "/rh_bin_3d", fluid_3d.particles.RHRef());
		BinaryDataIO::Write_Vector_Array(frame_dir + "/norm_bin_3d", fluid_3d.particles.SNRef());
		//{std::string file_name = frame_dir + "/point_force";
		//Write_Segments_To_File_3d_Fast<d, real>(fluid_3d.particles.XRef(), fluid_3d.particles.FRef(), file_name); }
		//{std::string file_name = frame_dir + "/point_velocity";
		//Write_Segments_To_File_3d_Fast<d, real>(fluid_3d.particles.XRef(), fluid_3d.particles.VRef(), file_name); }
	}

	BinaryDataIO::Write_Scalar_Array<real>(frame_dir + "/sa_bin", fluid.particles.SARef());

	if (frame % snapshot_stride == (int)0) Save_Snapshot(frame);

	now_time = omp_get_wtime();
	if (verbose) std::cout << "#  Time: "
		<< std::setw(9) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << (now_time - last_time) << "s / "
		<< std::setw(9) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << (now_time - begin_time) << "s" << std::endl;
	if (verbose) std::cout << "   ETA:  "
		<< std::setw(9) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << (now_time - begin_time) / frame * (last_frame - frame) << "s / "
		<< std::setw(9) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << (now_time - begin_time) / frame * last_frame << "s"
		<< std::endl << std::endl;
	std::cout << std::defaultfloat;
	last_time = now_time;

}

template<int d>
std::function<Vector<real,d>(const Vector<real, d>&)> FluidSPHBubbleDriver<d>::Corridor_Flow_Func(real R, real line_vel)
{
	std::function<VectorD(const Vector<real, d>&)> vel_func = nullptr;
	if constexpr (d == 3) {
		const static int start_frame = 0, end_frame = 100;
		const static VectorD &unit_y = VectorD::Unit(1);
		const static real half_h = R / 6.0;
		const static real omega_rate = line_vel / R;
		vel_func = [&](const VectorD &r)->VectorD {
			if (start_frame <= current_frame && current_frame < end_frame) {
				//if (fabs(r[1]) > half_h || r[2] < 0) return VectorD::Zero();
				if (fabs(r[1]) > half_h) return VectorD::Zero();
				real partition_val = (current_frame - start_frame + 0.0) / (end_frame - start_frame + 0.0) * pi;
				real w = std::min(sin(partition_val) * 2.0, 1.0);
				real wh = Center_Smooth_Kernel(r[1], half_h);
				VectorD omega = omega_rate * unit_y * wh;
				return omega.cross(r) * w;
			}
			else return VectorD::Zero();
		};
	}
	return vel_func;
}

template<int d>
void FluidSPHBubbleDriver<d>::Seed_Vortex(int seed_num, real max_strength, real sample_r)
{
	static std::shared_ptr<RandomInt> rand_id = std::make_shared<RandomInt>(0, fluid.particles.Size() - 1);
	static std::shared_ptr<RandomNumber> rand_w = std::make_shared<RandomNumber>(-1, 1);
	for (int i = 0; i < 10000; i++) rand_id->Value();
#pragma omp parallel for
	for (int i = 0; i < fluid.particles.Size(); i++) {
		fluid.particles.Vrt(i) = 0;
	}
	for (int sid = 0; sid < seed_num; sid++) {
		int id = rand_id->Value();
		const real vorticity = max_strength * rand_w->Value();
		const VectorD& pos = fluid.particles.X(id);
#pragma omp parallel for
		for (int i = 0; i < fluid.particles.Size(); i++) {
			real r = (fluid.particles.X(i) - pos).norm();
			real w = Center_Smooth_Kernel(r, sample_r, sqrt(2.0));
			fluid.particles.Vrt(i) += w * vorticity;
		}
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Seed_Vortex_Rand(int seed_num, real max_strength, real sample_r)
{
	static std::shared_ptr<RandomInt> rand_id = std::make_shared<RandomInt>(0, fluid.particles.Size() - 1);
	static std::shared_ptr<RandomNumber> rand_w = std::make_shared<RandomNumber>(-1, 1);
#pragma omp parallel for
	for (int i = 0; i < fluid.particles.Size(); i++) {
		fluid.particles.Vrt(i) = 0;
	}
	for (int sid = 0; sid < seed_num; sid++) {
		real h = (0.8 * Rand_Number() + 0.2) * sample_r;
		int id = rand_id->Value();
		const real vorticity = max_strength * rand_w->Value();
		const VectorD& pos = fluid.particles.X(id);
#pragma omp parallel for
		for (int i = 0; i < fluid.particles.Size(); i++) {
			real r = (fluid.particles.X(i) - pos).norm();
			real w = Center_Smooth_Kernel(r, h, sqrt(2.0)+0.3*rand_w->Value());
			fluid.particles.Vrt(i) += w * vorticity;
		}
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Initialize_Center_Oil(void)
{
	for (int i = 0; i < fluid.particles.Size(); i++) {
		Vector2 tmp;
		tmp << fluid.particles.X(i)[0], fluid.particles.X(i)[2];
		if (tmp.norm() < 0.5) {
			fluid.particles.Conc(i) = 0.8;
		}
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Set_Fixed_Gamma(real gamma)
{
	//gamma_water - gamma_soap * concentration = gamma
	real conc = (fluid.gamma_water - gamma) / fluid.gamma_soap;
	for (int i = 0; i < fluid.particles.Size(); i++) {
		fluid.particles.Conc(i) = conc;
	}
	fluid.diffuse_soap = false;
}

template<int d>
void FluidSPHBubbleDriver<d>::Set_Physical_Parameters(void)
{
	//Conc in [0,1]
	fluid.gamma_water = 72e-3;//72mN/m
	//fluid.gamma_soap = 42e-3;//limitation of surface concentration of soap is 30mN/m, with Conc=1
	fluid.gamma_soap = 70e-3;//limitation of surface concentration of soap is 30mN/m, with Conc=1
	fluid.n_pressure_params.atmo_pressure = 1.01e5;//101 kpa
	fluid.viscosity_water = 1.005e-3;
	fluid.g = VectorD::Unit(1) * (-9.8);
}

template<int d>
real FluidSPHBubbleDriver<d>::Perlin_Noise(VectorD pos, std::uint32_t seed, real perlin_freq, real perlin_scale, std::uint32_t octaves)
{
	/// 
	/// Return a perlin noise value.
	/// 
	/// pos : position of a particle.
	/// seed : A random seed.
	/// perlin_freq : The larger this parameter, the more frequent the noise.
	/// perlin_scale : The scale of offset from the original position. 
	/// octaves : Number of octaves in perlin noise, default value is 4.
	/// noise_value : Return value of perlin noise generator, range in [-1, 1].
	/// 
	const siv::PerlinNoise perlin(seed);
	real noise_value;
	if constexpr (d == 3) {
		noise_value = perlin.accumulatedOctaveNoise3D(pos[0] * perlin_freq, pos[1] * perlin_freq, pos[2] * perlin_freq, octaves);
	}
	else {
		noise_value = perlin.accumulatedOctaveNoise2D(pos[0] * perlin_freq, pos[1] * perlin_freq, octaves);
	}
	return noise_value * perlin_scale;
}


template<int d>
void FluidSPHBubbleDriver<d>::Case_1(void) {//Perlin Noise Test, Same as case 2
	int sphere_sub_num = 5;//5:10242, 6:40962, 7: 163842
	real R = 0.133;//5:0.133, 7:0.5320

	std::cout << "Enter initialization Case 2\n";

	max_iter_per_frame = -1;
	cfl = 0.1;
	frame_rate = 50;
	real dx = 0.005;
	int number_2d = 1000;//number of particles for 2d
	VectorD ctr = VectorD::Zero();
	if constexpr (d == 2) {
		R = dx * number_2d / (2 * pi);
		dx = PointSetFunc::Initialize_Circle_Points(ctr, R, number_2d, fluid.particles);
	}
	if constexpr (d == 3) { dx = PointSetFunc::Initialize_Sphere_Points(ctr, R, sphere_sub_num, fluid.particles); }
	fluid.simulation_scale = R;
	std::cout << "Initialize with " << fluid.particles.Size() << " particles and dx = " << dx << " and R = " << R << "\n";

	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "init R: " << R << "\n";
		ext_force = [&](const int i)->VectorD {
			const VectorD& pos = fluid.particles.X(i);
			const real r = fluid.simulation_scale * 0.25;
			//const real r = 0.198944;
			if (pos[0] * pos[0] <= r * r && current_frame < 10) {
				return -fluid.particles.X(i) * fluid.particles.M(i);
			}
			else return VectorD::Zero();
		};
	}
	else if constexpr (d == 3) {
		ext_force = [&](const int i)->VectorD {
			const VectorD& pos = fluid.particles.X(i);
			const real r = fluid.simulation_scale * 0.25;
			real sqr_tan_len = pos[0] * pos[0] + pos[2] * pos[2];
			VectorD force = VectorD::Zero();

			if (sqr_tan_len <= r * r && current_frame < 50) {
				VectorD ext_norm = -fluid.particles.X(i) * fluid.particles.M(i);
				force += ext_norm * Center_Smooth_Kernel(sqrt(sqr_tan_len), r) * 20;
			}
			return force;
		};
	}

	Set_Physical_Parameters();
	real rho = 1e3; fluid.thickness = 5e-7;//10um
	real default_conc = 1;//"extreme soap"
	real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
	real total_vol = total_surface * fluid.thickness, total_mass = total_vol * rho;
	int n = fluid.particles.Size();
	std::cout << "Particle Vol: " << total_vol / n << ", Particle Mass: " << total_mass / n << "\n";
	for (int i = 0; i < n; i++) {
		Initialize_Particle(i, total_mass / n, fluid.thickness);
		fluid.particles.B(i) = 0;
		fluid.particles.Vol(i) = total_vol / n;
		fluid.particles.Conc(i) = default_conc;
		//if (i == 0) { fluid.particles.V(i) = -fluid.particles.Normal(i)*100; }	

		//Perlin noise test
		fluid.particles.X(i) += fluid.particles.Normal(i) * Perlin_Noise(fluid.particles.X(i),123,4.0,0.05);
	}

	real V = (d == 2) ? R * R * pi : 4.0 / 3.0 * pi * R * R * R;
	std::cout << "Analytical volume: " << V << "\n";
	real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
	real alpha_vel = 0.005 * 50;//move 0.005 in one frame with FPS=50

	//std::ifstream fin(scene_file_name);
	//if (!fin.is_open()) { std::cerr << "Driver Case_2 error: cannot open " << scene_file_name << "\n"; exit(0); }

	//world forces
	fluid.gravity_coeff = 0;
	//fluid.friction_coeff = 1, fluid.air_velocity_func = Corridor_Flow_Func(R, alpha_vel);
	fluid.external_force_func = ext_force;

	//tangential force
	fluid.t_pressure_params.Set_Baseline(rho, fluid.gamma_water, fluid.thickness, frame_rate);
	fluid.t_pressure_params.boundary_pressure_coeff = 0;
	fluid.divergence_params.calculate_field = "all";
	fluid.viscosity_coeff = 0.2 * alpha_vis;
	fluid.marangoni_coeff = 0;
	//real vis_term; fin >> vis_term; fluid.viscosity_coeff = vis_term * alpha_vis; std::cout << "vis term: " << vis_term << "\n";

	//vorticity params
	//fluid.vorticity_params.weight_kernel = KernelType::QUINTIC;
	//fluid.vorticity_params.Initialize_From_File(fin);
	//fluid.vorticity_params = VorticityParams(KernelType::QUINTIC, 1, true, 1e-4);
	fluid.vorticity_params = VorticityParams(KernelType::QUINTIC, 0, true, 1e-4);
	//Seed_Vortex(10, 0, dx * 0.5);

	fluid.rh_params = RenderHeightParams("laplacian", true, false, 0);

	fluid.normal_viscosity_coeff = 1 * 0;
	fluid.boundary_params.Set_No_Boundary();

	fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
	fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);

	//fin.close();
	fluid.max_vel = alpha_vel;
	fluid.Initialize(dx);

	for (int i = 0; i < fluid.particles.Size(); i++) fluid.particles.RH(i) = fluid.particles.H(i);

	real V_num = fluid.Compute_Enclosed_Volume();
	fluid.n_pressure_params.Set_Sphere_Baseline(d, R, fluid.Surface_Tension_Coefficient(default_conc), V_num, rho, fluid.thickness);


}

template<int d>
void FluidSPHBubbleDriver<d>::Case_2(void)
{


	std::cout << "Enter initialization Case 2\n";

	int sphere_sub_num = 7;//5:10242, 6:40962, 7: 163842
	real R = 0.5320;//5:0.133, 7:0.5320

	real ext_term = 0, vor_term = 0;
	std::ifstream fin(scene_file_name);
	if (!fin.is_open()) { std::cerr << "Driver Case_2 error: cannot open " << scene_file_name << "\n"; exit(0); }
	fin >> ext_term >> vor_term;
	fin.close();

	max_iter_per_frame = -1;
	cfl = 0.1;
	frame_rate = 50;
	real dx=0.005;
	int number_2d = 1000;//number of particles for 2d
	VectorD ctr = VectorD::Zero();
	if constexpr (d == 2) {
		R = dx * number_2d / (2 * pi);
		dx = PointSetFunc::Initialize_Circle_Points(ctr, R, number_2d, fluid.particles);
	}
	if constexpr (d == 3) { dx = PointSetFunc::Initialize_Sphere_Points(ctr, R, sphere_sub_num, fluid.particles); }
	fluid.simulation_scale = R;
	std::cout << "Initialize with " << fluid.particles.Size() << " particles and dx = " << dx << " and R = " << R << "\n";

	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "init R: " << R << "\n";
		ext_force = [&](const int i)->VectorD {
			const VectorD& pos = fluid.particles.X(i);
			const real r = fluid.simulation_scale * 0.25;
			//const real r = 0.198944;
			if (pos[0] * pos[0] <= r * r && current_frame < 10) {
				return -fluid.particles.X(i) * fluid.particles.M(i);
			}
			else return VectorD::Zero();
		};
	}
	else if constexpr (d == 3) {
		const static real ext_coeff = ext_term;
		ext_force = [&](const int i)->VectorD {
			const VectorD& pos = fluid.particles.X(i);
			const real r = fluid.simulation_scale * 0.25;
			real sqr_tan_len = pos[0] * pos[0] + pos[2] * pos[2];
			VectorD force = VectorD::Zero();

			if (sqr_tan_len <= r * r && current_frame < 50) {
				VectorD ext_norm = -fluid.particles.X(i) * fluid.particles.M(i);
				force += ext_norm * Center_Smooth_Kernel(sqrt(sqr_tan_len), r) * ext_coeff;
			}
			return force;
		};
	}

	Set_Physical_Parameters();
	real rho = 1e3; fluid.thickness = 5e-7;//10um
	real default_conc = 1;//"extreme soap"
	real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
	real total_vol = total_surface * fluid.thickness, total_mass = total_vol * rho;
	int n = fluid.particles.Size();
	std::cout << "Particle Vol: " << total_vol / n << ", Particle Mass: " << total_mass / n << "\n";
	for (int i = 0; i < n; i++) {
		Initialize_Particle(i, total_mass / n, fluid.thickness);
		fluid.particles.B(i) = 0;
		fluid.particles.Vol(i) = total_vol / n;
		fluid.particles.Conc(i) = default_conc;
	}
	
	real V = (d == 2) ? R * R * pi : 4.0 / 3.0 * pi * R * R * R;
	std::cout << "Analytical volume: " << V << "\n";
	real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
	real alpha_vel = 0.005 * 50;//move 0.005 in one frame with FPS=50

	//world forces
	fluid.gravity_coeff = 0;
	//fluid.friction_coeff = 1, fluid.air_velocity_func = Corridor_Flow_Func(R, alpha_vel);
	fluid.external_force_func = ext_force;

	//tangential force
	fluid.t_pressure_params.Set_Baseline2(rho, fluid.gamma_water, fluid.thickness, frame_rate, R);
	fluid.t_pressure_params.boundary_pressure_coeff = 0;
	fluid.divergence_params.calculate_field = "all";
	fluid.viscosity_coeff = 0.2 * alpha_vis;
	fluid.marangoni_coeff = 0;

	//vorticity params
	fluid.vorticity_params = VorticityParams(KernelType::QUINTIC, 1, true, 1e-4);
	Seed_Vortex(100, vor_term, dx * 3);

	fluid.rh_params = RenderHeightParams("laplacian", true, false, 0);

	fluid.normal_viscosity_coeff = 1*0;
	fluid.boundary_params.Set_No_Boundary();

	fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
	fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);

	//fin.close();
	fluid.max_vel = alpha_vel;
	fluid.Initialize(dx);

	for (int i = 0; i < fluid.particles.Size(); i++) fluid.particles.RH(i) = fluid.particles.H(i);

	real V_num = fluid.Compute_Enclosed_Volume();
	fluid.n_pressure_params.Set_Sphere_Baseline(d, R, fluid.Surface_Tension_Coefficient(default_conc), V_num, rho, fluid.thickness);
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_3(void)	////Load obj test
{
	std::cout << "Enter initialization Case 3\n";

	max_iter_per_frame = -1;
	cfl = 0.1;
	frame_rate = 50;
	real dx = 0.005, R;

	if constexpr (d == 3) {

		std::string bubble_surface_mesh_file_name = Path::Data() + "/meshes/bubbles/bubble1.obj";
		Array<std::shared_ptr<TriangleMesh<3>>> bubble_obj;
		Obj::Read_From_Obj_File<TriangleMesh<3>>(bubble_surface_mesh_file_name, bubble_obj);
		if (!bubble_obj.size()) { std::cerr << "case 3 cannot read\n"; exit(0); }

		auto bubble_mesh = bubble_obj[0];
		MeshFunc::Rescale<3>(bubble_mesh->Vertices(), (real)1);
		real mesh_dx = MeshFunc::Average_Edge_Length<3>(bubble_mesh->Vertices(), bubble_mesh->Elements());
		real multiplier = dx / mesh_dx;
		fluid.particles.Resize((int)bubble_mesh->Vertices().size());
		fluid.particles.XRef() = bubble_mesh->Vertices();
		//std::cout << "normal: " << (*(bubble_mesh->Normals())).empty() << std::endl;
		for (int i = 0; i < fluid.particles.Size(); i++) {
			Vector3 normal = (*(bubble_mesh->Normals()))[i];
			Vector3 t1 = -PointSetFunc::Orthogonal_Vector(normal);
			Vector3 t2 = t1.cross(normal);
			fluid.particles.E(i).col(0) = t1;
			fluid.particles.E(i).col(1) = t2;
			fluid.particles.E(i).col(2) = normal;
			fluid.particles.X(i) *= multiplier;
		}
		VectorD x_center = AuxFunc::Mean(fluid.particles.XRef());
		real sum_norm = 0; for (int i = 0; i < fluid.particles.Size(); i++) sum_norm += (fluid.particles.X(i) - x_center).norm();
		R = sum_norm / (fluid.particles.Size() + 0.0);//mean norm
		fluid.simulation_scale = R;
		std::cout << "Initialize with " << fluid.particles.Size() << " particles and dx = " << dx << " and R = " << R << "\n";
	}

	Set_Physical_Parameters();
	real rho = 1e3; fluid.thickness = 5e-7;//10um
	real default_conc = 1;//"extreme soap"
	real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
	real total_vol = total_surface * fluid.thickness, total_mass = total_vol * rho;
	int n = fluid.particles.Size();
	std::cout << "Particle Vol: " << total_vol / n << ", Particle Mass: " << total_mass / n << "\n";
	for (int i = 0; i < n; i++) {
		Initialize_Particle(i, total_mass / n, fluid.thickness);
		fluid.particles.B(i) = 0;
		fluid.particles.Vol(i) = total_vol / n;
		fluid.particles.Conc(i) = default_conc;
		//if (i == 0) { fluid.particles.V(i) = -fluid.particles.Normal(i)*100; }	
	}

	real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
	real alpha_vel = 0.005 * 50;//move 0.005 in one frame with FPS=50

	//world forces
	fluid.gravity_coeff = 0;

	//tangential force
	fluid.t_pressure_params.Set_Baseline2(rho, fluid.gamma_water, fluid.thickness, frame_rate, R);
	fluid.t_pressure_params.height_pressure_coeff *= 10;
	fluid.t_pressure_params.boundary_pressure_coeff = 0;
	fluid.divergence_params.calculate_field = "all";
	fluid.viscosity_coeff = 0.2 * alpha_vis;
	fluid.marangoni_coeff = 0;

	//vorticity params
	fluid.vorticity_params = VorticityParams(KernelType::QUINTIC, 1, true, 1e-4);
	Seed_Vortex(10, 2e-6, dx * 3);

	fluid.rh_params = RenderHeightParams("divergence", true, true, 1e-4);

	fluid.normal_viscosity_coeff = 1 * 0;
	fluid.boundary_params.Set_No_Boundary();
	fluid.delete_solitary = false;

	fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
	fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);

	fluid.max_vel = alpha_vel;
	fluid.Initialize(dx);

	for (int i = 0; i < n; i++)fluid.particles.RH(i) = fluid.particles.H(i);
	real V_num = fluid.Compute_Enclosed_Volume();
	fluid.n_pressure_params.Set_Irregular(d, fluid.Surface_Tension_Coefficient(default_conc), V_num, rho, fluid.thickness);
	fluid.n_pressure_params.pressure_constant *= 1.1;
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_4(void)
{
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_5(void)	////dome
{
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_6(void)
{
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_7(void) ////Improved 2D Disk
{
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_8(void)	//capillary wave on a segment
{
	std::cout << "Enter initialization Case 9\n";
	max_iter_per_frame = -1;
	cfl = 0.1;
	frame_rate = 50;
	real dx = 0.005, R;
	int number_2d = 1000;//number of particles for 2d
	VectorD ctr = VectorD::Zero();
	if constexpr (d == 2) {
		R = dx * number_2d;
		PointSetFunc::Initialize_Segment_Points(ctr - VectorD::Unit(0) * 0.5 * R, ctr + VectorD::Unit(0) * 0.5 * 5, number_2d, fluid.particles);
	}
	fluid.simulation_scale = R;
	std::cout << "Initialize with " << fluid.particles.Size() << " particles and dx = " << dx << " and R = " << R << "\n";

	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "init R: " << R << "\n";
		ext_force = [&](const int i)->VectorD {
			const VectorD& pos = fluid.particles.X(i);
			const real r = fluid.simulation_scale * 0.25;
			if (pos[0] * pos[0] <= r * r && current_frame < 10&&pos[1]>0) {
				return -fluid.particles.X(i) * fluid.particles.M(i);
			}
			else return VectorD::Zero();
		};
	}

	Set_Physical_Parameters();
	real rho = 1e3; fluid.thickness = 5e-7;//10um
	real default_conc = 1;//"extreme soap"
	real total_vol = R * fluid.thickness, total_mass = total_vol * rho;
	int n = fluid.particles.Size();
	std::cout << "Particle Vol: " << total_vol / n << ", Particle Mass: " << total_mass / n << "\n";
	for (int i = 0; i < n; i++) {
		Initialize_Particle(i, total_mass / n, fluid.thickness);
		if (i == 0 || i == n - 1) fluid.particles.B(i) = 1;
		else fluid.particles.B(i) = 0;
		fluid.particles.Vol(i) = total_vol / n;
		fluid.particles.Conc(i) = default_conc;
	}

	real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
	real alpha_vel = 0.005 * 50;//move 0.005 in one frame with FPS=50

	//world forces
	fluid.gravity_coeff = 0;
	//fluid.friction_coeff = 1, fluid.air_velocity_func = Corridor_Flow_Func(R, alpha_vel);
	fluid.external_force_func = ext_force;

	//tangential force
	fluid.t_pressure_params.Set_Baseline2(rho, fluid.gamma_water, fluid.thickness, frame_rate, R);
	fluid.t_pressure_params.boundary_pressure_coeff = 0;
	fluid.divergence_params.calculate_field = "all";
	fluid.viscosity_coeff = 0.2 * alpha_vis;
	fluid.marangoni_coeff = 0;

	//vorticity params
	fluid.vorticity_params = VorticityParams(KernelType::QUINTIC, 0, false, 0);

	fluid.rh_params = RenderHeightParams("laplacian", true, true, 1e-4);

	fluid.normal_viscosity_coeff = 1 * 0;
	fluid.boundary_params.Set_No_Boundary();

	fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
	fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);

	fluid.max_vel = alpha_vel;
	fluid.Initialize(dx);

	fluid.n_pressure_params.capillary_coeff = 1;
	fluid.n_pressure_params.air_pressure_mode = "none";
	fluid.n_pressure_params.opt_mode = OperatorParams("only_fluid", KernelType::SPIKY);

	for (int i = 0; i < fluid.particles.Size(); i++) fluid.particles.RH(i) = fluid.particles.H(i);
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_9(void)
{
	std::cout << "Enter initialization Case 9\n";
	max_iter_per_frame = -1;
	cfl = 0.1;
	frame_rate = 50;
	real dx = 0.005, R;
	int number_2d = 1000;//number of particles for 2d
	VectorD ctr = VectorD::Zero();
	if constexpr (d == 2) {
		R = dx * number_2d / (2 * pi);
		dx = PointSetFunc::Initialize_Circle_Points(ctr, R, number_2d, fluid.particles);
	}
	fluid.simulation_scale = R;
	std::cout << "Initialize with " << fluid.particles.Size() << " particles and dx = " << dx << " and R = " << R << "\n";

	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "init R: " << R << "\n";
		ext_force = [&](const int i)->VectorD {
			const VectorD& pos = fluid.particles.X(i);
			real omega = atan2(pos[1], pos[0]) - pi / 2;
			//const real r = fluid.simulation_scale * 0.25;
			//const real r = 0.198944;
			//if (pos[0] * pos[0] <= r * r && current_frame < 10 && pos[1]>0) {
			if(pos[1]>0&&current_frame<10&&fabs(omega)<pi/8){
				return -fluid.particles.X(i) * fluid.particles.M(i)*0.2;
			}
			else return VectorD::Zero();
		};
	}

	Set_Physical_Parameters();
	real rho = 1e3; fluid.thickness = 5e-7;//10um
	real default_conc = 1;//"extreme soap"
	real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
	real total_vol = total_surface * fluid.thickness, total_mass = total_vol * rho;
	int n = fluid.particles.Size();
	std::cout << "Particle Vol: " << total_vol / n << ", Particle Mass: " << total_mass / n << "\n";
	for (int i = 0; i < n; i++) {
		Initialize_Particle(i, total_mass / n, fluid.thickness);
		fluid.particles.B(i) = 0;
		fluid.particles.Vol(i) = total_vol / n;
		fluid.particles.Conc(i) = default_conc;
	}

	real V = (d == 2) ? R * R * pi : 4.0 / 3.0 * pi * R * R * R;
	std::cout << "Analytical volume: " << V << "\n";
	real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
	real alpha_vel = 0.005 * 50;//move 0.005 in one frame with FPS=50

	//world forces
	fluid.gravity_coeff = 0;
	//fluid.friction_coeff = 1, fluid.air_velocity_func = Corridor_Flow_Func(R, alpha_vel);
	fluid.external_force_func = ext_force;

	//tangential force
	fluid.t_pressure_params.Set_Baseline2(rho, fluid.gamma_water, fluid.thickness, frame_rate, R);
	fluid.t_pressure_params.boundary_pressure_coeff = 0;
	fluid.t_pressure_params.tangential_pressure_coeff = 0;//
	fluid.divergence_params.calculate_field = "all";
	fluid.viscosity_coeff = 0.2 * alpha_vis;
	fluid.marangoni_coeff = 0;

	//vorticity params
	fluid.vorticity_params = VorticityParams(KernelType::QUINTIC, 0, true, 1e-4);

	fluid.rh_params = RenderHeightParams("laplacian", true, false, 0);

	fluid.normal_viscosity_coeff = 1 * 0;
	fluid.boundary_params.Set_No_Boundary();

	fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
	fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);

	//fin.close();
	fluid.max_vel = alpha_vel;
	fluid.Initialize(dx);

	for (int i = 0; i < fluid.particles.Size(); i++) fluid.particles.RH(i) = fluid.particles.H(i);

	real V_num = fluid.Compute_Enclosed_Volume();
	fluid.n_pressure_params.Set_Sphere_Baseline(d, R, fluid.Surface_Tension_Coefficient(default_conc), V_num, rho, fluid.thickness, 1.0 / 9);
	std::cout << "surface tension: " << fluid.n_pressure_params.capillary_coeff * fluid.Surface_Tension_Coefficient(default_conc) << "\n";

	fluid.output_file.open(output_dir + "/peaks.txt");


	real predicted_wave_speed = sqrt(fluid.n_pressure_params.capillary_coeff * 4 * fluid.gamma_soap / (rho * fluid.thickness));
	std::cout << "predicted wave speed: " << predicted_wave_speed << "\n";
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_10(void) ////2D Disk in 3D but with grid particles
{
	if constexpr (d == 2) {
		std::cerr << "Case_10 error: d==2 not supported\n" << std::endl;
	}
	max_iter_per_frame = -1;
	frame_rate = 5.;
	cfl = 0.1;
	fluid.max_vel = 0.1;

	real dx = 0.005;
	//int side_num = 100;
	int side_num = 240;
	//int side_num = 200;
	fluid.simulation_scale = side_num * dx;
	fluid.np_on_h = 6;
	//int cycle = 100.;	

	std::vector<int> is_boundary;
	if constexpr (d == 3) {
		VectorD domain_min = AuxFunc::V<d>(-fluid.simulation_scale, 0, -fluid.simulation_scale) * 0.5;
		is_boundary = PointSetFunc::Initialize_Lattice_Points(VectorD::Zero(), AuxFunc::Vi<2>(side_num, side_num), VectorD::Unit(0), VectorD::Unit(2), dx, fluid.particles);
	}

	Set_Physical_Parameters();

	real rho = 1e3;
	fluid.thickness = 5e-7;//10um
	real default_conc = 0;//"extreme soap"
	real total_surface = pow(fluid.simulation_scale, 2);
	real total_vol = total_surface * fluid.thickness;
	real total_mass = total_vol * rho;

	int n = fluid.particles.Size();
	for (int i = 0; i < n; i++) {
		Initialize_Particle(i, total_mass / n, fluid.thickness);
		fluid.particles.B(i) = is_boundary[i];
		fluid.particles.Vol(i) = total_vol / n;
		fluid.particles.Conc(i) = default_conc;
		fluid.particles.RH(i) = fluid.thickness;
	}

	//main params
	fluid.default_sim_params = DefaultSimParams(dx, rho, fluid.viscosity_water, fluid.gamma_water, fluid.thickness, fluid.simulation_scale);

	std::function<VectorD(const int)> ext_force = [&](const int i)->VectorD {
		const real g_strength = 1 * 1. * 9.8 * 0.33 / fluid.simulation_scale * 0.6 * 10 * pow(.1 / 0.15, 2) * fluid.default_sim_params.gravity_coeff;
		const real cycle = 400.;
		real alpha = (current_frame + 0.0) / cycle * 2 * pi; // 100 is the cycle
		return AuxFunc::V<d>(cos(alpha), sin(alpha), 0) * g_strength * fluid.particles.M(i);
	};
	//std::function<VectorD(const int)> ext_force = [&](const int i)->VectorD {
	//	const real g_strength = 1 * 1 * 9.8 * 0.33 / fluid.simulation_scale * 0.6 * 10 * pow(.1 / 0.15, 2) * fluid.default_sim_params.gravity_coeff;
	//	const real cycle = 100;
	//	real alpha = (current_frame + 0.0) / cycle * 2 * pi; // 100 is the cycle
	//	return -VectorD::Unit(1) * g_strength * fluid.particles.M(i);
	//};
	//world force
	fluid.gravity_coeff = 0;
	fluid.external_force_func = ext_force;
	//tangential
	fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, 50);
	fluid.t_pressure_params.height_pressure_coeff *= 2;
	fluid.t_pressure_params.laplacian_pressure_coeff *= 2;
	fluid.viscosity_coeff = 2 * fluid.default_sim_params.viscosity_coeff;
	fluid.marangoni_coeff = 0;
	//normal
	fluid.n_pressure_params.Set_Circle_Baseline(fluid.simulation_scale, fluid.gamma_water, rho, fluid.thickness);
	fluid.n_pressure_params.capillary_coeff *= 0.0015;
	fluid.normal_viscosity_coeff = 0;
	//boundary
	fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
	fluid.boundary_params.replenish = false;
	fluid.analytical_boundary.Add_Obstacle(std::make_shared<Plane<d>>(VectorD::Unit(0), VectorD::Unit(0) * (-0.5 * fluid.simulation_scale)));
	fluid.analytical_boundary.Add_Obstacle(std::make_shared<Plane<d>>(-VectorD::Unit(0), VectorD::Unit(0) * (0.5 * fluid.simulation_scale)));
	fluid.analytical_boundary.Add_Obstacle(std::make_shared<Plane<d>>(VectorD::Unit(2), VectorD::Unit(2) * (-0.5 * fluid.simulation_scale)));
	fluid.analytical_boundary.Add_Obstacle(std::make_shared<Plane<d>>(-VectorD::Unit(2), VectorD::Unit(2) * (0.5 * fluid.simulation_scale)));
	//rh
	fluid.rh_params = RenderHeightParams("divergence", false, true, 1e-4);
	//operators
	fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
	fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);
	fluid.geometry_params = OperatorParams("only_fluid", KernelType::SPIKY);

	////add perlin noise to stuff
	std::uint32_t seed = 123;
	real perlin_freq;
	std::uint32_t octaves = 2;
	const siv::PerlinNoise perlin(seed);
	real perlin_inc;

	if (first_frame < 1) {
		for (int i = 0; i < fluid.particles.Size(); i++) {
			perlin_freq = 20.;
			perlin_inc = perlin.accumulatedOctaveNoise3D(fluid.particles.X(i)[0] * perlin_freq, fluid.particles.X(i)[1] * perlin_freq, fluid.particles.X(i)[2] * perlin_freq, octaves);
			real perlin_multiplier_V = std::max(0.1, 1 + (perlin_inc * 1.));
			//real perlin_multiplier_V = 1.;
			fluid.particles.Vol(i) *= perlin_multiplier_V;
			fluid.particles.RH(i) *= perlin_multiplier_V;
			//std::cout << "Multiplier Perlina: " << std::max(0.1, 1 + (perlin_inc * 1)) << std::endl;
			//add mass perlin noise
			perlin_freq = 5.;
			perlin_inc = perlin.accumulatedOctaveNoise3D(fluid.particles.X(i)[0] * perlin_freq, fluid.particles.X(i)[1] * perlin_freq, fluid.particles.X(i)[2] * perlin_freq, octaves);
			real perlin_multiplier_density = std::max(0.1, 1 + (perlin_inc * 1.));
			//std::cout << "Multiplier Perlina: " << std::max(0.7, 1 + (perlin_inc * 0.2)) << std::endl;
			//add mass perlin noise
			fluid.particles.M(i) *= perlin_multiplier_V * perlin_multiplier_density;
		}
	}

	fluid.vorticity_params = VorticityParams(KernelType::QUINTIC, 1, true, 1e-4);
	if (first_frame < 1) {
		Seed_Vortex_Rand(20, 0.000001, 0.2 * fluid.simulation_scale);
	}
	//Seed_Vortex_Rand(20, 0.000001, 0.2*fluid.simulation_scale);
	fluid.dynamic_seed_vortex = true;
	fluid.dynamic_seed_rate = 1. / 10.;
	fluid.Initialize(dx);
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_11(void)
{
	std::cout << "Enter initialization Case 11\n";

	max_iter_per_frame = -1;
	cfl = 0.1;
	frame_rate = 50;
	real dx = 0.005;
	if constexpr (d == 3) { PointSetFunc::Initialize_Box_Points(10, 10, 10, dx, VectorD::Zero(), fluid_3d.particles); }

	real rho = 1e3;
	fluid.simulation_scale = dx * 10;
	fluid.particles.Add_Element();
	Initialize_Particle(0, 1.0, 1.0);
	fluid.particles.E(0) = MatrixD::Identity();
	fluid.particles.X(0) = AuxFunc::V<d>(0, 1e4, 0);
	Set_Physical_Parameters();

	for (int i = 0; i < fluid_3d.particles.Size(); i++) {
		fluid_3d.particles.M(i) = 3e-8;
		for (int axis = 0; axis < d; axis++) {
			real offset = (Rand_Number() - 0.5) * dx * 0.5;
			fluid_3d.particles.X(i)[axis] += offset;
		}
	}

	//world forces
	fluid.gravity_coeff = 0;// 0.02;

	fluid.max_vel = 0.005 * 50;//move 0.005 in one frame with FPS=50
	bool use_3d = true;
	if (use_3d) {
		fluid_3d.Initialize(dx, rho, fluid.simulation_scale);
		fluid_3d.use_surface_tension = true;
		fluid_3d.nden_0 = 3e6;
		fluid_3d.pressure_multiplier = 5e3;
		fluid_3d.curvature_multiplier = 1.;
		fluid_3d.viscosity_multiplier = 1;
		fluid_3d.analytical_boundary.Add_Obstacle(std::make_shared<Plane<d>>(VectorD::Unit(1), -VectorD::Unit(1) * 6 * dx));
		fluid.Initialize(dx, nullptr, &fluid_3d);
	}

}

//template<int d>
//void FluidSPHBubbleDriver<d>::Case_12(void)
//{
//	if constexpr (d !=3) {
//		std::cout << "case_12 error: d!=3" << std::endl;
//	}
//	if constexpr (d == 3) {
//
//		real vel_term = 0.05, capillary_term = 0.05;
//		std::ifstream fin(scene_file_name);
//		if (!fin.is_open()) { std::cerr << "Driver Case_3 error: cannot open " << scene_file_name << "\n"; exit(0); }
//		fin >> vel_term >> capillary_term;
//		fin.close();
//		//real vel_term = 0.5, capillary_term = 0.5;
//
//		max_iter_per_frame = -1;
//		cfl = 0.1;
//		frame_rate = 50;
//		fluid.max_vel = 0.1;
//
//		real R = (real)0.1;//0.1:1k, 0.3:11k, 1.0:126k
//		fluid.simulation_scale = R;
//		real dx = 0.005;
//		Array<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);
//
//		Set_Physical_Parameters();
//		real rho = 1e3;
//		fluid.thickness = 5e-7;//10um
//		real default_conc = 0;//"extreme soap"
//		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
//		real total_vol = total_surface * fluid.thickness;
//		real total_mass = total_vol * rho;
//
//		VectorD ctr = VectorD::Unit(0) * (-2 * R);
//		MatrixD Rot; Rot << 0, 1, 0,
//			0, 0, 1,
//			1, 0, 0;
//		int n = fluid.particles.Size();
//		for (int i = 0; i < n; i++) {
//			Initialize_Particle(i, total_mass / n, fluid.thickness);
//			fluid.particles.B(i) = is_boundary[i];
//			fluid.particles.Vol(i) = total_vol / n;
//			fluid.particles.Conc(i) = default_conc;
//			fluid.particles.RH(i) = fluid.thickness;
//			fluid.particles.Apply_Rotation(i, Rot);
//			fluid.particles.Apply_Translation(i, ctr);
//		}
//
//		real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
//
//		//world forces
//		fluid.gravity_coeff = 0;
//		//tangential force
//		fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, frame_rate);
//		fluid.t_pressure_params.height_pressure_coeff *= 5;
//		fluid.divergence_params.calculate_field = "all";
//		fluid.viscosity_coeff = 0.4 * alpha_vis;
//		fluid.marangoni_coeff = 0;
//		//render height
//		fluid.rh_params = RenderHeightParams("divergence", true, true, 0.01);
//		//normal
//		fluid.normal_viscosity_coeff = 0;
//		fluid.n_pressure_params.Set_IB(R, fluid.gamma_water, rho, fluid.thickness);
//		fluid.n_pressure_params.capillary_coeff *= capillary_term;
//		fluid.n_pressure_params.ib_force_coeff = 2;
//		EulerInitializer<d> perimeter;
//		real source_speed = vel_term;
//		perimeter.Set_Domain(R, 10, AuxFunc::Vi<d>(6, 3, 3));
//		perimeter.Set_Boundary_Width(1, -1, 1, 1, 1, 1);
//		perimeter.Set_Boundary_Value(0, 0, 0, 0, 0, 0);
//		perimeter.Generate_Parameters();
//		//boundary
//		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
//		fluid.boundary_params.replenish = true;
//		fluid.boundary_params.replenish_interval = 1e-3;
//		fluid.boundary_params.keep_xz_plane = false;
//		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Tube<d>>(ctr, R, VectorD::Unit(0), 0.5 * R));
//		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Plane<d>>(VectorD::Unit(0), ctr - VectorD::Unit(0) * (0.5 * dx)));
//
//		fluid.delete_solitary = false;
//
//		fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
//		fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);
//
//		fluid.Initialize(dx, &perimeter);
//
//		real nozzle_radius = R * 0.5;
//		VectorD nozzle_ctr = AuxFunc::V<d>(-3, 0, 0) * R;
//		Sphere<d> sphere(nozzle_ctr, nozzle_radius);
//		for (auto p : fluid.air_solver.bc.psi_N_values) {
//			int axis = p.first[0];
//			int face_index = p.first[1];
//			VectorDi face = fluid.air_solver.mac_grid.Face_Coord(axis, face_index);
//			VectorD pos = fluid.air_solver.mac_grid.Face_Center(axis, face);
//			if (axis == 0 && sphere.Inside(pos)) {
//				fluid.air_solver.bc.Set_Psi_N(axis, face, source_speed);
//			}
//		}
//		fluid.air_solver.kernel_coeff = 0.9;
//	}
//}

template<int d>
void FluidSPHBubbleDriver<d>::Case_12(void)
{
	if constexpr (d != 3) {
		std::cout << "case_12 error: d!=3" << std::endl;
	}
	if constexpr (d == 3) {

		//real vel_term = 0.05, capillary_term = 0.05;
		//std::ifstream fin(scene_file_name);
		//if (!fin.is_open()) { std::cerr << "Driver Case_3 error: cannot open " << scene_file_name << "\n"; exit(0); }
		//fin >> vel_term >> capillary_term;
		//fin.close();
		real vel_term = 0.5, capillary_term = 0.2;//0.5;

		max_iter_per_frame = -1;
		cfl = 1;
		frame_rate = 50;
		fluid.max_vel = 0.1;

		real R = (real)0.1;//0.1:1k, 0.3:11k, 1.0:126k
		fluid.simulation_scale = R;
		real dx = 0.005;
		Array<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();
		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		VectorD ctr = VectorD::Unit(0) * (-2 * R);
		MatrixD Rot; Rot << 0, 1, 0,
			0, 0, 1,
			1, 0, 0;
		int n = fluid.particles.Size();
		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			fluid.particles.RH(i) = fluid.thickness;
			fluid.particles.Apply_Rotation(i, Rot);
			fluid.particles.Apply_Translation(i, ctr);
		}

		real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;

		//world forces
		fluid.g = 9.8 * VectorD::Unit(0);
		fluid.gravity_coeff = 0.04;
		//tangential force
		fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, frame_rate);
		fluid.t_pressure_params.height_pressure_coeff *= 5;
		fluid.divergence_params.calculate_field = "all";
		fluid.viscosity_coeff = 0.4 * alpha_vis;
		fluid.marangoni_coeff = 0;
		//render height
		fluid.rh_params = RenderHeightParams("divergence", true, true, 0.01);
		//normal
		fluid.normal_viscosity_coeff = 0;
		fluid.n_pressure_params.Set_Circle_Baseline(0.1, fluid.Surface_Tension_Coefficient(0.), 1000, fluid.thickness);
		fluid.n_pressure_params.capillary_coeff *= capillary_term;
		//fluid.n_pressure_params.ib_force_coeff = 2;
		//EulerInitializer<d> perimeter;
		//real source_speed = vel_term;
		//perimeter.Set_Domain(R, 10, AuxFunc::Vi<d>(6, 3, 3));
		//perimeter.Set_Boundary_Width(1, -1, 1, 1, 1, 1);
		//perimeter.Set_Boundary_Value(0, 0, 0, 0, 0, 0);
		//perimeter.Generate_Parameters();
		//boundary
		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
		fluid.boundary_params.replenish = true;
		fluid.boundary_params.replenish_interval = 1e-3;
		fluid.boundary_params.keep_xz_plane = false;
		fluid.boundary_params.replenish_dx_num = 1.5;
		//fluid.analytical_boundary.Add_Obstacle(std::make_shared<Tube<d>>(ctr, R, VectorD::Unit(0), 0.5 * R));
		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Plane<d>>(VectorD::Unit(0), ctr - VectorD::Unit(0) * (0.5 * dx)));

		fluid.delete_solitary = false;

		fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
		fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);

		//fluid.Initialize(dx, &perimeter);
		fluid.Initialize(dx);

		//real nozzle_radius = R * 0.5;
		//VectorD nozzle_ctr = AuxFunc::V<d>(-3, 0, 0) * R;
		//Sphere<d> sphere(nozzle_ctr, nozzle_radius);
		//for (auto p : fluid.air_solver.bc.psi_N_values) {
		//	int axis = p.first[0];
		//	int face_index = p.first[1];
		//	VectorDi face = fluid.air_solver.mac_grid.Face_Coord(axis, face_index);
		//	VectorD pos = fluid.air_solver.mac_grid.Face_Center(axis, face);
		//	if (axis == 0 && sphere.Inside(pos)) {
		//		fluid.air_solver.bc.Set_Psi_N(axis, face, source_speed);
		//	}
		//}
		//fluid.air_solver.kernel_coeff = 0.9;
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_13(void)
{
	std::cout << "Enter initialization Case 13\n";

	int fine_level = 1;
	int sphere_sub_num; real R;
	if (fine_level == 1) {
		sphere_sub_num = 5;
		R = 0.133;
	}
	else if (fine_level == 2) {
		sphere_sub_num = 7;
		R = 0.5320;
	}

	real ext_term = 5.0, vor_term = 5e-6;

	snapshot_stride = 25;
	max_iter_per_frame = -1;
	cfl = 0.1;
	frame_rate = 50;
	real dx = 0.005;
	int number_2d = 1000;//number of particles for 2d
	VectorD ctr = AuxFunc::V<d>(0, R*2.5, 0);
	SPHBubbleParticles<d> full_particles;
	if constexpr (d == 2) {
		R = dx * number_2d / (2 * pi);
		dx = PointSetFunc::Initialize_Circle_Points(ctr, R, number_2d, full_particles);
	}
	if constexpr (d == 3) { dx = PointSetFunc::Initialize_Sphere_Points(ctr, R, sphere_sub_num, full_particles); }
	fluid.simulation_scale = R;
	fluid.particles.Resize(0);
	for (int i = 0; i < full_particles.Size(); i++) {
		if (full_particles.X(i)[1] - ctr[1] < dx * 0.1) {
			int idx = fluid.particles.Add_Element();
			fluid.particles.Copy_Element_From(idx, full_particles, i);
		}
	}
	std::cout << "Initialize with " << fluid.particles.Size() << " particles and dx = " << dx << " and R = " << R << "\n";
	std::cout << "ctr: "<<ctr.transpose() << "\n";


	Set_Physical_Parameters();
	real rho = 1e3; fluid.thickness = 5e-7;//10um
	real default_conc = 1;//"extreme soap"
	real total_surface = pi * 0.5 * pow(2 * R, 2);//works for both 2D and 3D
	real total_vol = total_surface * fluid.thickness, total_mass = total_vol * rho;
	int n = fluid.particles.Size();
	std::cout << "Particle Vol: " << total_vol / n << ", Particle Mass: " << total_mass / n << "\n";
	for (int i = 0; i < n; i++) {
		Initialize_Particle(i, total_mass / n, fluid.thickness);
		fluid.particles.B(i) = 0;
		fluid.particles.Vol(i) = total_vol / n;
		fluid.particles.Conc(i) = default_conc;
		if (fluid.particles.X(i)[1] - ctr[1] >= -dx*0.1) {
			fluid.particles.B(i) = 1;
		}
	}

	real V = (d == 2) ? R * R * pi * 0.5 : 2.0 / 3.0 * pi * R * R * R;
	std::cout << "Analytical volume: " << V << "\n";
	real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
	real alpha_vel = 0.005 * 50;//move 0.005 in one frame with FPS=50

	//world forces
	fluid.g = -9.8 * VectorD::Unit(1);
	fluid.gravity_coeff = 0.1;
	//fluid.gravity_coeff = 0.00088 / R;
	//tangential force
	fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, frame_rate);
	fluid.t_pressure_params.height_pressure_coeff *= 10;

	fluid.divergence_params.calculate_field = "all";
	fluid.viscosity_coeff = 0.4 * alpha_vis;
	fluid.marangoni_coeff = 0;
	//vorticity params
	fluid.vorticity_params = VorticityParams(KernelType::QUINTIC, 0, false, 0);
	//render height
	fluid.rh_params = RenderHeightParams("divergence", false, true, 1e-4);
	//normal
	fluid.normal_viscosity_coeff = 0;
	fluid.n_pressure_params.Set_IB(R, fluid.Surface_Tension_Coefficient(default_conc), rho, fluid.thickness);
	fluid.n_pressure_params.air_pressure_mode = "none";
	EulerInitializer<d> perimeter;
	perimeter.Set_Domain(R, 10, AuxFunc::Vi<d>(3, 6, 3));
	perimeter.Set_Boundary_Width(0, 0, 0, 0, 0, 0);
	perimeter.Set_Boundary_Value(0, 0, 0, 0, 0, 0);
	perimeter.Generate_Parameters();
	//boundary
	fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
	fluid.boundary_params.keep_xz_plane = false;
	fluid.boundary_params.replenish = true;
	fluid.boundary_params.replenish_interval = 1e-3;
	fluid.boundary_params.replenish_proportion = 1.0;
	fluid.boundary_params.replenish_dx_num = 2;
	fluid.boundary_params.replenish_criteria = "particle";
	fluid.delete_solitary = true;
	fluid.exp_mode = "dripping";
	fluid.analytical_boundary.Add_Obstacle(std::make_shared<Plane<d>>(-VectorD::Unit(1), ctr + 0.5 * dx * VectorD::Unit(1)));
	fluid.analytical_boundary.Add_Obstacle(std::make_shared<Tube<d>>(ctr, R + dx * 0.5, VectorD::Unit(1), 0.5 * R));
	

	fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
	fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);

	fluid.max_vel = alpha_vel;

	VectorD floor_ctr = VectorD::Zero() - VectorD::Unit(1) * R * 3;
	real floor_side = R * 3;
	int floor_num = int(floor_side / dx);
	VectorD floor_min = AuxFunc::V<d>(floor_ctr[0] - floor_side * 0.5, floor_ctr[1] + 0.5 * dx, floor_ctr[2] - floor_side * 0.5);
	/*if constexpr (d == 3) {
		std::vector<int> is_boundary = PointSetFunc::Initialize_Lattice_Points(floor_ctr, Vector2i(floor_num, floor_num), VectorD::Unit(0), VectorD::Unit(2), dx, fluid_3d.particles);
	}
	for (int i = 0; i < fluid_3d.particles.Size(); i++) fluid_3d.particles.B(i) = 1;*/
	fluid_3d.Initialize(dx, rho, R, fluid.g * fluid.gravity_coeff);
	fluid_3d.use_central_gravity = false;
	fluid_3d.use_surface_tension = true;
	fluid_3d.prune_far_away = false;
	fluid_3d.low_bound = -3 * R;
	fluid_3d.nden_0 = 5e6;
	fluid_3d.pressure_multiplier = 5e2;
	fluid_3d.curvature_multiplier = 0.1;
	fluid_3d.cohesion_multiplier = 1;
	fluid_3d.viscosity_multiplier = 0.5;
	fluid_3d.prune_low = true;
	//fluid_3d.analytical_boundary.Add_Obstacle(std::make_shared<Plane<d>>(VectorD::Unit(1), floor_ctr));
	
	

	//fluid.Initialize(dx, &perimeter, &fluid_3d);
	fluid.Initialize(dx, nullptr, &fluid_3d);

	for (int i = 0; i < fluid.particles.Size(); i++) fluid.particles.RH(i) = fluid.particles.H(i);
	/*real V_num = fluid.Compute_Enclosed_Volume();
	fluid.n_pressure_params.Set_Sphere_Weak(d, R, fluid.Surface_Tension_Coefficient(default_conc), V_num, rho, fluid.thickness);
	fluid.n_pressure_params.pressure_constant_0 = fluid.n_pressure_params.pressure_constant;
	std::cout << "pressure constant: " << fluid.n_pressure_params.pressure_constant << "\n";*/
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_14(void)
{
	if constexpr (d == 2) {
		std::string numerical_test = "circle";
		if (numerical_test == "circle") {
			real R = 1, dx;
			int number_2d = 2000;//number of particles for 2d
			int sphere_sub_num = 5;//7: 163842
			VectorD ctr = VectorD::Zero();
			if constexpr (d == 2) { dx = PointSetFunc::Initialize_Circle_Points(ctr, R, number_2d, fluid.particles); }
			if constexpr (d == 3) { dx = PointSetFunc::Initialize_Sphere_Points(ctr, R, sphere_sub_num, fluid.particles); }
			std::cout << "Initialize with " << fluid.particles.Size() << " particles\n";
			Set_Physical_Parameters();
			real rho = 1e3; fluid.thickness = 5e-7;//10um
			real default_conc = 1;//"extreme soap"
			real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
			real total_vol = total_surface * fluid.thickness, total_mass = total_vol * rho;
			int n = fluid.particles.Size();
			for (int i = 0; i < n; i++) {
				Initialize_Particle(i, total_mass / n, 5e-7);
				fluid.particles.B(i) = 0;
				fluid.particles.Vol(i) = total_vol / n;
				fluid.particles.Conc(i) = default_conc;
			}
			fluid.boundary_params.particle_geometry_mode = "mirror_compensate";//mirror compensate for d=3
			fluid.Initialize(dx);

			for (int i = 0; i < n; i++) {
				if (10 <= i && i <= 20) fluid.particles.V(i) = fluid.particles.X(i) * (-1);// *cos((i - 15.0) / 10.0 * pi / 2);
			}

			/*for (int i = 0; i <= 30; i++) {
				auto nominal_vel_i = [&](const int idx)->VectorD {return fluid.Nominal_Vector(fluid.particles.X(i), fluid.particles.Normal(i), fluid.particles.X(idx), fluid.particles.Normal(idx), fluid.particles.V(idx));  };
				VectorD vis_force = fluid.Surface_Laplacian<VectorD>(i, nominal_vel_i, fluid.viscosity_params);
				VectorD norm = fluid.particles.Normal(i);
				VectorD vis_normal = norm * vis_force.dot(norm) / norm.squaredNorm();//normal component
				VectorD vis_tangential = vis_force - vis_normal;//tangential component
				std::cout << "particle " << i << " vis normal: " << vis_normal.transpose() << " vis_tan: " << vis_tangential.transpose() << "\n";
				//particles.F(i) += vis_normal * normal_viscosity_coeff + vis_tangential * viscosity_coeff;
			}*/

		}
		else if (numerical_test == "kernel") {
			real length = 1, total_mass = 1, total_vol = 1;
			int n = 200;
			real h = length / n;
			VectorD normal = VectorD::Unit(1);
			fluid.Initialize(h);
			std::cout << "initialize with h=" << h << " kernel h=" << (fluid.kernel)->h << std::endl;
			real kernel_sum = 0;
			int rs = 6;
			for (int i = -rs; i <= rs; i++) {
				kernel_sum += fluid.Surface_Kernel_Weight(fabs(i * h), KernelType::QUINTIC) * h;
			}
			std::cout << "kernel discrete integration: " << kernel_sum << std::endl;
		}
		else if (numerical_test == "segment") {
			real length = 1, total_mass = 1, total_vol = 1;
			int n = 100 + 1;
			real h = length / n;
			VectorD normal = VectorD::Unit(1);
			VectorD start = AuxFunc::V<d>(-0.5 * length), end = AuxFunc::V<d>(0.5 * length);
			PointSetFunc::Initialize_Segment_Points(start, end, n, fluid.particles);
			for (int i = 0; i < fluid.particles.Size(); i++) {
				Initialize_Particle(i, total_mass / n, 5e-7);
				fluid.particles.Vol(i) = total_vol / n;
			}
			for (int i = 0; i < fluid.particles.Size(); i++) {
				normal = normal.normalized();
				VectorD t1 = VectorD::Unit(0);
				fluid.particles.E(i).col(0) = t1;
				fluid.particles.E(i).col(1) = normal;
			}
			fluid.boundary_params.particle_geometry_mode = "mirror_compensate";
			fluid.Initialize(h);
			/*for (int i = 0; i < n; i++) {
				if (i < n / 2) fluid.particles.Conc(i) = 0;
				else fluid.particles.Conc(i) = fluid.particles.X(i)[0];
			}*/
			//for (int i = 0; i < n; i++) fluid.particles.Conc(i) = pow(fluid.particles.X(i)[0], 2) * 0.5;//its laplacian should be 1
			for (int i = 0; i < n; i++) fluid.particles.Conc(i) = pow(fluid.particles.X(i)[0],1);
			//std::cout << "pos:\n";for (int i = 0; i < n; i++) std::cout << fluid.particles.X(i).transpose() << " , "; std::cout << "\n";
			//std::cout << "Vol:\n"; for (int i = 0; i < n; i++) std::cout << fluid.particles.Vol(i) << " , "; std::cout << "\n";
			std::cout << "SA:\n"; for (int i = 0; i < n; i++) std::cout << fluid.particles.SA(i) << " , "; std::cout << "\n";
			std::cout << "H:\n"; for (int i = 0; i < n; i++) std::cout << fluid.particles.H(i) << " , "; std::cout << "\n";
			std::cout << "symm grad of conc:\n";
			for (int i = 0; i < n; i++) {
				VectorT grad = fluid.Surface_Gradient_Symmetric(i, Index_Function(fluid.particles.ConcRef()), OperatorParams("all", KernelType::QUINTIC));
				std::cout << i << " : " << grad.transpose() << std::endl;
			}
			std::cout << "diff grad of conc:\n";
			for (int i = 0; i < n; i++) {
				VectorT grad = fluid.Surface_Gradient_Difference(i, Index_Function(fluid.particles.ConcRef()), OperatorParams("all", KernelType::QUINTIC));
				std::cout << i << " : " << grad.transpose() << std::endl;
			}
			/*std::cout << "lap of conc:\n";
			for (int i = 0; i < n; i++) {
				real lap = fluid.Surface_Laplacian(i, Index_Function(fluid.particles.ConcRef()), OperatorParams("all", KernelType::QUINTIC));
				std::cout << i << " conc: " << fluid.particles.Conc(i) << " lap: " << lap << std::endl;
			}*/
		}
	}
	else {
		std::cerr << "initialization error: dimension wrong\n";
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_15(void) ////2D Disk in 3D but with grid particles
{
	if constexpr (d == 3) {
		real R = (real)1;

		real dx = 0.06; //this can be changed

		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		int num_fluid_particles = 0;
		//for (int i = 0; i < fluid.particles.Size(); i++) {
		//	fluid.particles.B(i) = is_boundary[i];
		//	if (fluid.particles.B(i) != 1) {
		//		num_fluid_particles++;
		//	}
		//}
		Set_Physical_Parameters();
		real rho = 1e3, thickness = 1e-7;
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * thickness, total_mass = total_vol * rho;
		int n = fluid.particles.Size();
		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
		}
		//fluid.particle_geometry_mode = (d == 2) ? "geometric_compensate" : "mirror_compensate";//mirror compensate for d=3
		fluid.Initialize(dx);
		//real total_sa = 0;
		//for (int i = 0; i < fluid.particles.Size(); i++) {
		//	total_sa += fluid.particles.SA(i);
		//}

		//std::cout << "total surface area: " << total_sa << std::endl;
		for (int i = 0; i < fluid.particles.Size();i++) {
			if (fluid.particles.X(i).norm() < 0.7) {
				fluid.particles.V(i) = -1 * fluid.particles.X(i);
				//std::cout << i << std::endl;
				//std::cout << particles.V(i) << std::endl;
			}
		}

		/*fluid.gravity_coeff = 0;
		fluid.tangential_pressure_coeff = 3000000;
		fluid.viscosity_coeff = 30000;
		fluid.marangoni_coeff = 0;
		fluid.capillary_coeff = 0;*/

		std::cout << "particles number=" << fluid.particles.Size() << std::endl;
		std::cout << "dx=" << dx << std::endl;
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_16(void) // test what's up with the kernels
{
	if constexpr (d == 3) {
		//real R = (real)1;
		//real dx = R * pi * 2 / 100;;
		////fluid.Initialize(dx);
		//fluid.surface.t_r = R * 2;

		////real dx=R*two_pi/(real)scale;
		//Vector3 normal = Vector3(0, -1, 0);
		//PointSetFunc::Initialize_Circle_Points(VectorD::Zero(), R, normal, dx, fluid.particles);
		//for (int i = 0; i < fluid.particles.Size(); i++) {
		//	Initialize_Particle(i);
		//}
		real R = (real)1; real dx = 0.06;
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		for (int i = 0; i < fluid.particles.Size(); i++) {
			Initialize_Particle(i, default_mass,5e-7);
		}
		for (int i = 0; i < fluid.particles.Size(); i++) {
			fluid.particles.B(i) = is_boundary[i];
		}

		//fluid.gravity_coeff = 0;
		//fluid.boundary_force_mode = "none";
		//fluid.capillary_coeff *= 0.1;
		//fluid.capillary_coeff *= 0;

		//fluid.tangential_pressure_coeff = 5e-2;
		//fluid.viscosity_coeff *= 0.5;


		fluid.Initialize(dx);
		std::cout << "particles number=" << fluid.particles.Size() << std::endl;
		std::cout << "dx=" << dx << std::endl;

		//fluid.Initialize(dx);
		for (int i = 0; i < fluid.particles.Size(); i++) {
			if (i != 814) {
				continue;
			}
			VectorT temp;
			temp << fluid.particles.X(i)[0], fluid.particles.X(i)[2];
			real value = fluid.kernel->W_Poly6(temp.norm());
			//real value = fluid.kernel->W_Spiky(temp.norm());
			//real value = fluid.kernel->W_Cubic(temp.norm());
			//real value = fluid.kernel->W_Intp4(temp.norm());
			//real value = fluid.kernel->W_Bell(temp.norm());
			//real value = fluid.kernel->W_Gaussian(temp.norm());
			//real value = fluid.kernel->W_New_Quadratic(temp.norm());
			fluid.particles.H(i) = 60*value;
			VectorD temp2;
			VectorT temp1 = fluid.kernel->Grad_Poly6(-temp);
			/*VectorT temp1 = fluid.kernel->Grad_Spiky(-temp);
			VectorT temp1 = fluid.kernel->Grad_Cubic(-temp);
			VectorT temp1 = fluid.kernel->Grad_Intp4(-temp);
			VectorT temp1 = fluid.kernel->Grad_Bell(-temp);
			VectorT temp1 = fluid.kernel->Grad_Gaussian(-temp);
			VectorT temp1 = fluid.kernel->Grad_New_Quadratic(-temp);*/
			temp2 << temp1[0], 0, temp1[2];
			if (i == 814) {
				std::cout << -temp << std::endl;
				std::cout << temp2 << std::endl;
			}
			fluid.particles.V(i) = temp2;
		}
		//std::cout << "particles number=" << fluid.particles.Size() << std::endl;
		//std::cout << "dx=" << dx << std::endl;
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_17(void) //unified 3d membrane with circle rim
{
	if constexpr (d == 3) {
		real R = (real)1;

		real dx = 0.06; //this can be changed

		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		int num_fluid_particles = 0;
		//for (int i = 0; i < fluid.particles.Size(); i++) {
		//	fluid.particles.B(i) = is_boundary[i];
		//	if (fluid.particles.B(i) != 1) {
		//		num_fluid_particles++;
		//	}
		//}
		Set_Physical_Parameters();
		real rho = 1e3, thickness = 1e-5;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_mass = total_surface * rho, total_vol = total_surface * thickness;
		int n = fluid.particles.Size();
		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n,5e-7);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
		}
		//fluid.particle_geometry_mode = (d == 2) ? "geometric_compensate" : "mirror_compensate";//mirror compensate for d=3
		fluid.Initialize(dx);
		//real total_sa = 0;
		//for (int i = 0; i < fluid.particles.Size(); i++) {
		//	total_sa += fluid.particles.SA(i);
		//}

		//std::cout << "total surface area: " << total_sa << std::endl;


		//fluid.gravity_coeff = 0;
		//fluid.tangential_pressure_coeff = 1;
		fluid.viscosity_coeff = 1;
		fluid.marangoni_coeff = 0;
		//fluid.closed = true;
		//fluid.pressure_constant = V * inner_p;
		//fluid.air_pressure_mode = "pv";
		//fluid.pv_air_force_coeff = 1;
		//fluid.capillary_coeff = 1000000;

		std::cout << "particles number=" << fluid.particles.Size() << std::endl;
		std::cout << "dx=" << dx << std::endl;
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_18(void) {//half sphere
	std::cout << "Enter initialization Case 18\n";

	int fine_level = 2;
	int sphere_sub_num; real R;
	if (fine_level == 1) {
		sphere_sub_num = 5;
		R = 0.133;
	}
	else if (fine_level == 2) {
		sphere_sub_num = 7;
		R = 0.5320;
	}
	
	real ext_term = 5.0, vor_term = 5e-6;

	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 3) {
		const static real ext_coeff = ext_term;
		ext_force = [&](const int i)->VectorD {
			const VectorD& pos = fluid.particles.X(i);
			const real r = fluid.simulation_scale * 0.25;
			VectorD axis = AuxFunc::V<d>(1, 2, 0).normalized();
			real sqr_tan_len = AuxFunc::Eliminate_Unit_Component(pos, axis).squaredNorm();
			VectorD force = VectorD::Zero();

			if (sqr_tan_len <= r * r && current_frame < 50) {
				VectorD ext_norm = -fluid.particles.X(i) * fluid.particles.M(i);
				force += ext_norm * Center_Smooth_Kernel(sqrt(sqr_tan_len), r) * ext_coeff;
			}
			return force;
		};
	}

	max_iter_per_frame = -1;
	cfl = 0.1;
	frame_rate = 50;
	real dx = 0.005;
	int number_2d = 1000;//number of particles for 2d
	VectorD ctr = VectorD::Zero();
	SPHBubbleParticles<d> full_particles;
	if constexpr (d == 2) {
		R = dx * number_2d / (2 * pi);
		dx = PointSetFunc::Initialize_Circle_Points(ctr, R, number_2d, full_particles);
	}
	if constexpr (d == 3) { dx = PointSetFunc::Initialize_Sphere_Points(ctr, R, sphere_sub_num, full_particles); }
	fluid.simulation_scale = R;
	fluid.particles.Resize(0);
	for (int i = 0; i < full_particles.Size(); i++) {
		if (full_particles.X(i)[1] >= -dx*5) {
			int idx = fluid.particles.Add_Element();
			fluid.particles.Copy_Element_From(idx, full_particles, i);
		}
	}
	std::cout << "Initialize with " << fluid.particles.Size() << " particles and dx = " << dx << " and R = " << R << "\n";


	Set_Physical_Parameters();
	real rho = 1e3; fluid.thickness = 5e-7;//10um
	real default_conc = 1;//"extreme soap"
	real total_surface = pi * 0.5 * pow(2 * R, 2);//works for both 2D and 3D
	real total_vol = total_surface * fluid.thickness, total_mass = total_vol * rho;
	int n = fluid.particles.Size();
	std::cout << "Particle Vol: " << total_vol / n << ", Particle Mass: " << total_mass / n << "\n";
	for (int i = 0; i < n; i++) {
		Initialize_Particle(i, total_mass / n, fluid.thickness);
		fluid.particles.B(i) = 0;
		fluid.particles.Vol(i) = total_vol / n;
		fluid.particles.Conc(i) = default_conc;
		if (fluid.particles.X(i)[1] <= dx) {
			fluid.particles.B(i) = 1;
		}
	}

	real V = (d == 2) ? R * R * pi * 0.5 : 2.0 / 3.0 * pi * R * R * R;
	std::cout << "Analytical volume: " << V << "\n";
	real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
	real alpha_vel = 0.005 * 50;//move 0.005 in one frame with FPS=50

	//world forces
	fluid.g = -9.8 * VectorD::Unit(1);
	fluid.gravity_coeff = 0.00088 / R;
	fluid.external_force_func = ext_force;
	//tangential force
	fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, frame_rate);
	fluid.t_pressure_params.height_pressure_coeff *= 10;

	fluid.divergence_params.calculate_field = "all";
	fluid.viscosity_coeff = 0.4 * alpha_vis;
	fluid.marangoni_coeff = 0;
	//vorticity params
	fluid.vorticity_params = VorticityParams(KernelType::QUINTIC, 1, true, 1e-4);
	Seed_Vortex(500, vor_term, dx * 3);
	//render height
	fluid.rh_params = RenderHeightParams("divergence", false, true, 1e-4);
	//normal viscosity
	fluid.normal_viscosity_coeff = 0;
	//boundary
	fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
	fluid.boundary_params.keep_xz_plane = false;
	fluid.boundary_params.replenish = true;
	fluid.boundary_params.replenish_interval = 1e-3;
	fluid.boundary_params.replenish_proportion = 1.0;
	fluid.analytical_boundary.Add_Obstacle(std::make_shared<Plane<d>>(VectorD::Unit(1), -1 * dx * VectorD::Unit(1)));
	fluid.delete_solitary = false;

	fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
	fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);

	fluid.max_vel = alpha_vel;
	fluid.Initialize(dx);

	for (int i = 0; i < fluid.particles.Size(); i++) fluid.particles.RH(i) = fluid.particles.H(i);
	real V_num = fluid.Compute_Enclosed_Volume();
	fluid.n_pressure_params.Set_Sphere_Baseline(d, R, fluid.Surface_Tension_Coefficient(default_conc), V_num, rho, fluid.thickness);
	std::cout << "pressure constant: " << fluid.n_pressure_params.pressure_constant << "\n";
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_19(void) {//test marangoni effect
	if constexpr (d == 3) {
		real R = (real)1;

		real dx = 0.06; //this can be changed

		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		int num_fluid_particles = 0;
		Set_Physical_Parameters();
		real rho = 1e3, thickness = 1e-5;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * thickness, total_mass = total_vol * rho;
		int n = fluid.particles.Size();
		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n,thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
		}
		//fluid.particle_geometry_mode = (d == 2) ? "geometric_compensate" : "mirror_compensate";//mirror compensate for d=3
		fluid.Initialize(dx);


		for (int i = 0; i < fluid.particles.Size();i++) {
			if (fluid.particles.X(i).norm() < 0.5) {
				fluid.particles.Conc(i) = 0.2;
			}
		}

		//fluid.gravity_coeff = 0;
		//fluid.thick_coeff = 0;
		//fluid.tangential_pressure_coeff = 30000;
		//fluid.viscosity_coeff = 3;
		//fluid.marangoni_coeff = 1;
		//fluid.capillary_coeff = 0;

		std::cout << "particles number=" << fluid.particles.Size() << std::endl;
		std::cout << "dx=" << dx << std::endl;
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_20(void) {//test planar effect
	cfl = 0.1;
	if constexpr (d == 3) {
		real R = (real)1;

		real dx = 0.06; //this can be changed

		VectorD ctr = VectorD::Zero();
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(ctr, R, dx, fluid.particles);
		int num_fluid_particles = 0;
		Set_Physical_Parameters();
		real rho = 1e3, thickness = 1e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * thickness, total_mass = total_vol * rho;
		int n = fluid.particles.Size();
		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n,thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
		}
		//fluid.particle_geometry_mode = (d == 2) ? "geometric_compensate" : "mirror_compensate";//mirror compensate for d=3
		fluid.Initialize(dx);


		for (int i = 0; i < fluid.particles.Size();i++) {
			fluid.particles.Conc(i) = 0.3 * Rand_Number();
		}

		//fluid.gravity_coeff = 0;
		//fluid.thick_coeff = 0;
		//fluid.tangential_pressure_coeff = 1000000;
		//fluid.viscosity_coeff = 300000;
		//fluid.marangoni_coeff = 0.02;
		//fluid.capillary_coeff = 0;
		//fluid.grav_boundary = 0;
		//fluid.vis_boundary = 0.01 * fluid.viscosity_coeff;
		//fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(ctr, R));
		//fluid.boundary_mode = "rotating_boundary";

		std::cout << "particles number=" << fluid.particles.Size() << std::endl;
		std::cout << "dx=" << dx << std::endl;
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_21(void)
{
	if(d!=3){
		std::cerr << "Case 21 must d==3\n"; exit(0);
	}

	std::cout << "Enter initialization Case 21\n";

	max_iter_per_frame = -1;
	cfl = 0.1;
	frame_rate = 50;
	real dx = 0.005;
	int nh = 20, nr = 100;
	VectorD ctr = VectorD::Zero();
	real R = dx * nr, height = dx * (nh - 1.0);
	Array<int> bnd_flags;
	if constexpr (d == 3) {
		bnd_flags = PointSetFunc::Initialize_Catenoid_Points(ctr, R, nr, height, nh, fluid.particles);
	}
	std::cout << "Initialize with " << fluid.particles.Size() << " particles and dx = " << dx << " and R = " << R << "\n";


	Set_Physical_Parameters();
	real rho = 1e3; fluid.thickness = 5e-7;//10um
	real default_conc = 1;//"extreme soap"
	real total_surface = pi * 0.5 * pow(2 * R, 2);//works for both 2D and 3D
	real total_vol = total_surface * fluid.thickness, total_mass = total_vol * rho;
	int n = fluid.particles.Size();
	std::cout << "Particle Vol: " << total_vol / n << ", Particle Mass: " << total_mass / n << "\n";
	for (int i = 0; i < n; i++) {
		Initialize_Particle(i, total_mass / n, fluid.thickness);
		fluid.particles.B(i) = bnd_flags[i];
		fluid.particles.Vol(i) = total_vol / n;
		fluid.particles.Conc(i) = default_conc;
	}

	real V = (d == 2) ? R * R * pi * 0.5 : 2.0 / 3.0 * pi * R * R * R;
	std::cout << "Analytical volume: " << V << "\n";
	real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
	real alpha_vel = 0.005 * 50;//move 0.005 in one frame with FPS=50

	//world forces
	fluid.g = -9.8 * VectorD::Unit(1);
	fluid.gravity_coeff = 0;
	//tangential force
	fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, frame_rate);
	fluid.t_pressure_params.height_pressure_coeff *= 10;

	fluid.divergence_params.calculate_field = "all";
	fluid.viscosity_coeff = 0.4 * alpha_vis;
	fluid.marangoni_coeff = 0;
	//vorticity params
	fluid.vorticity_params = VorticityParams(KernelType::QUINTIC, 0, false, 0);
	//render height
	fluid.rh_params = RenderHeightParams("divergence", true, false, 0);
	//normal viscosity
	fluid.normal_viscosity_coeff = 0;
	//boundary
	fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
	//fluid.boundary_params.Set_No_Boundary();
	fluid.boundary_params.boundary_mode = "none";
	//fluid.analytical_boundary.Add_Obstacle(std::make_shared<Plane<d>>(VectorD::Unit(1), -1 * dx * VectorD::Unit(1)));

	fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
	fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);

	fluid.max_vel = alpha_vel;
	fluid.Initialize(dx);

	for (int i = 0; i < fluid.particles.Size(); i++) fluid.particles.RH(i) = fluid.particles.H(i);
	real V_num = fluid.Compute_Enclosed_Volume();
	std::cout << "Numerical Volume: " << V_num << "\n";
	fluid.n_pressure_params.Set_Sphere_Baseline(d, R, fluid.Surface_Tension_Coefficient(default_conc), V_num, rho, fluid.thickness);
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_22(void) //unified 3d membrane with circle rim
{
	if constexpr (d == 3) {
		real R = (real)1;

		real dx = 0.04; //this can be changed

		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//int num_fluid_particles = 0;

		//Set_Physical_Parameters();
		real rho = 1e3;
		real thickness = 1e-7;//10um

		int n = fluid.particles.Size();
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * thickness, total_mass = total_vol * rho;
		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n,thickness);
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = 1.0;
		}

		fluid.Initialize(dx);

		for (int i = 0; i < n; i++) {
			fluid.particles.H(i) = (3*R - fluid.particles.X(i)[0] - fluid.particles.X(i)[2]) * 2 * thickness;
		}
		

		//for (int i = 0; i < n; i++) {
		//	fluid.particles.H(i) = (1.5 - sin(pi*fluid.particles.X(i)[0])) * 4 * thickness;
		//}

		////center sin
		//for (int i = 0; i < n; i++) {
		//	fluid.particles.H(i) = (1.5 + sin(2*pi * fluid.particles.X(i).norm())) * 4 * thickness;
		//}


		//fluid.gravity_coeff = 0;
		//fluid.tangential_pressure_coeff = 0;
		//fluid.viscosity_coeff = 0;
		//fluid.marangoni_coeff = 0;
		//fluid.capillary_coeff = 0;

		std::cout << "particles number=" << fluid.particles.Size() << std::endl;
		std::cout << "dx=" << dx << std::endl;
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_23(void)
{
	if constexpr (d == 3) {
		cfl = 0.1;
		frame_rate = 50;
		real length = 1.0; int ny = 32;
		Vector2i cell_counts = AuxFunc::Vi<2>(2, 1);
		VectorD domain_len = AuxFunc::V<d>(2, 1) * length;
		real dx = length * 1 / ny;
		std::vector<int> is_boundary = PointSetFunc::Initialize_Lattice_Points2(-VectorD::Unit(0) - .5 * VectorD::Unit(2), cell_counts * ny, VectorD::Unit(0), VectorD::Unit(1), dx, fluid.particles);
		Set_Physical_Parameters();
		real rho = 1e3, thickness = 2e-4;
		real total_vol = length * (2 * 1) * thickness, total_mass = total_vol * rho;
		int n = fluid.particles.Size();
		std::cout << "particle number " << n << "\n";
		real y0 = domain_len[1] * (-0.5), y1 = domain_len[1] * 0.5;
		real target_vel = 1e5;
		real total_force = (2 * thickness) * fluid.viscosity_water * target_vel / 1.0;
		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n,thickness);
			//real y = fluid.particles.X(i)[1];
			//if (y < y0 + dx / 2) {
			//	fluid.particles.B(i) = 1;
			//	fluid.particles.F(i) = VectorD::Zero();
			//}
			//else fluid.particles.B(i) = 0;
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = 1.0;
		}
		//fluid.replenish_boundary = true;
		//fluid.particle_geometry_mode = "no_compensate";
		fluid.Initialize(dx);
		//fluid.gravity_coeff = 0;
		//fluid.tangential_pressure_coeff = 450000;
		//fluid.thick_coeff = 0;
		//fluid.viscosity_coeff = 30000;
		//fluid.viscosity_coeff = 0;
		//fluid.marangoni_coeff = 0;
		fluid.n_pressure_params.Set_Planar();
		//fluid.grav_boundary = 0;
		//fluid.boundary_force_mode = "binary";
		//fluid.ignore_boundary_thick = true;
		real all_sa = 0; for (int i = 0; i < n; i++) all_sa += fluid.particles.SA(i);
		std::cout << "all sa: " << all_sa << "\n";

		for (int i = 0; i < n; i++) {
			if (fluid.particles.Is_Boundary(i) == 1)continue;
			if (fluid.particles.X(i)[0] < (-1. + 15 * dx) && !fluid.particles.Is_Boundary(i)) {
				fluid.particles.V(i) += 2 * VectorD::Unit(0);
			}
		}

		////left&right bounding planes
		//VectorD left_end = AuxFunc::V<d>(-1.2, 0) * length, right_end = -left_end;
		//fluid.analytical_boundary.Add_Obstacle(std::make_shared<Plane<d>>(VectorD::Unit(0), left_end));
		//fluid.analytical_boundary.Add_Obstacle(std::make_shared<Plane<d>>(-VectorD::Unit(0), right_end));
		////up&down bounding planes
		//VectorD up_end = AuxFunc::V<d>(0, 0.6) * length, down_end = -up_end;
		//fluid.analytical_boundary.Add_Obstacle(std::make_shared<Plane<d>>(-VectorD::Unit(1), up_end));
		//fluid.analytical_boundary.Add_Obstacle(std::make_shared<Plane<d>>(VectorD::Unit(1), down_end));
		VectorD ctr = 0.5 * VectorD::Unit(0);
		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(ctr, 0.3));

		//fluid.boundary_mode = "slippery";
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_24(void)
{
	if constexpr (d == 2) {
		cfl = 0.01;
		frame_rate = 50;
		int number_2d = 30;//number of particles for 2d
		int number_bnd_2d = (int)(number_2d / 10);
		real R = 1.5, dx = 2 * R / number_2d;
		VectorD start = AuxFunc::V<d>(-R, 0, 0), end = -start;

		if constexpr (d == 2) { PointSetFunc::Initialize_Segment_Points(start, end, number_2d, fluid.particles); }

		Set_Physical_Parameters();
		real rho = 1e3, thickness = 1.0;//e - 7;
			//5e-7;//500nm
		real default_conc = 0;//"extreme soap"
		real total_surface = (d == 2) ? 2 * R : pi * R * R;
		real total_vol = total_surface * thickness, total_mass = total_vol * rho;
		int n = fluid.particles.Size();
		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, thickness);
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			if constexpr (d == 2) {
				int bnd_len = number_bnd_2d;//6
				if (i < bnd_len || i >= n - bnd_len) {
					fluid.particles.B(i) = 1;
				}
				else fluid.particles.B(i) = 0;
			}
		}
		fluid.geometry_params = OperatorParams("all", KernelType::POLY6);
		/*for (int i = 0; i < n; i++) {
			std::cout << "SAi" << fluid.particles.SA(i) << std::endl;
		}

		std::cout << "H:\n"; for (int i = 0; i < n; i++) std::cout << fluid.particles.H(i) << " , "; std::cout << "\n";*/

		fluid.g = VectorD::Unit(0) * 9.8;
		fluid.gravity_coeff = 1.0/9.8;
		fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
		fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);

		//fluid.viscosity_params.calculate_field = "all";
		fluid.marangoni_coeff = 0;
		fluid.n_pressure_params.Set_Planar();

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R- (number_bnd_2d-1.) *dx));

		real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
		real alpha_height = 1.0 * rho / thickness * 100;
		std::cout << "alpha vis: " << alpha_vis << "\n";
		std::cout << "alpha height: " << alpha_height << "\n";

		fluid.viscosity_coeff = alpha_vis;
		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
		fluid.boundary_params.replenish = true;
		//fluid.vis_boundary = 100;
		//fluid.viscosity_coeff = 400;
		fluid.t_pressure_params.Set_Zero();
		fluid.t_pressure_params.height_pressure_coeff = alpha_height;
		fluid.t_pressure_params.tangential_pressure_coeff = 1;

		fluid.Initialize(dx);
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_25(void) //unified 3d membrane with circle rim
{
	if constexpr (d == 3) {
		real R = (real)1;

		real dx = 0.06; //this can be changed

		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		int num_fluid_particles = 0;

		Set_Physical_Parameters();
		real rho = 1e3;
		fluid.thickness = 1e-5;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_mass = total_surface * rho, total_vol = total_surface * fluid.thickness;
		int n = fluid.particles.Size();
		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
		}
		fluid.Initialize(dx);


		fluid.g = VectorD::Unit(0) * 9.8;
		fluid.gravity_coeff = 1.0 / 9.8;
		fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
		fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);

		//fluid.viscosity_params.calculate_field = "all";
		fluid.marangoni_coeff = 0;
		fluid.n_pressure_params.Set_Planar();

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R));

		real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
		real alpha_height = 1.0 * rho;
		std::cout << "alpha vis: " << alpha_vis << "\n";
		std::cout << "alpha height: " << alpha_height << "\n";

		fluid.viscosity_coeff = 100000 * alpha_vis;

		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
		fluid.boundary_params.replenish = true;		

		/*fluid.laplacian_pressure_coeff = 1000000;
		fluid.divergence_pressure_coeff = 200;
		fluid.height_pressure_coeff = 0.5 * alpha_height;
		fluid.tangential_pressure_coeff = 4000000;*/
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_26(void) //unified 3d membrane with circle rim
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
	}
	else if constexpr (d == 3) {
		ext_force = [&](const int i)->VectorD {
			VectorD pos;
			real r;
			VectorD result = VectorD::Zero();
			if (fluid.current_frame > 20) return result;
			pos = fluid.particles.X(i) - (0.4 * VectorD::Unit(0) + 0.4 * VectorD::Unit(2));
			r = 0.3;
			if (pos.norm() < r) {
				result += (r - pos.norm()) * 10 * VectorD::Unit(0);
			}

			pos = fluid.particles.X(i) - (-0.3 * VectorD::Unit(0) - 0.3 * VectorD::Unit(2));
			r = 0.2;
			if (pos.norm() < r) {
				result += (r - pos.norm()) * 10 * -VectorD::Unit(2);
			}

			//return result * 1e-5;
			return result * fluid.particles.M(i);
		};
	}
	if constexpr (d == 3) {
		real R = (real)1;

		real dx = 0.005; //this can be changed

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);
		int num_fluid_particles = 0;

		Set_Physical_Parameters();
		
		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << "\n";
		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
		}

		real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
		real alpha_height = 1.0 * rho;
		real alpha_height_laplacian = 1 * rho / (2 * fluid.gamma_water) / (fluid.thickness / (pow((fluid.thickness + 1), 2. / 3.)));
		real alpha_divergence = rho * frame_rate;
		std::cout << "alpha vis: " << alpha_vis << "\n";
		std::cout << "alpha height: " << alpha_height << "\n";
		std::cout << "alpha height laplacian: " << alpha_height << "\n";
		std::cout << "alpha divergence: " << alpha_height << "\n";

		//world forces
		fluid.gravity_coeff = 0;
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;
		fluid.external_force_func = ext_force;

		//tangential force
		//fluid.laplacian_pressure_coeff = 1 * alpha_height_laplacian;
		//fluid.divergence_pressure_coeff = 0.01 * alpha_divergence;
		//fluid.height_pressure_coeff = 0 * alpha_height;
		//fluid.tangential_pressure_coeff = 1;
		//fluid.viscosity_coeff = 0.3 * alpha_vis;

		//Before changing to alpha_height laplacian and etc.
		fluid.t_pressure_params.laplacian_pressure_coeff = 1 * 15000;
		fluid.t_pressure_params.divergence_pressure_coeff = 1 * 2;
		fluid.t_pressure_params.height_pressure_coeff = 1;
		fluid.t_pressure_params.tangential_pressure_coeff = 10;
		fluid.viscosity_coeff = 2 * alpha_vis;

		//fluid.laplacian_pressure_coeff = 15 * alpha_height_laplacian;
		//fluid.divergence_pressure_coeff = 0.02 * alpha_divergence;
		//fluid.height_pressure_coeff = 1;
		//fluid.tangential_pressure_coeff = 40;
		//fluid.viscosity_coeff = 0.3 * alpha_vis;

		fluid.marangoni_coeff = 0;

		//normal force
		fluid.n_pressure_params.Set_Planar();

		//operator parameters
		fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
		fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::QUINTIC);

		
		fluid.boundary_params.Set_Pure_Analytical();

		fluid.Initialize(dx);

		



		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0.5 * dx));
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_27(void) //unified 3d membrane with circle rim
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		ext_force = [&](const int i)->VectorD {
			VectorD pos;
			real r;
			VectorD result = VectorD::Zero();
			if (current_frame <= 10) {
				pos = fluid.particles.X(i) - (0.4 * VectorD::Unit(0));
				r = 0.3;
				if (pos.norm() < r) {
					result += (r - pos.norm()) * 10 * VectorD::Unit(0);
				}
			}
			/*pos = fluid.particles.X(i) - (-0.3 * VectorD::Unit(0));
			r = 0.2;
			if (pos.norm() < r) {
				result += (r - pos.norm()) * 10 * -VectorD::Unit(0);
			}*/
			return result * fluid.particles.M(i);
		};
	}
	else if constexpr (d == 3) {
		ext_force = [&](const int i)->VectorD {
			VectorD pos;
			real r;
			VectorD result = VectorD::Zero();
			pos = fluid.particles.X(i) - (0.4 * VectorD::Unit(0) + 0.4 * VectorD::Unit(2));
			r = 0.3;
			if (pos.norm() < r) {
				result += (r - pos.norm()) * 10 * VectorD::Unit(0);
			}

			pos = fluid.particles.X(i) - (-0.3 * VectorD::Unit(0) - 0.3 * VectorD::Unit(2));
			r = 0.2;
			if (pos.norm() < r) {
				result += (r - pos.norm()) * 10 * -VectorD::Unit(2);
			}

			//return result * 1e-5;
			return result * fluid.particles.M(i);
		};
	}
	real R = (real)1;

	real dx = 0.002; //this can be changed
	if constexpr (d == 2) {
		int num = int(R * 2.0 / dx);
		VectorD v1 = AuxFunc::V<d>(-R, 0), v2 = AuxFunc::V<d>(R, 0);
		PointSetFunc::Initialize_Segment_Points(v1, v2, num, fluid.particles);
		for (int i = 0; i < fluid.particles.Size(); i++) fluid.particles.B(i) = 0;
		fluid.particles.B(0) = fluid.particles.B(fluid.particles.Size() - 1) = 1;
	}
	else if constexpr (d == 3) {
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		for (int i = 0; i < fluid.particles.Size(); i++) fluid.particles.B(i) = is_boundary[i];
	}

	Set_Physical_Parameters();
	real rho = 1e3;
	fluid.thickness = 5e-7;
	real default_conc = 0;//"extreme soap"
	real total_surface = (d == 3) ? pi * pow(2 * R, 2) : 2 * R;
	//real total_mass = total_surface * rho, total_vol = total_surface * fluid.thickness;
	real total_vol = total_surface * fluid.thickness, total_mass = total_vol * rho;
	int n = fluid.particles.Size();
	std::cout << "single particle mass: " << total_mass / n << "\n";
	for (int i = 0; i < n; i++) {
		Initialize_Particle(i, total_mass / n, fluid.thickness);
		fluid.particles.Vol(i) = total_vol / n;
		fluid.particles.Conc(i) = default_conc;
	}	

	real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
	real alpha_height = 1.0 * rho;
	std::cout << "alpha vis: " << alpha_vis << "\n";
	std::cout << "alpha height: " << alpha_height << "\n";

	//world forces
	fluid.gravity_coeff = 0;
	fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;
	fluid.external_force_func = ext_force;

	//tangential force
	fluid.t_pressure_params.laplacian_pressure_coeff = 1 * 15000;
	fluid.t_pressure_params.divergence_pressure_coeff = 1 * 20*0;
	fluid.t_pressure_params.height_pressure_coeff = 1;
	fluid.t_pressure_params.tangential_pressure_coeff = 40;
	fluid.viscosity_coeff = 0.3 * alpha_vis;
	fluid.marangoni_coeff = 0;

	//normal force
	fluid.n_pressure_params.Set_Planar();

	fluid.boundary_params.Set_Pure_Analytical();
	fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0.5 * dx));

	//operator parameters
	fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
	fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);

	fluid.Initialize(dx);
}


template<int d>
void FluidSPHBubbleDriver<d>::Case_28(void) //unified 3d membrane with circle rim, scale up SIZE instead of FINENESS
{	
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
	}
	else if constexpr (d == 3) {
		ext_force = [&](const int i)->VectorD {
			VectorD pos;
			real r;
			VectorD result = VectorD::Zero();
			//if (fluid.current_frame > 100) return result;
			pos = fluid.particles.X(i) - (0.4 * fluid.simulation_scale *  VectorD::Unit(0) + 0.4 * fluid.simulation_scale * VectorD::Unit(2));
			r = 0.3 * fluid.simulation_scale;
			if (pos.norm() < r) {
				result += (r - pos.norm()) * 3 * VectorD::Unit(0);
			}

			pos = fluid.particles.X(i) - (-0.3 * fluid.simulation_scale * VectorD::Unit(0) - 0.3 * fluid.simulation_scale * VectorD::Unit(2));
			r = 0.3 * fluid.simulation_scale;
			if (pos.norm() < r) {
				result += (r - pos.norm()) * 3 * -VectorD::Unit(2);
			}

			//return result * 1e-5;
			//std::cout << "Fuck" << r << std::endl;
			return result * fluid.particles.M(i);
		};
	}
	if constexpr (d == 3) {
		real R = (real)0.2; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings
		fluid.simulation_scale = R;
		real dx = 0.005; //this can be changed
		fluid.np_on_h = 6;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";
		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
		}
		//basic geometry
		//fluid.replenish_boundary = false;
		//fluid.particle_geometry_mode = (d == 2) ? "geometric_compensate" : "mirror_compensate";//mirror compensate for d=3

		real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
		real alpha_height = 1.0 * rho;
		real alpha_height_laplacian = 1 * rho / (2 * fluid.gamma_water) / (fluid.thickness / (pow((fluid.thickness + 1), 2. / 3.)));
		real alpha_divergence = rho * frame_rate;
		std::cout << "alpha vis: " << alpha_vis << "\n";
		std::cout << "alpha height: " << alpha_height << "\n";
		std::cout << "alpha height laplacian: " << alpha_height << "\n";
		std::cout << "alpha divergence: " << alpha_height << "\n";

		//world forces
		fluid.gravity_coeff = 0;
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;
		fluid.external_force_func = ext_force;

		//tangential force
		fluid.t_pressure_params.Set_Weak1(rho, fluid.gamma_water, fluid.thickness, frame_rate);
		fluid.t_pressure_params.boundary_pressure_coeff = 0;
		fluid.divergence_params.calculate_field = "all";
		fluid.viscosity_coeff = 0.2 * alpha_vis;

		//Before changing to alpha_height laplacian and etc.
		//fluid.laplacian_pressure_coeff = 1 * 15000;
		//fluid.divergence_pressure_coeff = 1 * 20;
		//fluid.height_pressure_coeff = 1;
		//fluid.tangential_pressure_coeff = 40;
		//fluid.viscosity_coeff = 0.3 * alpha_vis;

		fluid.marangoni_coeff = 0;

		//normal force
		fluid.n_pressure_params.Set_Planar();

		//operator parameters
		fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
		fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);
		fluid.geometry_params = OperatorParams("only_fluid", KernelType::SPIKY);

		fluid.Initialize(dx);

		//boundary force
		//fluid.boundary_force_mode = "none";

		//control
		//fluid.boundary_mode = "analytical";

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0.5 * dx));
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_29(void) //unified 3d membrane with circle rim, scale up SIZE instead of FINENESS
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "Not serving d = 2, sorry!" << std::endl;
	}
	if constexpr (d == 3) {
		cfl = 0.01;

		real R = (real)0.1; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings

		fluid.simulation_scale = R;
		//frame_rate = (int)(5 / fluid.simulation_scale);
		real fineness = 1.0;
		real dx = 0.005/fineness; //this can be changed
		fluid.np_on_h = 6*fineness;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";

		//fluid.vorticity_params = VorticityParams(KernelType::QUINTIC, 1, true, 1e-4, true);
		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			//if (i == 363) {
			//	fluid.particles.Vrt(i) = 1 * fluid.vorticity_params.default_vort;
			//}
			//else if (i == 897) {
			//	fluid.particles.Vrt(i) = -1 * fluid.vorticity_params.default_vort;
			//}
		}

		real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
		real alpha_height = 1.0 * rho;
		real alpha_height_laplacian = 1 * rho / (2 * fluid.gamma_water) / (fluid.thickness / (pow((fluid.thickness + 1), 2. / 3.)));
		real alpha_divergence = rho * frame_rate;
		std::cout << "alpha vis: " << alpha_vis << "\n";
		std::cout << "alpha height: " << alpha_height << "\n";
		std::cout << "alpha height laplacian: " << alpha_height << "\n";
		std::cout << "alpha divergence: " << alpha_height << "\n";


		//world forces
		fluid.gravity_coeff = 0.1 * 0.1;
		//fluid.g = 0.1 * 9.8 * VectorD::Unit(0);
		fluid.g = 0.1/fluid.simulation_scale * 1 * 0.1 * 9.8 * VectorD::Unit(0);

		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;
		//fluid.external_force_func = ext_force;

		//tangential force
		fluid.t_pressure_params.Set_Baseline(rho, fluid.gamma_water, fluid.thickness, frame_rate);
		fluid.divergence_params.calculate_field = "all";
		fluid.viscosity_coeff = 0.1 * 2 * alpha_vis; //2
		fluid.normal_viscosity_coeff = 0;

		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(R, fluid.Surface_Tension_Coefficient(default_conc), rho, fluid.thickness);

		fluid.marangoni_coeff = 0;

		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);

		//operator parameters
		fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
		fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);
		fluid.geometry_params = OperatorParams("only_fluid", KernelType::SPIKY);

		

		fluid.Initialize(dx);

		

		//std::cout << "total surface area: " << total_sa << std::endl;
		//for (int i = 0; i < fluid.particles.Size();i++) {
		//	if (fluid.particles.X(i).norm() < 0.7) {
		//		fluid.particles.V(i) = -1 * fluid.particles.X(i);
		//		//std::cout << i << std::endl;
		//		//std::cout << particles.V(i) << std::endl;
		//	}
		//}

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0. * dx));
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_30(void)
{
	if constexpr (d == 3) {
		real R = (real)0.1; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings
		fluid.simulation_scale = R;
		real dx = 0.005; //this can be changed
		fluid.np_on_h = 6;

		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";
		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;

			//if (i == 631) {
			if(i==363){
				fluid.particles.Vrt(i) = 1e-3;
			}
			else if (i == 897) {
				fluid.particles.Vrt(i) = -1e-3;
			}
			/*if (fluid.particles.X(i).norm() <= R * 0.1) {
				fluid.particles.Vrt(i) = 1e-4;
			}*/
		}

		real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
		real alpha_height = 1.0 * rho;
		real alpha_height_laplacian = 1 * rho / (2 * fluid.gamma_water) / (fluid.thickness / (pow((fluid.thickness + 1), 2. / 3.)));
		real alpha_divergence = rho * frame_rate;
		std::cout << "alpha vis: " << alpha_vis << "\n";
		std::cout << "alpha height: " << alpha_height << "\n";
		std::cout << "alpha height laplacian: " << alpha_height << "\n";
		std::cout << "alpha divergence: " << alpha_height << "\n";

		//world forces
		fluid.gravity_coeff = 0;
		fluid.g = 0.3 * 9.8 * VectorD::Unit(0);
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;
		//fluid.external_force_func = ext_force;

		//tangential force
		fluid.vorticity_params = VorticityParams(KernelType::QUINTIC, 1, true, 1e-4);
		fluid.t_pressure_params.Set_Weak1(rho, fluid.gamma_water, fluid.thickness, frame_rate);
		fluid.divergence_params.calculate_field = "all";
		fluid.viscosity_coeff = 0.2 * alpha_vis;
		


		fluid.marangoni_coeff = 0;

		//normal force
		fluid.n_pressure_params.Set_Planar();

		fluid.boundary_params.Set_Pure_Analytical();

		//operator parameters
		fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
		fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);
		fluid.geometry_params = OperatorParams("only_fluid", KernelType::SPIKY);



		fluid.Initialize(dx);

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0.5 * dx));
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_31(void) //unified 3d membrane with circle rim, scale up SIZE instead of FINENESS
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
	}
	else if constexpr (d == 3) {
		ext_force = [&](const int i)->VectorD {
			VectorD pos;
			real r;
			VectorD result = VectorD::Zero();
			//if (fluid.current_frame > 100) return result;
			pos = fluid.particles.X(i) - (0.4 * fluid.simulation_scale * VectorD::Unit(0) + 0.4 * fluid.simulation_scale * VectorD::Unit(2));
			r = 0.3 * fluid.simulation_scale;
			if (pos.norm() < r) {
				result += (r - pos.norm()) * 3 * VectorD::Unit(0);
			}

			pos = fluid.particles.X(i) - (-0.3 * fluid.simulation_scale * VectorD::Unit(0) - 0.3 * fluid.simulation_scale * VectorD::Unit(2));
			r = 0.3 * fluid.simulation_scale;
			if (pos.norm() < r) {
				result += (r - pos.norm()) * 3 * -VectorD::Unit(2);
			}

			//return result * 1e-5;
			//std::cout << "Fuck" << r << std::endl;
			return result * fluid.particles.M(i);
		};
	}
	if constexpr (d == 3) {
		real R = (real)0.1; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings
		fluid.simulation_scale = R;
		real dx = 0.005; //this can be changed
		fluid.np_on_h = 6;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";
		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
		}

		real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
		real alpha_height = 1.0 * rho;
		real alpha_height_laplacian = 1 * rho / (2 * fluid.gamma_water) / (fluid.thickness / (pow((fluid.thickness + 1), 2. / 3.)));
		real alpha_divergence = rho * frame_rate;
		std::cout << "alpha vis: " << alpha_vis << "\n";
		std::cout << "alpha height: " << alpha_height << "\n";
		std::cout << "alpha height laplacian: " << alpha_height << "\n";
		std::cout << "alpha divergence: " << alpha_height << "\n";

		//world forces
		fluid.gravity_coeff = 0.1 * fluid.simulation_scale;
		//fluid.g = 0.1 * 9.8 * VectorD::Unit(0);
		fluid.g = 0.3 * -9.8 * VectorD::Unit(1);

		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;
		//fluid.external_force_func = ext_force;

		//tangential force
		fluid.t_pressure_params.Set_Zero();
		fluid.divergence_params.calculate_field = "all";
		fluid.viscosity_coeff = 1 * 2 * alpha_vis;

		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(R, fluid.Surface_Tension_Coefficient(default_conc), rho, fluid.thickness);
		fluid.n_pressure_params.capillary_coeff *= 1.5;

		fluid.marangoni_coeff = 0;

		//operator parameters
		fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
		fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);
		fluid.geometry_params = OperatorParams("only_fluid", KernelType::SPIKY);

		fluid.boundary_params.Set_Pure_Gravity(0.0000005 / (total_mass / n));
		fluid.boundary_params.boundary_mode = "slippery";

		fluid.Initialize(dx);


		//fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0.5 * dx));
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_32(void) //unified 3d membrane with circle rim, scale up SIZE instead of FINENESS
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
	}
	else if constexpr (d == 3) {
		ext_force = [&](const int i)->VectorD {
			VectorD pos;
			real r;
			VectorD result = VectorD::Zero();
			//if (fluid.current_frame > 100) return result;
			pos = fluid.particles.X(i) - (0.4 * fluid.simulation_scale * VectorD::Unit(0) + 0.4 * fluid.simulation_scale * VectorD::Unit(2));
			r = 0.3 * fluid.simulation_scale;
			if (pos.norm() < r) {
				result += (r - pos.norm()) * 3 * VectorD::Unit(0);
			}

			pos = fluid.particles.X(i) - (-0.3 * fluid.simulation_scale * VectorD::Unit(0) - 0.3 * fluid.simulation_scale * VectorD::Unit(2));
			r = 0.3 * fluid.simulation_scale;
			if (pos.norm() < r) {
				result += (r - pos.norm()) * 3 * -VectorD::Unit(2);
			}

			//return result * 1e-5;
			//std::cout << "Fuck" << r << std::endl;
			return result * fluid.particles.M(i);
		};
	}
	if constexpr (d == 3) {
		real R = (real)0.15; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings
		fluid.simulation_scale = R;
		real dx = 0.005; //this can be changed
		fluid.np_on_h = 6;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";
		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
		}

		real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
		real alpha_height = 1.0 * rho;
		real alpha_height_laplacian = 1 * rho / (2 * fluid.gamma_water) / (fluid.thickness / (pow((fluid.thickness + 1), 2. / 3.)));
		real alpha_divergence = rho * frame_rate;
		std::cout << "alpha vis: " << alpha_vis << "\n";
		std::cout << "alpha height: " << alpha_height << "\n";
		std::cout << "alpha height laplacian: " << alpha_height << "\n";
		std::cout << "alpha divergence: " << alpha_height << "\n";

		//world forces
		fluid.gravity_coeff = 0.1 * fluid.simulation_scale;
		//fluid.g = 0.1 * 9.8 * VectorD::Unit(0);
		fluid.g = 1 * 0.3 * -9.8 * VectorD::Unit(1);

		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;
		//fluid.external_force_func = ext_force;

		//tangential force
		fluid.t_pressure_params.Set_Weak2(rho, fluid.gamma_water, fluid.thickness, frame_rate);
		fluid.divergence_params.calculate_field = "all";
		
		fluid.viscosity_coeff = 1 * 2 * alpha_vis;

		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(R, fluid.Surface_Tension_Coefficient(default_conc), rho, fluid.thickness);
		fluid.n_pressure_params.capillary_coeff *= 1.5;

		fluid.marangoni_coeff = 0;

		fluid.boundary_params.Set_Pure_Gravity(1 * 0.0000005 * 1. / (total_mass / n));


		//operator parameters
		fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
		fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);
		fluid.geometry_params= OperatorParams("only_fluid", KernelType::SPIKY);

		fluid.Initialize(dx);

		//std::cout << "total surface area: " << total_sa << std::endl;
		//for (int i = 0; i < fluid.particles.Size();i++) {
		//	if (fluid.particles.X(i).norm() < 0.7) {
		//		fluid.particles.V(i) = -1 * fluid.particles.X(i);
		//		//std::cout << i << std::endl;
		//		//std::cout << particles.V(i) << std::endl;
		//	}
		//}

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0.5 * dx));
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_33(void) //unified 3d membrane with circle rim, scale up SIZE instead of FINENESS
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
	}
	else if constexpr (d == 3) {
		ext_force = [&](const int i)->VectorD {
			VectorD pos;
			real r;
			VectorD result = VectorD::Zero();
			//if (fluid.current_frame > 100) return result;
			pos = fluid.particles.X(i) - (0.4 * fluid.simulation_scale * VectorD::Unit(0) + 0.4 * fluid.simulation_scale * VectorD::Unit(2));
			r = 0.3 * fluid.simulation_scale;
			if (pos.norm() < r) {
				result += (r - pos.norm()) * 3 * VectorD::Unit(0);
			}

			pos = fluid.particles.X(i) - (-0.3 * fluid.simulation_scale * VectorD::Unit(0) - 0.3 * fluid.simulation_scale * VectorD::Unit(2));
			r = 0.3 * fluid.simulation_scale;
			if (pos.norm() < r) {
				result += (r - pos.norm()) * 3 * -VectorD::Unit(2);
			}

			//return result * 1e-5;
			//std::cout << "Fuck" << r << std::endl;
			return result * fluid.particles.M(i);
		};
	}
	if constexpr (d == 3) {
		real R = (real)0.1; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings
		fluid.simulation_scale = R;
		real dx = 0.005; //this can be changed
		fluid.np_on_h = 6;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";
		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
		}

		real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
		real alpha_height = 1.0 * rho;
		real alpha_height_laplacian = 1 * rho / (2 * fluid.gamma_water) / (fluid.thickness / (pow((fluid.thickness + 1), 2. / 3.)));
		real alpha_divergence = rho * frame_rate;
		std::cout << "alpha vis: " << alpha_vis << "\n";
		std::cout << "alpha height: " << alpha_height << "\n";
		std::cout << "alpha height laplacian: " << alpha_height << "\n";
		std::cout << "alpha divergence: " << alpha_height << "\n";

		//world forces
		fluid.gravity_coeff = 0.1 * fluid.simulation_scale;
		//fluid.g = 1 * 0.3 * -9.8 * VectorD::Unit(1);
		fluid.g = 1 * 0.3 * 9.8 * VectorD::Unit(0);
		//fluid.g = 0 * VectorD::Unit(1);

		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;
		//fluid.external_force_func = ext_force;

		//tangential force
		fluid.t_pressure_params.Set_Weak2(rho, fluid.gamma_water, fluid.thickness, frame_rate);
		fluid.divergence_params.calculate_field = "all";
		fluid.viscosity_coeff = 1 * 2 * alpha_vis;

		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(R, fluid.Surface_Tension_Coefficient(default_conc), rho, fluid.thickness);
		fluid.n_pressure_params.capillary_coeff *= 30;

		fluid.marangoni_coeff = 0;

		fluid.boundary_params.Set_Pure_Gravity(1 * 0.0000008 * 1. / (total_mass / n));

		//operator parameters
		fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
		fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);
		fluid.geometry_params = OperatorParams("only_fluid", KernelType::SPIKY);


		fluid.Initialize(dx);

		////std::cout << "total surface area: " << total_sa << std::endl;
		//for (int i = 0; i < fluid.particles.Size();i++) {
		//	if (fluid.particles.X(i).norm() < 0.7 * fluid.simulation_scale) {
		//		fluid.particles.V(i) = -1 * fluid.particles.X(i);
		//		//std::cout << i << std::endl;
		//		//std::cout << particles.V(i) << std::endl;
		//	}
		//}

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0.5 * dx));
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_34(void) //unified 3d membrane with circle rim, scale up SIZE instead of FINENESS
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "Not serving d = 2, sorry!" << std::endl;
	}
	if constexpr (d == 3) {
		max_iter_per_frame = -1;
		cfl = 0.1;
		fluid.max_vel = 0.1;

		real R = (real)0.1; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings

		fluid.simulation_scale = R;
		frame_rate = 50.;
		real fineness = 1.0;
		real dx = 0.005 / fineness; //this can be changed
		fluid.np_on_h = 6 * fineness;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";

		//fluid.vorticity_params = VorticityParams(KernelType::QUINTIC, 1, true, 1e-4, true);
		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			//if (i == 363) {
			//	fluid.particles.Vrt(i) = 1 * fluid.vorticity_params.default_vort;
			//}
			//else if (i == 897) {
			//	fluid.particles.Vrt(i) = -1 * fluid.vorticity_params.default_vort;
			//}
		}
		//basic geometry
		//fluid.replenish_boundary = true;

		// grav--unscaled
		fluid.g = 9.8 * VectorD::Unit(0);

		//external force
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;

		//main params
		fluid.default_sim_params = DefaultSimParams(dx, rho, fluid.viscosity_water, fluid.gamma_water, fluid.thickness, fluid.simulation_scale);

		//world forces
		fluid.gravity_coeff = pow(.1/R, 2) * fluid.default_sim_params.gravity_coeff;

		//tangential force
		fluid.t_pressure_params.Set_Baseline(rho, fluid.gamma_water, fluid.thickness, 50);
		fluid.viscosity_coeff = 1*fluid.default_sim_params.viscosity_coeff;
		fluid.normal_viscosity_coeff = fluid.default_sim_params.normal_viscosity_coeff;

		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(R, fluid.gamma_water, rho, fluid.thickness);

		fluid.marangoni_coeff = fluid.default_sim_params.marangoni_coeff;

		fluid.n_pressure_params.capillary_coeff *= pow(.1 / R, 2); /// to adapt to the reduced G

		//operator parameters
		fluid.grad_force_params = fluid.default_sim_params.grad_force_params;
		fluid.height_laplacian_params = fluid.default_sim_params.height_laplacian_params;
		fluid.geometry_params = fluid.default_sim_params.geometry_params;

		//fluid.t_pressure_params.divergence_pressure_coeff = 0.;
		fluid.t_pressure_params.boundary_pressure_coeff *= 1000.;
		fluid.boundary_params.vis_boundary = 0.;

		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0. * dx));

		fluid.Initialize(dx);
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_35(void) //play with marangoni
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "Not serving d = 2, sorry!" << std::endl;
	}
	if constexpr (d == 3) {
		max_iter_per_frame = -1;
		cfl = 0.1;
		fluid.max_vel = 0.1;

		real R = (real)0.3; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings

		fluid.simulation_scale = R;
		frame_rate = (50. / 8.);
		real fineness = 1.0;
		real dx = 0.005 / fineness; //this can be changed
		fluid.np_on_h = 6 * fineness;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";

		//fluid.vorticity_params = VorticityParams(KernelType::QUINTIC, 1, true, 1e-4, true);
		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			//if (i == 363) {
			//	fluid.particles.Vrt(i) = 1 * fluid.vorticity_params.default_vort;
			//}
			//else if (i == 897) {
			//	fluid.particles.Vrt(i) = -1 * fluid.vorticity_params.default_vort;
			//}
		}
		//basic geometry
		//fluid.replenish_boundary = true;

		// grav--unscaled
		fluid.g = 0 * 9.8 * VectorD::Unit(0);

		//external force
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;

		//main params
		fluid.default_sim_params = DefaultSimParams(dx, rho, fluid.viscosity_water, fluid.gamma_water, fluid.thickness, fluid.simulation_scale);

		//world forces
		fluid.gravity_coeff = pow(.1 / R, 2) * fluid.default_sim_params.gravity_coeff;

		//tangential force
		fluid.t_pressure_params.Set_Baseline(rho, fluid.gamma_water, fluid.thickness, 50);
		fluid.viscosity_coeff = 0*fluid.default_sim_params.viscosity_coeff;
		fluid.normal_viscosity_coeff = fluid.default_sim_params.normal_viscosity_coeff;

		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(R, fluid.gamma_water, rho, fluid.thickness);

		fluid.marangoni_coeff = fluid.default_sim_params.marangoni_coeff;

		fluid.n_pressure_params.capillary_coeff *= pow(.1 / R, 2); /// to adapt to the reduced G

		//operator parameters
		fluid.grad_force_params = fluid.default_sim_params.grad_force_params;
		fluid.height_laplacian_params = fluid.default_sim_params.height_laplacian_params;
		fluid.geometry_params = fluid.default_sim_params.geometry_params;

		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0. * dx));

		fluid.Initialize(dx);

		for (int i = 0; i < fluid.particles.Size(); i++) {
			if (fluid.particles.X(i).norm() < R) {
				fluid.particles.Conc(i) = Rand_Number();
			}
		}
		fluid.diffuse_soap = true;
		fluid.marangoni_coeff = 0.0001;
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_36(void) //play with different init volume
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "Not serving d = 2, sorry!" << std::endl;
	}
	if constexpr (d == 3) {
		max_iter_per_frame = -1;
		cfl = 0.1;
		fluid.max_vel = 0.1;

		real R = (real)0.3; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings

		fluid.simulation_scale = R;
		frame_rate = 50.;
		real fineness = 1.0;
		real dx = 0.005 / fineness; //this can be changed
		fluid.np_on_h = 6 * fineness;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";

		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n,fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			if (fluid.particles.X(i)[0] < 0.) {
				fluid.particles.M(i) *= 2;
				fluid.particles.Vol(i) *= 2;
			}
		}
		//basic geometry
		//fluid.replenish_boundary = true;

		// grav--unscaled
		fluid.g = 9.8 * VectorD::Unit(0);
		fluid.g *= 10;

		//external force
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;

		//main params
		fluid.default_sim_params = DefaultSimParams(dx, rho, fluid.viscosity_water, fluid.gamma_water, fluid.thickness, fluid.simulation_scale);

		//world forces
		fluid.gravity_coeff = pow(.1 / R, 2) * fluid.default_sim_params.gravity_coeff;

		//tangential force
		fluid.t_pressure_params.Set_Baseline(rho, fluid.gamma_water, fluid.thickness, 50);
		fluid.viscosity_coeff = 0.01 * fluid.default_sim_params.viscosity_coeff;
		fluid.t_pressure_params.divergence_pressure_coeff = 0.;
		fluid.normal_viscosity_coeff = fluid.default_sim_params.normal_viscosity_coeff;

		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(R, fluid.gamma_water, rho, fluid.thickness);

		fluid.marangoni_coeff = fluid.default_sim_params.marangoni_coeff;

		fluid.n_pressure_params.capillary_coeff *= pow(.1 / R, 2); /// to adapt to the reduced G

		//operator parameters
		fluid.grad_force_params = fluid.default_sim_params.grad_force_params;
		fluid.height_laplacian_params = fluid.default_sim_params.height_laplacian_params;
		fluid.geometry_params = fluid.default_sim_params.geometry_params;

		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0. * dx));

		fluid.Initialize(dx);
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_37(void) //play with different init density
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "Not serving d = 2, sorry!" << std::endl;
	}
	if constexpr (d == 3) {
		max_iter_per_frame = -1;
		cfl = 0.1;
		fluid.max_vel = 0.1;

		real R = (real)0.15; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings

		fluid.simulation_scale = R;
		frame_rate = 50.;
		real fineness = 1;
		real dx = 0.005 / fineness; //this can be changed
		fluid.np_on_h = 6 * fineness;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";

		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			if (fluid.particles.X(i)[0] < 0.) {
				fluid.particles.M(i) *= 4;
				//fluid.particles.Vol(i) *= (1 + 3*Rand_Number());
				fluid.particles.RH(i) = 8e-7;
			}
			else {
				fluid.particles.M(i) *= 1;
				//fluid.particles.Vol(i) *= (1 + 1 * Rand_Number());
				fluid.particles.RH(i) = 3e-7;
			}

		}
		//basic geometry
		//fluid.replenish_boundary = true;

		// grav--unscaled
		fluid.g = 9.8 * VectorD::Unit(0);
		fluid.g *= 10;

		//external force
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;

		//main params
		fluid.default_sim_params = DefaultSimParams(dx, rho, fluid.viscosity_water, fluid.gamma_water, fluid.thickness, fluid.simulation_scale);

		//world forces
		fluid.gravity_coeff = pow(.1 / R, 2) * fluid.default_sim_params.gravity_coeff;

		//tangential force
		fluid.t_pressure_params.Set_Baseline(rho, fluid.gamma_water, fluid.thickness, 50);

		fluid.t_pressure_params.boundary_pressure_coeff *= 20;
		fluid.t_pressure_params.height_pressure_coeff = 30;
		fluid.t_pressure_params.laplacian_pressure_coeff = 0;
		fluid.t_pressure_params.divergence_pressure_coeff = 50;

		fluid.viscosity_coeff = 1.0 * fluid.default_sim_params.viscosity_coeff;
		fluid.normal_viscosity_coeff = fluid.default_sim_params.normal_viscosity_coeff;

		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(R, fluid.gamma_water, rho, fluid.thickness);

		fluid.marangoni_coeff = fluid.default_sim_params.marangoni_coeff;

		fluid.n_pressure_params.capillary_coeff *= pow(.1 / R, 2); /// to adapt to the reduced G

		//operator parameters
		fluid.grad_force_params = fluid.default_sim_params.grad_force_params;
		fluid.height_laplacian_params = fluid.default_sim_params.height_laplacian_params;
		fluid.geometry_params = fluid.default_sim_params.geometry_params;

		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
		fluid.boundary_params.replenish = true;

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0. * dx));

		fluid.rh_params = RenderHeightParams("divergence", true, false, 0);

		fluid.Initialize(dx);
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_38(void) //play with different init density
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "Not serving d = 2, sorry!" << std::endl;
	}
	if constexpr (d == 3) {
		max_iter_per_frame = -1;
		cfl = 0.1;
		fluid.max_vel = 0.1;

		real R = (real)1.; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings

		fluid.simulation_scale = R;
		//frame_rate = 50.;
		frame_rate = 10.;
		real fineness = 1;
		real dx = 0.005 / fineness; //this can be changed
		fluid.np_on_h = 6 * fineness;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";

		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			fluid.particles.RH(i) = fluid.thickness;
		}

		//basic geometry
		//fluid.replenish_boundary = true;

		// grav--unscaled
		fluid.g = 9.8 * VectorD::Unit(0);

		//external force
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;

		//main params
		fluid.default_sim_params = DefaultSimParams(dx, rho, fluid.viscosity_water, fluid.gamma_water, fluid.thickness, fluid.simulation_scale);

		//These three must go together
		fluid.gravity_coeff = 0.33 / R * 0.6 * 10 * pow(.1 / 0.15, 2) * fluid.default_sim_params.gravity_coeff;
		fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, 50);
		fluid.viscosity_coeff = 2 * fluid.default_sim_params.viscosity_coeff;
		//

		fluid.normal_viscosity_coeff = fluid.default_sim_params.normal_viscosity_coeff;

		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(0.1, fluid.gamma_water, rho, fluid.thickness);
		fluid.n_pressure_params.capillary_coeff *= 0.33;

		fluid.marangoni_coeff = fluid.default_sim_params.marangoni_coeff;

		//fluid.n_pressure_params.capillary_coeff *= R / 0.15; /// to adapt to the reduced G

		//operator parameters
		fluid.grad_force_params = fluid.default_sim_params.grad_force_params;
		fluid.height_laplacian_params = fluid.default_sim_params.height_laplacian_params;
		fluid.geometry_params = fluid.default_sim_params.geometry_params;

		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
		fluid.boundary_params.replenish = false;

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0. * dx));

		fluid.rh_params = RenderHeightParams("divergence", false, true, 0);

		fluid.Initialize(dx);
		fluid.clip_velocity = true;
		fluid.vel_threshold = 0.1; //50 frames per second, one frame max 

		//fluid.g_init = 9.8 * VectorD::Unit(0);
		//fluid.g_final = 9.8 * -VectorD::Unit(1);
		//fluid.interp_g = true;
		//fluid.first_frame = first_frame;
		//fluid.last_frame = last_frame;


		if (first_frame < 1 && init_snapshot_name != "none") {
			bool loaded;
			loaded = fluid.particles.Load_Snapshot(init_snapshot_name);
			if (loaded) {
				int n = fluid.particles.Size();
				//for (int i = 0; i < n; i++) {
				//	fluid.particles.RH(i) = fluid.thickness;
				//	fluid.particles.RH_V(i) = 0.;
				//	fluid.particles.Div(i) = 0.;
				//	fluid.particles.V(i) *= 0.;
				//}
			}
			else {
				std::cout << "Unable to Load From Snapshot!" << std::endl;
				exit(0);
			}
		}
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_39(void) //play with different init density
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "Not serving d = 2, sorry!" << std::endl;
	}
	if constexpr (d == 3) {
		max_iter_per_frame = -1;
		cfl = 0.1;
		fluid.max_vel = 0.1;

		real R = (real).33; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings

		fluid.simulation_scale = R;
		frame_rate = 50.;
		real fineness = 1;
		real dx = 0.005 / fineness; //this can be changed
		fluid.np_on_h = 6 * fineness;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";

		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			fluid.particles.RH(i) = fluid.thickness;
		}

		//basic geometry
		//fluid.replenish_boundary = true;

		// grav--unscaled
		fluid.g = 0 * 9.8 * VectorD::Unit(0);

		//external force
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;

		//main params
		fluid.default_sim_params = DefaultSimParams(dx, rho, fluid.viscosity_water, fluid.gamma_water, fluid.thickness, fluid.simulation_scale);

		//These three must go together
		fluid.gravity_coeff = 0.33 / R * 0.6 * 10 * pow(.1 / 0.15, 2) * fluid.default_sim_params.gravity_coeff;
		fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, 50);
		fluid.viscosity_coeff = 2 * fluid.default_sim_params.viscosity_coeff;
		//

		fluid.normal_viscosity_coeff = fluid.default_sim_params.normal_viscosity_coeff;

		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(0.33, fluid.gamma_water, rho, fluid.thickness);

		fluid.marangoni_coeff = fluid.default_sim_params.marangoni_coeff;

		//fluid.n_pressure_params.capillary_coeff *= R / 0.15; /// to adapt to the reduced G

		//operator parameters
		fluid.grad_force_params = fluid.default_sim_params.grad_force_params;
		fluid.height_laplacian_params = fluid.default_sim_params.height_laplacian_params;
		fluid.geometry_params = fluid.default_sim_params.geometry_params;

		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
		fluid.boundary_params.replenish = false;

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0. * dx));

		fluid.rh_params = RenderHeightParams("divergence", false, true, 0);

		fluid.Initialize(dx);
		fluid.clip_velocity = true;
		fluid.vel_threshold = 0.1; //50 frames per second, one frame max 

		if (first_frame < 1 && init_snapshot_name != "none") {
			bool loaded;
			loaded = fluid.particles.Load_Snapshot(init_snapshot_name);
			if (loaded) {
				int n = fluid.particles.Size();
				for (int i = 0; i < n; i++) {
					fluid.particles.RH(i) = fluid.thickness;
					fluid.particles.RH_V(i) = 0.;
					fluid.particles.Div(i) = 0.;
					fluid.particles.V(i) *= 0.;
				}
			}
			else {
				std::cout << "Unable to Load From Snapshot!" << std::endl;
				exit(0);
			}
		}
	}
}


template<int d>
void FluidSPHBubbleDriver<d>::Case_40(void)
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "Not serving d = 2, sorry!" << std::endl;
	}
	if constexpr (d == 3) {
		max_iter_per_frame = -1;
		cfl = 0.1;
		fluid.max_vel = 0.1;

		real R = (real).15; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings

		fluid.simulation_scale = R;
		frame_rate = 50.;
		real fineness = 1;
		real dx = 0.005 / fineness; //this can be changed
		fluid.np_on_h = 6 * fineness;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";

		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			fluid.particles.RH(i) = fluid.thickness;
		}

		//basic geometry
		//fluid.replenish_boundary = true;

		// grav--unscaled
		fluid.g = 9.8 * VectorD::Unit(0);

		//external force
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;

		//main params
		fluid.default_sim_params = DefaultSimParams(dx, rho, fluid.viscosity_water, fluid.gamma_water, fluid.thickness, fluid.simulation_scale);

		//These three must go together
		fluid.gravity_coeff = 0.33 / R * 0.6 * 10 * pow(.1 / 0.15, 2) * fluid.default_sim_params.gravity_coeff;
		fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, 50);
		fluid.viscosity_coeff = 2 * fluid.default_sim_params.viscosity_coeff;
		//

		fluid.normal_viscosity_coeff = fluid.default_sim_params.normal_viscosity_coeff;

		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(0.33, fluid.gamma_water, rho, fluid.thickness);

		fluid.marangoni_coeff = fluid.default_sim_params.marangoni_coeff;

		//fluid.n_pressure_params.capillary_coeff *= R / 0.15; /// to adapt to the reduced G

		//operator parameters
		fluid.grad_force_params = fluid.default_sim_params.grad_force_params;
		fluid.height_laplacian_params = fluid.default_sim_params.height_laplacian_params;
		fluid.geometry_params = fluid.default_sim_params.geometry_params;

		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
		fluid.boundary_params.replenish = false;

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0. * dx));

		fluid.rh_params = RenderHeightParams("divergence", false, true, 0);

		fluid.Initialize(dx);
		fluid.clip_velocity = true;
		fluid.vel_threshold = 0.1; //50 frames per second, one frame max 

		if (first_frame < 1 && init_snapshot_name != "none") {
			bool loaded;
			loaded = fluid.particles.Load_Snapshot(init_snapshot_name);
			if (loaded) {
				int n = fluid.particles.Size();
				for (int i = 0; i < n; i++) {
					fluid.particles.RH(i) = fluid.thickness;
					fluid.particles.RH_V(i) = 0.;
					fluid.particles.Div(i) = 0.;
					fluid.particles.V(i) *= 0.;
				}
			}
			else {
				std::cout << "Unable to Load From Snapshot!" << std::endl;
				exit(0);
			}
		}

		////add perlin noise to stuff
		std::uint32_t seed = 123;
		real perlin_freq;
		std::uint32_t octaves = 2;
		const siv::PerlinNoise perlin(seed);
		real perlin_inc;

		if (first_frame < 1) {
			for (int i = 0; i < fluid.particles.Size(); i++) {
				perlin_freq = 5.;
				perlin_inc = perlin.accumulatedOctaveNoise3D(fluid.particles.X(i)[0] * perlin_freq, fluid.particles.X(i)[1] * perlin_freq, fluid.particles.X(i)[2] * perlin_freq, octaves);
				real perlin_multiplier_V = std::max(0.1, 1 + (perlin_inc * 1.));
				//real perlin_multiplier_V = 1.;
				fluid.particles.Vol(i) *= perlin_multiplier_V;
				fluid.particles.RH(i) *= perlin_multiplier_V;
				std::cout << "Multiplier Perlina: " << std::max(0.1, 1 + (perlin_inc * 1)) << std::endl;
								//add mass perlin noise
				perlin_freq = 20.;
				perlin_inc = perlin.accumulatedOctaveNoise3D(fluid.particles.X(i)[0] * perlin_freq, fluid.particles.X(i)[1] * perlin_freq, fluid.particles.X(i)[2] * perlin_freq, octaves);
				real perlin_multiplier_density = std::max(0.1, 1 + (perlin_inc * 1.));
				std::cout << "Multiplier Perlina: " << std::max(0.7, 1 + (perlin_inc * 2.)) << std::endl;
				//add mass perlin noise
				fluid.particles.M(i) *= perlin_multiplier_V * perlin_multiplier_density;
			}
		}
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_41(void) //random init
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "Not serving d = 2, sorry!" << std::endl;
	}
	if constexpr (d == 3) {
		max_iter_per_frame = -1;
		cfl = 0.1;
		fluid.max_vel = 0.1;

		real R = (real)1.; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings

		fluid.simulation_scale = R;
		frame_rate = 50.;
		real fineness = 1;
		real dx = 0.005 / fineness; //this can be changed
		fluid.np_on_h = 6 * fineness;

		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";

		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			fluid.particles.RH(i) = fluid.thickness;
		}

		//basic geometry
		//fluid.replenish_boundary = true;

		// grav--unscaled
		fluid.g = 9.8 * VectorD::Unit(0);

		//external force
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;

		//main params
		fluid.default_sim_params = DefaultSimParams(dx, rho, fluid.viscosity_water, fluid.gamma_water, fluid.thickness, fluid.simulation_scale);

		//These three must go together
		fluid.gravity_coeff = 0.33 / R * 0.6 * 10 * pow(.1 / 0.15, 2) * fluid.default_sim_params.gravity_coeff;
		fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, 50);
		fluid.viscosity_coeff = 2 * fluid.default_sim_params.viscosity_coeff;
		//

		fluid.normal_viscosity_coeff = fluid.default_sim_params.normal_viscosity_coeff;

		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(0.33, fluid.gamma_water, rho, fluid.thickness);

		fluid.marangoni_coeff = fluid.default_sim_params.marangoni_coeff;

		//fluid.n_pressure_params.capillary_coeff *= R / 0.15; /// to adapt to the reduced G

		//operator parameters
		fluid.grad_force_params = fluid.default_sim_params.grad_force_params;
		fluid.height_laplacian_params = fluid.default_sim_params.height_laplacian_params;
		fluid.geometry_params = fluid.default_sim_params.geometry_params;

		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
		fluid.boundary_params.replenish = false;

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0. * dx));

		fluid.rh_params = RenderHeightParams("divergence", false, true, 0);

		fluid.Initialize(dx);
		fluid.clip_velocity = true;
		fluid.vel_threshold = 0.1; //50 frames per second, one frame max 
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_42(void) //play with different init density
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "Not serving d = 2, sorry!" << std::endl;
	}
	if constexpr (d == 3) {
		max_iter_per_frame = -1;
		cfl = 0.1;
		fluid.max_vel = 0.1;

		real R = (real)1.; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings

		fluid.simulation_scale = R;
		//frame_rate = 50.;
		frame_rate = 10.;
		real fineness = 1;
		real dx = 0.005 / fineness; //this can be changed
		fluid.np_on_h = 6 * fineness;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";

		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			fluid.particles.RH(i) = fluid.thickness;
			if (fluid.particles.X(i)[0] < 0.) {
				fluid.particles.M(i) *= 4;
			}
			else {
				fluid.particles.Vol(i) *= 1.5;
				fluid.particles.RH(i) *= 4;
			}
		}

		//basic geometry
		//fluid.replenish_boundary = true;
		//fluid.delete_idle = true;

		// grav--unscaled
		fluid.g = 0*1.5 * 9.8 * VectorD::Unit(0);

		//external force
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;

		//main params
		fluid.default_sim_params = DefaultSimParams(dx, rho, fluid.viscosity_water, fluid.gamma_water, fluid.thickness, fluid.simulation_scale);

		//std::function<VectorD(const int)> ext_force = [&](const int i)->VectorD {
		//	const real g_strength = 0.3 * 9.8 * 0.33 / fluid.simulation_scale * 0.6 * 10 * pow(.1 / 0.15, 2) * fluid.default_sim_params.gravity_coeff;
		//	const real cycle = 200;
		//	real alpha = (current_frame + 0.0) / cycle * 2 * pi; // 100 is the cycle
		//	return (AuxFunc::V<d>(cos(alpha), sin(alpha), 0) + AuxFunc::V<d>(0, cos(2*alpha), sin(2*alpha) * 0.6)) * g_strength * fluid.particles.M(i);
		//};
		std::function<VectorD(const int)> ext_force = [&](const int i)->VectorD {
			real along_normal = fluid.particles.V(i).dot(fluid.particles.Normal(i));
			return -0.1 * fluid.particles.M(i) * along_normal * fluid.particles.Normal(i);
		};
		fluid.external_force_func = ext_force;

		//These three must go together
		fluid.gravity_coeff = 0.33 / R * 0.6 * 10 * pow(.1 / 0.15, 2) * fluid.default_sim_params.gravity_coeff;
		fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, 50);
		fluid.viscosity_coeff = 0.6 * fluid.default_sim_params.viscosity_coeff;
		fluid.t_pressure_params.height_pressure_coeff *= 5;
		fluid.t_pressure_params.laplacian_pressure_coeff *= 5;
		fluid.t_pressure_params.divergence_pressure_coeff *= 0.3;
		//

		fluid.normal_viscosity_coeff = fluid.default_sim_params.normal_viscosity_coeff;

		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(0.33, fluid.gamma_water, rho, fluid.thickness);
		fluid.n_pressure_params.capillary_coeff *= 0.33;

		fluid.marangoni_coeff = fluid.default_sim_params.marangoni_coeff;

		//fluid.n_pressure_params.capillary_coeff *= R / 0.15; /// to adapt to the reduced G

		//operator parameters
		fluid.grad_force_params = fluid.default_sim_params.grad_force_params;
		fluid.height_laplacian_params = fluid.default_sim_params.height_laplacian_params;
		fluid.geometry_params = fluid.default_sim_params.geometry_params;

		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
		fluid.boundary_params.replenish = false;
		//fluid.boundary_params.replenish = true;
		//fluid.boundary_params.inherit_rh = true;
		//fluid.boundary_params.replenish_interval = 1e-3;
		//fluid.boundary_params.replenish_dx_num = 3.;

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0. * dx));

		fluid.rh_params = RenderHeightParams("divergence", false, true, 0);

		fluid.Initialize(dx);
		fluid.clip_velocity = true;
		fluid.vel_threshold = 0.1; //50 frames per second, one frame max 

		//fluid.g_init = 9.8 * VectorD::Unit(0);
		//fluid.g_final = 9.8 * -VectorD::Unit(1);
		//fluid.interp_g = true;
		//fluid.first_frame = first_frame;
		//fluid.last_frame = last_frame;


		if (first_frame < 1 && init_snapshot_name != "none") {
			bool loaded;
			loaded = fluid.particles.Load_Snapshot(init_snapshot_name);
			if (loaded) {
				int n = fluid.particles.Size();
				for (int i = 0; i < n; i++) {
					//fluid.particles.RH(i) = fluid.thickness;
					//fluid.particles.RH_V(i) = 0.;
					//fluid.particles.Div(i) = 0.;
					//fluid.particles.V(i) *= 0.;
				}
			}
			else {
				std::cout << "Unable to Load From Snapshot!" << std::endl;
				exit(0);
			}
		}
	}
}


template<int d>
void FluidSPHBubbleDriver<d>::Case_43(void)
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "Not serving d = 2, sorry!" << std::endl;
	}
	if constexpr (d == 3) {
		max_iter_per_frame = -1;
		cfl = 0.1;
		fluid.max_vel = 0.1;

		real R = (real)1.; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings

		fluid.simulation_scale = R;
		frame_rate = 50.;
		real fineness = 1;
		real dx = 0.005 / fineness; //this can be changed
		fluid.np_on_h = 6 * fineness;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";

		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			fluid.particles.RH(i) = fluid.thickness;
		}

		//basic geometry
		//fluid.replenish_boundary = true;

		// grav--unscaled
		fluid.g = 9.8 * VectorD::Unit(0);

		//external force
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;

		//main params
		fluid.default_sim_params = DefaultSimParams(dx, rho, fluid.viscosity_water, fluid.gamma_water, fluid.thickness, fluid.simulation_scale);

		//These three must go together
		fluid.gravity_coeff = 0.33 / R * 0.6 * 10 * pow(.1 / 0.15, 2) * fluid.default_sim_params.gravity_coeff;
		fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, 50);
		fluid.viscosity_coeff = 2 * fluid.default_sim_params.viscosity_coeff;
		//

		fluid.normal_viscosity_coeff = fluid.default_sim_params.normal_viscosity_coeff;

		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(0.33, fluid.gamma_water, rho, fluid.thickness);

		fluid.marangoni_coeff = fluid.default_sim_params.marangoni_coeff;

		//fluid.n_pressure_params.capillary_coeff *= R / 0.15; /// to adapt to the reduced G

		//operator parameters
		fluid.grad_force_params = fluid.default_sim_params.grad_force_params;
		fluid.height_laplacian_params = fluid.default_sim_params.height_laplacian_params;
		fluid.geometry_params = fluid.default_sim_params.geometry_params;

		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
		fluid.boundary_params.replenish = false;

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0. * dx));

		fluid.rh_params = RenderHeightParams("divergence", false, true, 0);

		fluid.Initialize(dx);
		fluid.clip_velocity = true;
		fluid.vel_threshold = 0.1; //50 frames per second, one frame max 

		if (first_frame < 1 && init_snapshot_name != "none") {
			bool loaded;
			loaded = fluid.particles.Load_Snapshot(init_snapshot_name);
			if (loaded) {
				int n = fluid.particles.Size();
				for (int i = 0; i < n; i++) {
					fluid.particles.RH(i) = fluid.thickness;
					fluid.particles.RH_V(i) = 0.;
					fluid.particles.Div(i) = 0.;
					fluid.particles.V(i) *= 0.;
				}
			}
			else {
				std::cout << "Unable to Load From Snapshot!" << std::endl;
				exit(0);
			}
		}
		fluid.vorticity_params = VorticityParams(KernelType::QUINTIC, 1, true, 1e-4);
		Seed_Vortex_Rand(5, 0.33/0.15 * 0.00001 * Rand_Number(), 0.33/3.);
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_44(void)
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "Not serving d = 2, sorry!" << std::endl;
	}
	if constexpr (d == 3) {
		max_iter_per_frame = 10;
		cfl = 0.1;
		fluid.max_vel = 0.1;

		real R = (real)1.; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings

		fluid.simulation_scale = R;
		frame_rate = 10.;
		real fineness = 1;
		real dx = 0.005 / fineness; //this can be changed
		fluid.np_on_h = 6 * fineness;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";

		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			fluid.particles.RH(i) = fluid.thickness;
		}

		//basic geometry

		// grav--unscaled
		fluid.g = 9.8 * VectorD::Unit(0);

		//external force
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;

		//main params
		fluid.default_sim_params = DefaultSimParams(dx, rho, fluid.viscosity_water, fluid.gamma_water, fluid.thickness, fluid.simulation_scale);

		std::function<VectorD(const int)> ext_force = [&](const int i)->VectorD {
			const real g_strength = 5 * 9.8 * 0.33 / fluid.simulation_scale * 0.6 * 10 * pow(.1 / 0.15, 2) * fluid.default_sim_params.gravity_coeff;
			if (fluid.particles.Phase(i) > 0) {
				return VectorD::Unit(0) * g_strength * fluid.particles.M(i);
			}
			else {
				return VectorD::Zero();
			}
		};

		fluid.external_force_func = ext_force;

		//These three must go together
		fluid.gravity_coeff = 0.33 / R * 0.6 * 10 * pow(.1 / 0.15, 2) * fluid.default_sim_params.gravity_coeff;
		fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, 50);
		fluid.t_pressure_params.divergence_pressure_coeff *= 0.33;
		fluid.viscosity_coeff = 2 * fluid.default_sim_params.viscosity_coeff;
		//

		fluid.normal_viscosity_coeff = fluid.default_sim_params.normal_viscosity_coeff;

		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(0.33, fluid.gamma_water, rho, fluid.thickness);

		fluid.marangoni_coeff = fluid.default_sim_params.marangoni_coeff;

		//fluid.n_pressure_params.capillary_coeff *= R / 0.15; /// to adapt to the reduced G

		//operator parameters
		fluid.grad_force_params = fluid.default_sim_params.grad_force_params;
		fluid.height_laplacian_params = fluid.default_sim_params.height_laplacian_params;
		fluid.geometry_params = fluid.default_sim_params.geometry_params;

		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
		fluid.boundary_params.replenish = false;
		//fluid.boundary_params.replenish = true;
		//fluid.boundary_params.inherit_rh = true;
		//fluid.boundary_params.replenish_interval = 1e-3;
		//fluid.boundary_params.replenish_dx_num = 3.;

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0. * dx));

		fluid.rh_params = RenderHeightParams("divergence", false, true, 0);

		fluid.Initialize(dx);
		fluid.clip_velocity = true;
		fluid.vel_threshold = 0.1; //50 frames per second, one frame max 

		if (first_frame < 1 && init_snapshot_name != "none") {
			bool loaded;
			loaded = fluid.particles.Load_Snapshot(init_snapshot_name);
			if (loaded) {
				n = fluid.particles.Size();
				for (int i = 0; i < n; i++) {
					//fluid.particles.RH(i) = fluid.thickness;
					//fluid.particles.RH_V(i) = 0.;
					//fluid.particles.Div(i) = 0.;
					//fluid.particles.V(i) *= 0.;
					fluid.particles.Phase(i) = 0;
				}
			}
			else {
				std::cout << "Unable to Load From Snapshot!" << std::endl;
				exit(0);
			}
		}
		int num_nozzles = 1;
		Array<VectorD> drip_centers;
		real nozzle_size = 1.2 * fluid.length_scale;
		for (int k = 0;k < num_nozzles;k++) {
			real tmp = Rand_Number() * pi / 2. + ((3. * pi) / 4.);
			VectorD drip_center;
			drip_center << fluid.simulation_scale * cos(tmp), 0, fluid.simulation_scale* sin(tmp);
			drip_centers.push_back(drip_center);
		}
		for (int i = 0; i < fluid.particles.Size(); i++) {
			bool already_in = false;
			if (!fluid.particles.Is_Boundary(i))continue;
			for (int k = 0; k < num_nozzles;k++) {
				if ((fluid.particles.X(i) - drip_centers[k]).norm() < nozzle_size) {
					fluid.droplet_nozzles.push_back(i);
					fluid.drops.push_back(0);
					already_in = true;
				}
			}
			if (already_in)continue;
		}
		
		fluid.max_droplets = 18;
		fluid.use_multiphase = true;
		fluid.droplet_m = 30. * total_mass / n;
		fluid.droplet_rh = 10. * fluid.thickness;

		fluid.seeding_droplet = false;

		if (first_frame < 1) {
			for (int k = 0; k < fluid.max_droplets; k++) {
				//real tmp = (real)k / (real)fluid.max_droplets * (pi - 2 * 1. / 6. * pi) + (1./2.*pi + 1./6. * pi);
				real tmp = Rand_Number() * (pi - 2 * 1. / 6. * pi) + (1. / 2. * pi + 1. / 6. * pi);
				VectorD drip_center;
				real drip_r = (0.7 + Rand_Number() * 0.6) * (fluid.simulation_scale/26.) * (0.8 + 0.0 * Rand_Number());
				drip_center << (fluid.simulation_scale - 0.5 * drip_r) * cos(tmp), 0, (fluid.simulation_scale - 0.5 * drip_r)* sin(tmp);
				for (int i = 0; i < fluid.particles.Size(); i++) {
					real dist = (fluid.particles.X(i) - drip_center).norm();
					if (dist < drip_r) {
						fluid.particles.M(i) = fluid.droplet_m;
						fluid.particles.V(i) = 0.1 * VectorD::Unit(0);
						//particles.RH(i) = thickness + droplet_rh * std::min(cos(dist / drip_r * pi / 2.0), 1.0);
						fluid.particles.RH(i) = fluid.droplet_rh;
						fluid.particles.Phase(i) = 1;
					}
				}
			}
		}
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_45(void)
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "Not serving d = 2, sorry!" << std::endl;
	}
	if constexpr (d == 3) {
		max_iter_per_frame = -1;
		cfl = 0.1;
		fluid.max_vel = 0.1;

		real R = (real)1.; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings

		fluid.simulation_scale = R;
		frame_rate = 50.;
		real fineness = 1;
		real dx = 0.005 / fineness; //this can be changed
		fluid.np_on_h = 6 * fineness;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";

		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			fluid.particles.RH(i) = fluid.thickness;
		}

		//basic geometry
		//fluid.replenish_boundary = true;

		// grav--unscaled
		fluid.g = 0*9.8 * VectorD::Unit(0);

		//external force
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;

		//main params
		fluid.default_sim_params = DefaultSimParams(dx, rho, fluid.viscosity_water, fluid.gamma_water, fluid.thickness, fluid.simulation_scale);

		//These three must go together
		fluid.gravity_coeff = 0.33 / R * 0.6 * 10 * pow(.1 / 0.15, 2) * fluid.default_sim_params.gravity_coeff;
		fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, 50);
		fluid.viscosity_coeff = 2 * fluid.default_sim_params.viscosity_coeff;
		//

		fluid.normal_viscosity_coeff = fluid.default_sim_params.normal_viscosity_coeff;

		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(0.33, fluid.gamma_water, rho, fluid.thickness);

		fluid.marangoni_coeff = 0.0001;
		fluid.diffuse_soap = true;

		//fluid.n_pressure_params.capillary_coeff *= R / 0.15; /// to adapt to the reduced G

		//operator parameters
		fluid.grad_force_params = fluid.default_sim_params.grad_force_params;
		fluid.height_laplacian_params = fluid.default_sim_params.height_laplacian_params;
		fluid.geometry_params = fluid.default_sim_params.geometry_params;

		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
		fluid.boundary_params.replenish = false;

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0. * dx));

		fluid.rh_params = RenderHeightParams("divergence", false, true, 0);

		fluid.Initialize(dx);
		fluid.clip_velocity = true;
		fluid.vel_threshold = 0.1; //50 frames per second, one frame max 

		if (first_frame < 1 && init_snapshot_name != "none") {
			bool loaded;
			loaded = fluid.particles.Load_Snapshot(init_snapshot_name);
			if (loaded) {
				n = fluid.particles.Size();
				//for (int i = 0; i < n; i++) {
				//	fluid.particles.RH(i) = fluid.thickness;
				//	fluid.particles.RH_V(i) = 0.;
				//	fluid.particles.Div(i) = 0.;
				//	fluid.particles.V(i) *= 0.;
				//}
			}
			else {
				std::cout << "Unable to Load From Snapshot!" << std::endl;
				exit(0);
			}
		}

		
		for (int i = 0; i < n; i++) {
			if (fluid.particles.X(i).norm() < 0.1 * R) {
				fluid.particles.Conc(i) = 1.;
			}
		}
	}
}


template<int d>
void FluidSPHBubbleDriver<d>::Case_46(void)
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "Not serving d = 2, sorry!" << std::endl;
	}
	if constexpr (d == 3) {
		max_iter_per_frame = 20;
		cfl = 0.1;
		fluid.max_vel = 0.1;

		real R = (real).22; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings

		fluid.simulation_scale = R;
		frame_rate = 300.;
		real fineness = 1;
		real dx = 0.005 / fineness; //this can be changed
		fluid.np_on_h = 6 * fineness;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Catenoid_Points(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Catenoid_Points_Random(-R / 3. * VectorD::Unit(1), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		real avg_height = 1.61987e-06;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";

		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			if (is_boundary[i] == 0) {
				fluid.particles.B(i) = 0;
			}
			else{
				fluid.particles.B(i) = 1;
			}
			//fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			fluid.particles.RH(i) = avg_height;
		}

		//basic geometry
		//fluid.replenish_boundary = true;

		// grav--unscaled
		fluid.g = 0 * 9.8 * VectorD::Unit(0);

		//external force
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;

		//main params
		fluid.default_sim_params = DefaultSimParams(dx, rho, fluid.viscosity_water, fluid.gamma_water, fluid.thickness, fluid.simulation_scale);

		//These three must go together
		fluid.gravity_coeff = 0.33 / R * 0.6 * 10 * pow(.1 / 0.15, 2) * fluid.default_sim_params.gravity_coeff;
		fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, 50);
		fluid.viscosity_coeff = 2 * fluid.default_sim_params.viscosity_coeff;
		
		//fluid.t_pressure_params.tangential_pressure_coeff *= 6.66;
		//fluid.t_pressure_params.boundary_pressure_coeff *= 0;
	    fluid.t_pressure_params.height_pressure_coeff *= 150;
		//fluid.t_pressure_params.divergence_pressure_coeff *= 0;
		//fluid.t_pressure_params.laplacian_pressure_coeff *= 0;


		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(0.33, fluid.gamma_water, rho, fluid.thickness);
		
		//fluid.n_pressure_params.capillary_coeff *= 0.1/0.3 * pow(R/0.1, 2);
		//fluid.n_pressure_params.capillary_coeff *= 0.2 * (R/0.1); // this is what used to generate previous example
		fluid.n_pressure_params.capillary_coeff *= 0.2 * 30 * (R / 0.1);

		//fluid.n_pressure_params.capillary_coeff *= R / 0.15; /// to adapt to the reduced G

		//operator parameters
		fluid.grad_force_params = fluid.default_sim_params.grad_force_params;
		fluid.height_laplacian_params = fluid.default_sim_params.height_laplacian_params;
		fluid.geometry_params = fluid.default_sim_params.geometry_params;

		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
		fluid.boundary_params.replenish = true;
		fluid.boundary_params.keep_xz_plane = false;

		//fluid.boundary_params.boundary_force_mode = "binary";
		//fluid.boundary_params.grav_boundary = 0.1;
		
		//fluid.boundary_params.boundary_mode = "none";

		fluid.normal_viscosity_coeff = 0 * fluid.default_sim_params.normal_viscosity_coeff;

		fluid.exp_mode = "catenoid"; fluid.catenoid_speed = 1./3. * R;


		fluid.catenoid_stop_frame = 500.;
		fluid.catenoid_stop_rate = 1.;

		VectorD catenoid_positive_center = VectorD::Unit(1)*R/3.;
		for (int i = 0; i < first_frame; i++) {
			real temp_speed = fluid.catenoid_speed;
			if (i > fluid.catenoid_stop_frame) temp_speed = fluid.catenoid_speed * pow(fluid.catenoid_stop_rate, i - fluid.catenoid_stop_frame);
			catenoid_positive_center += temp_speed * VectorD::Unit(1) * 1. / frame_rate;
		}
		//VectorD catenoid_positive_center = (R / 3. + fluid.catenoid_speed * Time_At_Frame(first_frame))* VectorD::Unit(1);
		fluid.analytical_boundary.Add_Obstacle(std::make_shared<HalfBowl<d>>(-catenoid_positive_center, R + 0. * dx, -VectorD::Unit(1)));
		fluid.analytical_boundary.Add_Obstacle(std::make_shared<HalfBowl<d>>(catenoid_positive_center, R + 0* dx, VectorD::Unit(1)));

		fluid.rh_params = RenderHeightParams("divergence", false, true, 0);

		bool use_3d = true;
		if (use_3d) {
			fluid_3d.Initialize(dx, rho, R, fluid.gravity_coeff* fluid.g);
			fluid_3d.use_surface_tension = true;
			fluid_3d.prune_far_away = true;
			fluid_3d.nden_0 = 3e6;
			fluid_3d.pressure_multiplier = 1*5e4;
			//fluid_3d.curvature_multiplier = 0.1;
			fluid_3d.curvature_multiplier = 1*15.;
			fluid_3d.cohesion_multiplier = 1*30.;
			fluid_3d.viscosity_multiplier = 3.;
			fluid.Initialize(dx, nullptr, &fluid_3d);
			//fluid_3d.Initialize(dx, rho, R, fluid.gravity_coeff * fluid.g);
			//fluid_3d.use_surface_tension = true;
			//fluid_3d.prune_far_away = true;
			//fluid_3d.nden_0 = 3e6;
			//fluid_3d.pressure_multiplier = 5e1;
			////fluid_3d.curvature_multiplier = 0.1;
			//fluid_3d.curvature_multiplier = .33;
			//fluid_3d.cohesion_multiplier = .33;
			//fluid_3d.viscosity_multiplier = .5;
			//fluid.Initialize(dx, nullptr, &fluid_3d);
		}
		else{ fluid.Initialize(dx); }

		fluid.clip_velocity = true;
		fluid.vel_threshold = 0.1; //50 frames per second, one frame max 

		fluid.delete_speedy_particles = true;
		fluid.quit_when_too_fast = false;

		fluid.rh_params.blend_h = true;
		fluid.rh_params.blend_constant = 0.002;
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_47(void)
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "Not serving d = 2, sorry!" << std::endl;
	}
	if constexpr (d == 3) {
		max_iter_per_frame = 20;
		cfl = 0.1;
		fluid.max_vel = 0.1;

		real R = (real).15; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings

		fluid.simulation_scale = R;
		frame_rate = 50.;
		real fineness = 1;
		real dx = 0.005 / fineness; //this can be changed
		fluid.np_on_h = 6 * fineness;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Catenoid_Points(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Catenoid_Points_Random(-R / 3. * VectorD::Unit(1), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";

		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			if (is_boundary[i] == 0) {
				fluid.particles.B(i) = 0;
			}
			else if (is_boundary[i] == 1) {
				fluid.particles.B(i) = 1;
			}
			else if (is_boundary[i] == 2) {
				fluid.particles.B(i) = 1;
			}
			//fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			fluid.particles.RH(i) = fluid.thickness;
		}

		//basic geometry
		//fluid.replenish_boundary = true;

		// grav--unscaled
		fluid.g = 0 * 9.8 * VectorD::Unit(0);

		//external force
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;

		//main params
		fluid.default_sim_params = DefaultSimParams(dx, rho, fluid.viscosity_water, fluid.gamma_water, fluid.thickness, fluid.simulation_scale);

		//These three must go together
		fluid.gravity_coeff = 0.33 / R * 0.6 * 10 * pow(.1 / 0.15, 2) * fluid.default_sim_params.gravity_coeff;
		fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, 50);
		fluid.viscosity_coeff = 2 * fluid.default_sim_params.viscosity_coeff;

		//fluid.t_pressure_params.tangential_pressure_coeff *= 6.66;
		//fluid.t_pressure_params.boundary_pressure_coeff *= 0;
		fluid.t_pressure_params.height_pressure_coeff *= 6.66;
		//fluid.t_pressure_params.divergence_pressure_coeff *= 0;
		//fluid.t_pressure_params.laplacian_pressure_coeff *= 0;


		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(0.33, fluid.gamma_water, rho, fluid.thickness);

		fluid.n_pressure_params.capillary_coeff *= 0.1 / 0.3 * pow(R / 0.1, 2);

		//fluid.n_pressure_params.capillary_coeff *= R / 0.15; /// to adapt to the reduced G

		//operator parameters
		fluid.grad_force_params = fluid.default_sim_params.grad_force_params;
		fluid.height_laplacian_params = fluid.default_sim_params.height_laplacian_params;
		fluid.geometry_params = fluid.default_sim_params.geometry_params;

		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
		fluid.boundary_params.replenish = false;
		fluid.boundary_params.boundary_force_mode = "binary";
		fluid.boundary_params.grav_boundary = 1;
		fluid.boundary_params.keep_xz_plane = false;

		//fluid.boundary_params.boundary_mode = "none";

		fluid.normal_viscosity_coeff = 0 * fluid.default_sim_params.normal_viscosity_coeff;

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<HalfBowl<d>>(-R / 3. * VectorD::Unit(1), R + 0. * dx, -VectorD::Unit(1)));
		fluid.analytical_boundary.Add_Obstacle(std::make_shared<HalfBowl<d>>(R / 3. * VectorD::Unit(1), R + 0 * dx, VectorD::Unit(1)));

		fluid.rh_params = RenderHeightParams("divergence", false, true, 0);

		fluid_3d.Initialize(dx);
		//fluid_3d.use_central_gravity = true;
		///////////////

		fluid.Initialize(dx, nullptr, &fluid_3d);
		fluid.clip_velocity = true;
		fluid.vel_threshold = 0.1; //50 frames per second, one frame max 

		fluid.exp_mode = "catenoid";
		fluid.delete_speedy_particles = true;
		snapshot_stride = 5;
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_48(void)
{


	std::cout << "Enter initialization Case 48\n";

	//int sphere_sub_num = 5;//5:10242, 6:40962, 7: 163842
	//real R = 0.133;//5:0.133, 7:0.5320
	int sphere_sub_num = 4;//5:10242, 6:40962, 7: 163842
	real R = 0.133/2;//5:0.133, 7:0.5320
	//int sphere_sub_num = 6;//5:10242, 6:40962, 7: 163842
	//real R = 0.133 * 2;//5:0.133, 7:0.5320
	//int sphere_sub_num = 7;//5:10242, 6:40962, 7: 163842
	//real R = 0.532;//5:0.133, 7:0.5320

	max_iter_per_frame = -1;
	cfl = 0.1;
	frame_rate = 50;
	real dx = 0.005;
	int number_2d = 1000;//number of particles for 2d
	VectorD ctr = VectorD::Zero();
	if constexpr (d == 2) {
		R = dx * number_2d / (2 * pi);
		dx = PointSetFunc::Initialize_Circle_Points(ctr, R, number_2d, fluid.particles);
	}
	if constexpr (d == 3) { dx = PointSetFunc::Initialize_Sphere_Points(ctr, R, sphere_sub_num, fluid.particles); }
	fluid.simulation_scale = R;
	std::cout << "Initialize with " << fluid.particles.Size() << " particles and dx = " << dx << " and R = " << R << "\n";


	Set_Physical_Parameters();
	real rho = 1e3; fluid.thickness = 5e-7;//10um
	real default_conc = 1;//"extreme soap"
	real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
	real total_vol = total_surface * fluid.thickness, total_mass = total_vol * rho;
	int n = fluid.particles.Size();
	std::cout << "Particle Vol: " << total_vol / n << ", Particle Mass: " << total_mass / n << "\n";
	for (int i = 0; i < n; i++) {
		Initialize_Particle(i, total_mass / n, fluid.thickness);
		fluid.particles.B(i) = 0;
		fluid.particles.Vol(i) = total_vol / n;
		fluid.particles.Conc(i) = default_conc;
	}

	real V = (d == 2) ? R * R * pi : 4.0 / 3.0 * pi * R * R * R;
	std::cout << "Analytical volume: " << V << "\n";
	real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
	real alpha_vel = 0.005 * 50;//move 0.005 in one frame with FPS=50

	//world forces
	fluid.gravity_coeff = 0;
	//fluid.friction_coeff = 1, fluid.air_velocity_func = Corridor_Flow_Func(R, alpha_vel);
	//fluid.external_force_func = ext_force;

	//tangential force
	fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, frame_rate);
	fluid.t_pressure_params.height_pressure_coeff *= 0;

	fluid.rh_params = RenderHeightParams("laplacian", true, false, 0);

	fluid.normal_viscosity_coeff = 1 * 0;
	fluid.boundary_params.Set_No_Boundary();

	fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
	fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);

	//fin.close();
	fluid.max_vel = alpha_vel;
	bool use_3d = true;
	if (use_3d) {
		fluid_3d.Initialize(0.2*dx, rho, R, fluid.gravity_coeff * fluid.g);
		fluid_3d.use_surface_tension = true;
		fluid_3d.nden_0 = 3e6;
		fluid_3d.pressure_multiplier = 0.01*5e2;
		fluid_3d.curvature_multiplier = 0.2 * 0.1;
		fluid_3d.cohesion_multiplier = 0.2 * 1.;
		fluid_3d.viscosity_multiplier = 0.2 * .5;
		fluid.Initialize(dx, nullptr, &fluid_3d);
	}
	else { fluid.Initialize(dx); }
	//fluid.Initialize(dx);

	for (int i = 0; i < fluid.particles.Size(); i++) fluid.particles.RH(i) = fluid.particles.H(i);

	real V_num = fluid.Compute_Enclosed_Volume();
	fluid.n_pressure_params.Set_Sphere_Baseline(d, R, fluid.Surface_Tension_Coefficient(default_conc), V_num, rho, fluid.thickness);

	if (first_frame < 1 && init_snapshot_name != "none") {
		bool loaded;
		loaded = fluid.particles.Load_Snapshot(init_snapshot_name);
		if (!loaded) {
			std::cout << "Unable to Load From Snapshot!" << std::endl;
			exit(0);
		}
	}
	
	int rightmost_idx=0;
	real rightmost = -1.;
	for (int i = 0; i < fluid.particles.Size(); i++) {
		if (fluid.particles.X(i)[0] > rightmost) {
			rightmost = fluid.particles.X(i)[0];
			rightmost_idx = i;
		}
	}


	if (first_frame == 0) {
		for (int i = 0; i < fluid.particles.Size(); i++) {
			if ((fluid.particles.X(i)-fluid.particles.X(rightmost_idx)).norm() < 0.15 * R) {
				fluid.to_explode.push_back(i);
				//fluid.particles.V(i) = -.2 * VectorD::Unit(0);
			}
			//fluid.to_explode.push_back(i);
		}
		std::cout << "num to explode: " << fluid.to_explode.size() << std::endl;
		fluid.init_explosion_time = 5. / 50.;
		fluid.explode_time = fluid.init_explosion_time + 5. / 50.;
	}
	else {
		fluid.to_explode.clear();
	}

	fluid.n_pressure_params.capillary_coeff *= 1;
	fluid.n_pressure_params.init_capillary_coeff = fluid.n_pressure_params.capillary_coeff;
	//fluid.n_pressure_params.capillary_coeff *= 1 * 0.0002;
	//fluid.n_pressure_params.capillary_coeff *= 0.3;
	fluid.n_pressure_params.closed = true;
	//fluid.t_pressure_params.height_pressure_coeff *= -3.;
	//fluid.viscosity_coeff = 0 * fluid.default_sim_params.viscosity_coeff;
	fluid.exp_mode = "bursting_bubble";
	fluid.clip_velocity = true;
	fluid.vel_threshold = 0.1; //50 frames per second, one frame max 
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_49(void) //
{
	std::cout << "Enter initialization Case 49\n";

	max_iter_per_frame = -1;
	cfl = 0.1;
	frame_rate = 50;
	real dx = 0.005;
	if constexpr (d == 3) { PointSetFunc::Initialize_Box_Points(7, 20, 7, dx, VectorD::Zero(), fluid_3d.particles); }
	real R = 0.1;
	//if constexpr (d == 3) dx = PointSetFunc::Initialize_Sphere_Points(VectorD::Zero(), R, 3, fluid_3d.particles);

	real rho = 1e3;
	fluid.simulation_scale = dx * 10;
	fluid.particles.Add_Element();
	Initialize_Particle(0, 1.0, 1.0);
	fluid.particles.E(0) = MatrixD::Identity();
	fluid.particles.X(0) = AuxFunc::V<d>(0, 1e4, 0);
	Set_Physical_Parameters();

	//static Array<bool> to_delete; to_delete.resize(fluid_3d.particles.Size()); AuxFunc::Fill(to_delete, false);
	for (int i = 0; i < fluid_3d.particles.Size(); i++) {
		fluid_3d.particles.M(i) = 8e-8;
		for (int axis = 0; axis < d; axis++) {
			real offset = (Rand_Number() - 0.5) * dx * 0.5;
			fluid_3d.particles.X(i)[axis] += offset;
		}
		//if (fluid_3d.particles.X(i)[0] > 2 * R) to_delete[i] = true;
	}
	//fluid_3d.particles.Delete_Elements(to_delete);

	//world forces
	fluid.gravity_coeff = 0;// 0.02;

	fluid.max_vel = 0.005 * 50;//move 0.005 in one frame with FPS=50
	bool use_3d = true;
	if (use_3d) {
		//fluid_3d.Initialize(dx, rho, R, fluid.gravity_coeff * fluid.g);
		//fluid_3d.use_surface_tension = true;
		//fluid_3d.nden_0 = 3e6;
		//fluid_3d.pressure_multiplier = 5e2;
		//fluid_3d.curvature_multiplier = 0.1;
		//fluid_3d.cohesion_multiplier = 0.1;
		//fluid_3d.viscosity_multiplier = .5;
		//fluid.Initialize(dx, nullptr, &fluid_3d);
		fluid_3d.Initialize(dx, rho, R, fluid.gravity_coeff * fluid.g);
		fluid_3d.use_surface_tension = true;
		fluid_3d.nden_0 = 3e6;
		fluid_3d.pressure_multiplier = 10e2;
		fluid_3d.curvature_multiplier = .3;
		fluid_3d.cohesion_multiplier = .4;
		fluid_3d.viscosity_multiplier = .2;
		fluid.Initialize(dx, nullptr, &fluid_3d);
	}


	//fluid_3d.analytical_boundary.Add_Obstacle(std::make_shared<Plane<d>>(VectorD::Unit(1), -VectorD::Unit(1) * 6 * dx));
	//if constexpr (d == 3){
	//	SPHBubbleParticles<3> boundary_particles;
	//	Vector2i cell_counts = AuxFunc::Vi<2>(30, 30);
	//	PointSetFunc::Initialize_Lattice_Points(-VectorD::Unit(1) * (6 * dx - 0.5*dx), cell_counts, Vector3::Unit(0), Vector3::Unit(2), dx, boundary_particles);
	//	for (int i = 0; i < boundary_particles.Size();i++) {
	//		int k = fluid_3d.particles.Add_Element();
	//		fluid_3d.particles.X(k) = boundary_particles.X(i);
	//		fluid_3d.particles.M(k) = 3e-8;
	//		fluid_3d.particles.B(k) = 1;
	//	}
	//}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_50(void) //play with different init density
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "Not serving d = 2, sorry!" << std::endl;
	}
	if constexpr (d == 3) {
		max_iter_per_frame = -1;
		cfl = 0.1;
		fluid.max_vel = 0.1;

		real R = (real).22; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings

		fluid.simulation_scale = R;
		frame_rate = 50.;
		real fineness = 1;
		real dx = 0.005 / fineness; //this can be changed
		fluid.np_on_h = 6 * fineness;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;
		real actual_thickness = 1.87e-6;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";

		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			fluid.particles.RH(i) = actual_thickness;
		}

		//basic geometry
		//fluid.replenish_boundary = true;

		// grav--unscaled
		fluid.g = 1*9.8 * -VectorD::Unit(1);

		//external force
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;

		//main params
		fluid.default_sim_params = DefaultSimParams(dx, rho, fluid.viscosity_water, fluid.gamma_water, fluid.thickness, fluid.simulation_scale);

		//These three must go together
		fluid.gravity_coeff = 0.33 / .1 * 0.6 * 10 * pow(.1 / 0.15, 2) * fluid.default_sim_params.gravity_coeff;
		fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, 50);
		fluid.viscosity_coeff = 2 * fluid.default_sim_params.viscosity_coeff;
		//

		fluid.normal_viscosity_coeff = fluid.default_sim_params.normal_viscosity_coeff;

		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(0.1, fluid.gamma_water, rho, fluid.thickness);
		fluid.n_pressure_params.capillary_coeff *= 2.5;

		fluid.marangoni_coeff = fluid.default_sim_params.marangoni_coeff;

		//fluid.n_pressure_params.capillary_coeff *= R / 0.15; /// to adapt to the reduced G

		//operator parameters
		fluid.grad_force_params = fluid.default_sim_params.grad_force_params;
		fluid.height_laplacian_params = fluid.default_sim_params.height_laplacian_params;
		fluid.geometry_params = fluid.default_sim_params.geometry_params;

		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
		fluid.boundary_params.replenish = false;

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0. * dx));

		fluid.rh_params = RenderHeightParams("divergence", false, true, 0);

		/*fluid.Initialize(dx);*/
		bool use_3d = true;

		if (use_3d) {
			//fluid_3d.Initialize(dx, rho, R, fluid.gravity_coeff * fluid.g);
			//fluid_3d.use_surface_tension = true;
			//fluid_3d.nden_0 = 3e6;
			//fluid_3d.pressure_multiplier = 5e2;
			//fluid_3d.curvature_multiplier = 0.1;
			//fluid_3d.cohesion_multiplier = 0.1;
			//fluid_3d.viscosity_multiplier = .5;
			//fluid.Initialize(dx, nullptr, &fluid_3d);
			fluid_3d.Initialize(dx, rho, R, fluid.gravity_coeff * fluid.g);
			fluid_3d.use_surface_tension = true;
			fluid_3d.nden_0 = 3e6;
			//fluid_3d.pressure_multiplier = 10e2;
			//fluid_3d.curvature_multiplier = .3;
			//fluid_3d.cohesion_multiplier = .4;
			//fluid_3d.viscosity_multiplier = .2;
			fluid_3d.pressure_multiplier = 10e2;
			//fluid_3d.curvature_multiplier = 0.1;
			fluid_3d.curvature_multiplier = 0.6;
			fluid_3d.cohesion_multiplier = 0.5;
			fluid_3d.viscosity_multiplier = 0.66;
			fluid.Initialize(dx, nullptr, &fluid_3d);


			PointSetFunc::Initialize_Box_Points(5, 20, 4, dx, 1.3*VectorD::Unit(1) * 0.6*R, fluid_3d.particles);
			for (int i = 0; i < fluid_3d.particles.Size(); i++) {
				fluid_3d.particles.M(i) = total_mass / n;
				fluid_3d.particles.B(i) = 0;
				fluid_3d.particles.Vol(i) = total_vol / n;
				fluid_3d.particles.Conc(i) = 1;
				fluid_3d.particles.RH(i) = actual_thickness;

			}
			SPHBubbleParticles<3> added_particles;
			//add second
			PointSetFunc::Initialize_Box_Points(18, 4, 6, dx, 1.3*VectorD::Unit(1) * 0.8 * R + VectorD::Unit(0) * 0.6* R, added_particles);
			for (int i = 0; i < added_particles.Size();i++) {
				int k = fluid_3d.particles.Add_Element();
				fluid_3d.particles.X(k) = added_particles.X(i);
				fluid_3d.particles.M(k) = total_mass / n;
				fluid_3d.particles.Vol(k) = total_vol / n;
				fluid_3d.particles.Conc(k) = 1;
				fluid_3d.particles.RH(k) = actual_thickness;
			}
			//add third
			PointSetFunc::Initialize_Box_Points(6, 17, 6, dx, 1.3*VectorD::Unit(1) * 0.9 * R + VectorD::Unit(0) * -0.6 * R + VectorD::Unit(2) * 0.4 * R, added_particles);
			for (int i = 0; i < added_particles.Size();i++) {
				int k = fluid_3d.particles.Add_Element();
				fluid_3d.particles.X(k) = added_particles.X(i);
				fluid_3d.particles.M(k) = total_mass / n;
				fluid_3d.particles.Vol(k) = total_vol / n;
				fluid_3d.particles.Conc(k) = 1;
				fluid_3d.particles.RH(k) = actual_thickness;
			}
			//add fourth
			PointSetFunc::Initialize_Box_Points(8, 8, 8, dx, 1.3*VectorD::Unit(1) * 1.2 * R + VectorD::Unit(0) * -0.1 * R + VectorD::Unit(2) * -0.7 * R, added_particles);
			for (int i = 0; i < added_particles.Size();i++) {
				int k = fluid_3d.particles.Add_Element();
				fluid_3d.particles.X(k) = added_particles.X(i);
				fluid_3d.particles.M(k) = total_mass / n;
				fluid_3d.particles.Vol(k) = total_vol / n;
				fluid_3d.particles.Conc(k) = 1;
				fluid_3d.particles.RH(k) = actual_thickness;
			}
		}
		else { fluid.Initialize(dx); }

		fluid.diffuse_soap = true;
		fluid.marangoni_coeff = 0.001;
		fluid.merge_3d = true;

		fluid.clip_velocity = true;
		fluid.vel_threshold = 0.1; //50 frames per second, one frame max 


		fluid.rh_params.blend_h = true;
		fluid.rh_params.blend_constant = 0.003;

		if (first_frame < 1 && init_snapshot_name != "none") {
			bool loaded;
			loaded = fluid.particles.Load_Snapshot(init_snapshot_name);
			if (loaded) {
				int n = fluid.particles.Size();
				for (int i = 0; i < n; i++) {
					fluid.particles.RH(i) = fluid.thickness;
					fluid.particles.RH_V(i) = 0.;
					fluid.particles.Div(i) = 0.;
					fluid.particles.V(i) *= 0.;
				}
			}
			else {
				std::cout << "Unable to Load From Snapshot!" << std::endl;
				exit(0);
			}
		}
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_51(void)
{
	if constexpr (d != 3) {
		std::cout << "case_51 error: d!=3" << std::endl;
	}
	if constexpr (d == 3) {

		real vel_term = 0.5, capillary_term = 0.5;
		//std::ifstream fin(scene_file_name);
		//if (!fin.is_open()) { std::cerr << "Driver Case_3 error: cannot open " << scene_file_name << "\n"; exit(0); }
		//fin >> vel_term >> capillary_term;
		//fin.close();

		max_iter_per_frame = -1;
		cfl = 1;
		frame_rate = 50;
		fluid.max_vel = 0.1;

		real R = (real)0.1;//0.1:1k, 0.3:11k, 1.0:126k
		fluid.simulation_scale = R;
		real dx = 0.005;
		Array<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();
		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		VectorD ctr = VectorD::Unit(0) * (-2 * R);
		MatrixD Rot; Rot << 0, 1, 0,
			0, 0, 1,
			1, 0, 0;
		int n = fluid.particles.Size();
		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			fluid.particles.RH(i) = fluid.thickness;
			fluid.particles.Apply_Rotation(i, Rot);
			fluid.particles.Apply_Translation(i, ctr);
		}

		real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;

		//world forces
		fluid.gravity_coeff = 0;
		//tangential force
		fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, frame_rate);
		fluid.t_pressure_params.height_pressure_coeff *= 5;
		fluid.divergence_params.calculate_field = "all";
		fluid.viscosity_coeff = 0.4 * alpha_vis;
		fluid.marangoni_coeff = 0;
		//render height
		fluid.rh_params = RenderHeightParams("divergence", true, true, 0.01);
		//normal
		fluid.normal_viscosity_coeff = 0;
		fluid.n_pressure_params.Set_IB(R, fluid.gamma_water, rho, fluid.thickness);
		fluid.n_pressure_params.capillary_coeff *= capillary_term;
		fluid.n_pressure_params.ib_force_coeff = 2;
		EulerInitializer<d> perimeter;
		real source_speed = vel_term;
		perimeter.Set_Domain(R, 10, AuxFunc::Vi<d>(6, 3, 3));
		perimeter.Set_Boundary_Width(1, -1, 1, 1, 1, 1);
		perimeter.Set_Boundary_Value(0, 0, 0, 0, 0, 0);
		perimeter.Generate_Parameters();
		//boundary
		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
		fluid.boundary_params.replenish = true;
		fluid.boundary_params.replenish_interval = 1e-3;
		fluid.boundary_params.replenish_dx_num = 1.5;
		fluid.boundary_params.keep_xz_plane = false;
		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Tube<d>>(ctr, R, VectorD::Unit(0), 0.5 * R));
		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Plane<d>>(VectorD::Unit(0), ctr - VectorD::Unit(0) * (0.5 * dx)));

		fluid.delete_solitary = false;

		fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
		fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);

		fluid.Initialize(dx, &perimeter);

		real nozzle_radius = R * 0.5;
		VectorD nozzle_ctr = AuxFunc::V<d>(-3, 0, 0) * R;
		Sphere<d> sphere(nozzle_ctr, nozzle_radius);
		for (auto p : fluid.air_solver.bc.psi_N_values) {
			int axis = p.first[0];
			int face_index = p.first[1];
			VectorDi face = fluid.air_solver.mac_grid.Face_Coord(axis, face_index);
			VectorD pos = fluid.air_solver.mac_grid.Face_Center(axis, face);
			if (axis == 0 && sphere.Inside(pos)) {
				fluid.air_solver.bc.Set_Psi_N(axis, face, source_speed);
			}
		}
		fluid.air_solver.kernel_coeff = 0.9;
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_52(void) //play with different init density
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "Not serving d = 2, sorry!" << std::endl;
	}
	if constexpr (d == 3) {
		max_iter_per_frame = -1;
		cfl = 0.1;
		fluid.max_vel = 0.1;

		real R = (real).2; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings

		fluid.simulation_scale = R;
		//frame_rate = 50.;
		frame_rate = 10.;
		real fineness = 1;
		real dx = 0.005 / fineness; //this can be changed
		fluid.np_on_h = 6 * fineness;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";

		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			fluid.particles.RH(i) = 1.93832e-06;
		}

		//basic geometry
		//fluid.replenish_boundary = true;

		// grav--unscaled
		fluid.g = 0.25 * 9.8 * VectorD::Unit(0);

		//external force
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;

		//main params
		fluid.default_sim_params = DefaultSimParams(dx, rho, fluid.viscosity_water, fluid.gamma_water, fluid.thickness, fluid.simulation_scale);

		//These three must go together
		fluid.gravity_coeff = 0.33 / R * 0.6 * 10 * pow(.1 / 0.15, 2) * fluid.default_sim_params.gravity_coeff;
		std::cout << "gravity constant: " << fluid.gravity_coeff * 0.25 * 9.8 << std::endl;
		fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, 50);
		fluid.viscosity_coeff = 4 * fluid.default_sim_params.viscosity_coeff;
		fluid.t_pressure_params.height_pressure_coeff *= 0.1;
		fluid.t_pressure_params.laplacian_pressure_coeff *= 2.5;
		//

		fluid.normal_viscosity_coeff = 0* fluid.default_sim_params.normal_viscosity_coeff;

		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(0.1, fluid.gamma_water, rho, fluid.thickness);
		fluid.n_pressure_params.capillary_coeff *= 0.;

		fluid.marangoni_coeff = fluid.default_sim_params.marangoni_coeff;

		//fluid.n_pressure_params.capillary_coeff *= R / 0.15; /// to adapt to the reduced G

		//operator parameters
		fluid.grad_force_params = fluid.default_sim_params.grad_force_params;
		fluid.height_laplacian_params = fluid.default_sim_params.height_laplacian_params;
		fluid.geometry_params = fluid.default_sim_params.geometry_params;

		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
		fluid.boundary_params.replenish = false;

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0. * dx));

		fluid.rh_params = RenderHeightParams("divergence", false, true, 0);

		fluid.Initialize(dx);
		fluid.clip_velocity = true;
		fluid.vel_threshold = 0.1; //50 frames per second, one frame max 

		//fluid.g_init = 9.8 * VectorD::Unit(0);
		//fluid.g_final = 9.8 * -VectorD::Unit(1);
		//fluid.interp_g = true;
		//fluid.first_frame = first_frame;
		//fluid.last_frame = last_frame;

		fluid.boundary_params.replenish = true;
		fluid.boundary_params.inherit_rh = false;
		fluid.boundary_params.replenish_interval = 1e-3;
		fluid.boundary_params.replenish_dx_num = 1.5;
		fluid.replenish_V_rate = 0.0;

		fluid.rh_params.blend_h = true;
		fluid.rh_params.blend_constant = 0.006;


		if (first_frame < 1 && init_snapshot_name != "none") {
			bool loaded;
			loaded = fluid.particles.Load_Snapshot(init_snapshot_name);
			if (loaded) {
				int n = fluid.particles.Size();
				//for (int i = 0; i < n; i++) {
				//	fluid.particles.RH(i) = fluid.thickness;
				//	fluid.particles.RH_V(i) = 0.;
				//	fluid.particles.Div(i) = 0.;
				//	fluid.particles.V(i) *= 0.;
				//}
			}
			else {
				std::cout << "Unable to Load From Snapshot!" << std::endl;
				exit(0);
			}
		}
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_53(void) //play with different init density
{
	std::function<VectorD(const int)> ext_force = nullptr;
	if constexpr (d == 2) {
		std::cout << "Not serving d = 2, sorry!" << std::endl;
	}
	if constexpr (d == 3) {
		max_iter_per_frame = -1;
		cfl = 0.1;
		fluid.max_vel = 0.1;

		real R = (real).2; //how radius of the rim centered at (0,0,0) --- Let's not give it any physical meanings

		fluid.simulation_scale = R;
		//frame_rate = 50.;
		frame_rate = 10.;
		real fineness = 1;
		real dx = 0.005 / fineness; //this can be changed
		fluid.np_on_h = 6 * fineness;

		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Random(VectorD::Zero(), R, dx, fluid.particles);
		//std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid(VectorD::Zero(), R, dx, fluid.particles);
		std::vector<int> is_boundary = PointSetFunc::Initialize_Circle_Points_Grid_Random(VectorD::Zero(), R, dx, fluid.particles);

		Set_Physical_Parameters();

		real rho = 1e3;
		fluid.thickness = 5e-7;//10um
		real default_conc = 0;//"extreme soap"
		real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
		real total_vol = total_surface * fluid.thickness;
		real total_mass = total_vol * rho;

		int n = fluid.particles.Size();
		std::cout << "single particle mass: " << total_mass / n << " vol: " << total_vol / n << "\n";

		for (int i = 0; i < n; i++) {
			Initialize_Particle(i, total_mass / n, fluid.thickness);
			fluid.particles.B(i) = is_boundary[i];
			fluid.particles.Vol(i) = total_vol / n;
			fluid.particles.Conc(i) = default_conc;
			fluid.particles.RH(i) = 1.93832e-06;
		}

		//basic geometry
		//fluid.replenish_boundary = true;

		// grav--unscaled
		fluid.g = 0. * 9.8 * VectorD::Unit(0);

		//external force
		fluid.friction_coeff = 0, fluid.air_velocity_func = nullptr;

		//main params
		fluid.default_sim_params = DefaultSimParams(dx, rho, fluid.viscosity_water, fluid.gamma_water, fluid.thickness, fluid.simulation_scale);

		//These three must go together
		fluid.gravity_coeff = 0.33 / R * 0.6 * 10 * pow(.1 / 0.15, 2) * fluid.default_sim_params.gravity_coeff;
		fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, 50);
		fluid.viscosity_coeff = .7 * fluid.default_sim_params.viscosity_coeff;
		fluid.t_pressure_params.divergence_pressure_coeff *= 0.;
		fluid.t_pressure_params.height_pressure_coeff *= .2;
		fluid.t_pressure_params.height_pressure_coeff *= .0;

		//fluid.t_pressure_params.laplacian_pressure_coeff *= 0;
		fluid.t_pressure_params.laplacian_pressure_coeff *= 3;
		//

		fluid.normal_viscosity_coeff = 0 * fluid.default_sim_params.normal_viscosity_coeff;

		//--normal direction
		fluid.n_pressure_params.Set_Circle_Baseline(0.1, fluid.gamma_water, rho, fluid.thickness);
		fluid.n_pressure_params.capillary_coeff *= 0.;

		fluid.marangoni_coeff = fluid.default_sim_params.marangoni_coeff;

		//fluid.n_pressure_params.capillary_coeff *= R / 0.15; /// to adapt to the reduced G

		//operator parameters
		fluid.grad_force_params = fluid.default_sim_params.grad_force_params;
		fluid.height_laplacian_params = fluid.default_sim_params.height_laplacian_params;
		fluid.geometry_params = fluid.default_sim_params.geometry_params;

		fluid.boundary_params.Set_Circle_Baseline(fluid.viscosity_coeff);
		fluid.boundary_params.replenish = false;

		fluid.analytical_boundary.Add_Obstacle(std::make_shared<Bowl<d>>(VectorD::Zero(), R + 0. * dx));

		fluid.rh_params = RenderHeightParams("divergence", false, true, 0);

		fluid.Initialize(dx);
		fluid.clip_velocity = true;
		fluid.vel_threshold = 0.1; //50 frames per second, one frame max 

		//fluid.g_init = 9.8 * VectorD::Unit(0);
		//fluid.g_final = 9.8 * -VectorD::Unit(1);
		//fluid.interp_g = true;
		//fluid.first_frame = first_frame;
		//fluid.last_frame = last_frame;

		//fluid.boundary_params.replenish = true;
		//fluid.boundary_params.inherit_rh = false;
		//fluid.boundary_params.replenish_interval = 1e-3;
		//fluid.boundary_params.replenish_dx_num = 1.5;
		//fluid.replenish_V_rate = 0.0;

		for (int i = 0; i < fluid.particles.Size(); i++) {
			if (fluid.particles.X(i).norm() < 0.7 * fluid.simulation_scale) {
				fluid.particles.V(i) = -.3 * fluid.particles.X(i);
			}
		}


		if (first_frame < 1 && init_snapshot_name != "none") {
			bool loaded;
			loaded = fluid.particles.Load_Snapshot(init_snapshot_name);
			if (loaded) {
				int n = fluid.particles.Size();
				//for (int i = 0; i < n; i++) {
				//	fluid.particles.RH(i) = fluid.thickness;
				//	fluid.particles.RH_V(i) = 0.;
				//	fluid.particles.Div(i) = 0.;
				//	fluid.particles.V(i) *= 0.;
				//}
			}
			else {
				std::cout << "Unable to Load From Snapshot!" << std::endl;
				exit(0);
			}
		}
	}
}

template<int d>
void FluidSPHBubbleDriver<d>::Case_54(void) //sphere curvature test
{


	std::cout << "Enter initialization Case 48\n";

	//int sphere_sub_num = 5;//5:10242, 6:40962, 7: 163842
	//real R = 0.133;//5:0.133, 7:0.5320
	int sphere_sub_num = 4;//5:10242, 6:40962, 7: 163842
	real R = 0.133/2.;//5:0.133, 7:0.5320
	//int sphere_sub_num = 6;//5:10242, 6:40962, 7: 163842
	//real R = 0.133 * 2;//5:0.133, 7:0.5320
	//int sphere_sub_num = 7;//5:10242, 6:40962, 7: 163842
	//real R = 0.532;//5:0.133, 7:0.5320

	max_iter_per_frame = -1;
	cfl = 0.1;
	frame_rate = 50;
	real dx = 0.005;
	int number_2d = 1000;//number of particles for 2d
	VectorD ctr = VectorD::Zero();
	if constexpr (d == 2) {
		R = dx * number_2d / (2 * pi);
		dx = PointSetFunc::Initialize_Circle_Points(ctr, R, number_2d, fluid.particles);
	}
	if constexpr (d == 3) { 
		SPHBubbleParticles<d> temp_ptc;
		dx = PointSetFunc::Initialize_Sphere_Points(ctr, R, sphere_sub_num, temp_ptc); 
		VectorD avg_pos = VectorD::Zero();
		for (int i = 0; i <temp_ptc.Size(); i++) {
			int k = fluid.particles.Add_Element();
			real theta = Rand_Number() * 2 * pi;
			real phi = Rand_Number() * 2 * pi;
			VectorD pos;
			pos << R * sin(theta) * cos(phi), R* sin(theta)* sin(phi), R* cos(theta);
			//std::cout << pos << std::endl;
			fluid.particles.X(k) = pos;
			fluid.particles.V(k) = Vector3::Zero();
			fluid.particles.F(k) = Vector3::Zero();
			fluid.particles.M(k) = (real)1;
			Vector3 normal, t1, t2;
			normal = fluid.particles.X(k).normalized();
			t1 = -AuxFunc::Orthogonal_Vector(normal);
			t2 = t1.cross(normal);
			fluid.particles.E(k).col(0) = t1;
			fluid.particles.E(k).col(1) = t2;
			fluid.particles.E(k).col(2) = normal;
			avg_pos += 1./ (real) temp_ptc.Size() * fluid.particles.X(k);
		}

		for (int i = 0; i < fluid.particles.Size();i++) {
			fluid.particles.X(i) -= avg_pos;
		}
	}
	fluid.simulation_scale = R;
	std::cout << "Initialize with " << fluid.particles.Size() << " particles and dx = " << dx << " and R = " << R << "\n";


	Set_Physical_Parameters();
	real rho = 1e3; fluid.thickness = 5e-7;//10um
	real default_conc = 1;//"extreme soap"
	real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
	real total_vol = total_surface * fluid.thickness, total_mass = total_vol * rho;
	int n = fluid.particles.Size();
	std::cout << "Particle Vol: " << total_vol / n << ", Particle Mass: " << total_mass / n << "\n";
	for (int i = 0; i < n; i++) {
		Initialize_Particle(i, total_mass / n, fluid.thickness);
		fluid.particles.B(i) = 0;
		fluid.particles.Vol(i) = total_vol / n;
		fluid.particles.Conc(i) = default_conc;
	}

	real V = (d == 2) ? R * R * pi : 4.0 / 3.0 * pi * R * R * R;
	std::cout << "Analytical volume: " << V << "\n";
	real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
	real alpha_vel = 0.005 * 50;//move 0.005 in one frame with FPS=50

	//world forces
	fluid.gravity_coeff = 0;
	//fluid.friction_coeff = 1, fluid.air_velocity_func = Corridor_Flow_Func(R, alpha_vel);
	//fluid.external_force_func = ext_force;

	//tangential force
	fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, frame_rate);
	fluid.t_pressure_params.tangential_pressure_coeff *= 0.;

	fluid.rh_params = RenderHeightParams("laplacian", true, false, 0);

	fluid.normal_viscosity_coeff = 1 * 0;
	fluid.boundary_params.Set_No_Boundary();

	fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
	fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);

	//fin.close();
	fluid.max_vel = alpha_vel;
	bool use_3d = false;
	if (use_3d) {
		fluid_3d.Initialize(0.2 * dx, rho, R, fluid.gravity_coeff * fluid.g);
		fluid_3d.use_surface_tension = true;
		fluid_3d.nden_0 = 3e6;
		fluid_3d.pressure_multiplier = 0.01 * 5e2;
		fluid_3d.curvature_multiplier = 0.2 * 0.1;
		fluid_3d.cohesion_multiplier = 0.2 * 1.;
		fluid_3d.viscosity_multiplier = 0.2 * .5;
		fluid.Initialize(dx, nullptr, &fluid_3d);
	}
	else { fluid.Initialize(dx); }
	//fluid.Initialize(dx);

	for (int i = 0; i < fluid.particles.Size(); i++) fluid.particles.RH(i) = fluid.particles.H(i);

	real V_num = fluid.Compute_Enclosed_Volume();
	fluid.n_pressure_params.Set_Sphere_Baseline(d, R, fluid.Surface_Tension_Coefficient(default_conc), V_num, rho, fluid.thickness);

	fluid.n_pressure_params.capillary_coeff *= 0;

}

template<int d>
void FluidSPHBubbleDriver<d>::Case_55(void) //analytical curvature test
{


	std::cout << "Enter initialization Case 48\n";

	//int sphere_sub_num = 5;//5:10242, 6:40962, 7: 163842
	//real R = 0.133;//5:0.133, 7:0.5320
	real R = 1.;//5:0.133, 7:0.5320
	real dx = 0.01;
	//int sphere_sub_num = 6;//5:10242, 6:40962, 7: 163842
	//real R = 0.133 * 2;//5:0.133, 7:0.5320
	//int sphere_sub_num = 7;//5:10242, 6:40962, 7: 163842
	//real R = 0.532;//5:0.133, 7:0.5320

	max_iter_per_frame = -1;
	cfl = 0.1;
	frame_rate = 50;

	if constexpr (d == 3) {
		SPHBubbleParticles<d> temp_ptc;
		dx = 1./100. * 2*pi * 1.8;
		VectorD avg_pos = VectorD::Zero();
		for (int i = 0; i < 10000; i++) {
			int k = fluid.particles.Add_Element();
			real x = Rand_Number() * 2 * pi;
			real y = Rand_Number() * 2 * pi;
			VectorD pos;
			pos << x, 0.1 * (3 * sin(x) + 2 * cos(y) + 4 * sin(2 * x + y)), y;
			//std::cout << pos << std::endl;
			fluid.particles.X(k) = pos;
			fluid.particles.V(k) = Vector3::Zero();
			fluid.particles.F(k) = Vector3::Zero();
			fluid.particles.M(k) = (real)1;
			fluid.particles.E(k).col(0) = -Vector3d::Unit(2);
			fluid.particles.E(k).col(1) = -Vector3d::Unit(0);
			fluid.particles.E(k).col(2) = -Vector3d::Unit(1);
			avg_pos += 1. / 10000. * fluid.particles.X(k);
		}

		for (int i = 0; i < fluid.particles.Size();i++) {
			fluid.particles.X(i) -= avg_pos;
		}
	}
	fluid.simulation_scale = R;
	std::cout << "Initialize with " << fluid.particles.Size() << " particles and dx = " << dx << " and R = " << R << "\n";


	Set_Physical_Parameters();
	real rho = 1e3; fluid.thickness = 5e-7;//10um
	real default_conc = 1;//"extreme soap"
	real total_surface = pi * pow(2 * R, 2);//works for both 2D and 3D
	real total_vol = total_surface * fluid.thickness, total_mass = total_vol * rho;
	int n = fluid.particles.Size();
	std::cout << "Particle Vol: " << total_vol / n << ", Particle Mass: " << total_mass / n << "\n";
	for (int i = 0; i < n; i++) {
		Initialize_Particle(i, total_mass / n, fluid.thickness);
		fluid.particles.B(i) = 0;
		fluid.particles.Vol(i) = total_vol / n;
		fluid.particles.Conc(i) = default_conc;
	}

	real V = (d == 2) ? R * R * pi : 4.0 / 3.0 * pi * R * R * R;
	std::cout << "Analytical volume: " << V << "\n";
	real alpha_vis = dx * dx * rho / (fluid.viscosity_water) * frame_rate * 2;
	real alpha_vel = 0.005 * 50;//move 0.005 in one frame with FPS=50

	//world forces
	fluid.gravity_coeff = 0;
	//fluid.friction_coeff = 1, fluid.air_velocity_func = Corridor_Flow_Func(R, alpha_vel);
	//fluid.external_force_func = ext_force;

	//tangential force
	fluid.t_pressure_params.Set_Baseline3(rho, fluid.gamma_water, fluid.thickness, frame_rate);
	fluid.t_pressure_params.tangential_pressure_coeff *= 0.;

	fluid.rh_params = RenderHeightParams("laplacian", true, false, 0);

	fluid.normal_viscosity_coeff = 1 * 0;
	fluid.boundary_params.Set_No_Boundary();

	fluid.grad_force_params = OperatorParams("only_fluid", KernelType::SPIKY);
	fluid.height_laplacian_params = OperatorParams("only_fluid", KernelType::GAUSSIAN);

	//fin.close();
	fluid.max_vel = alpha_vel;
	bool use_3d = false;
	if (use_3d) {
		fluid_3d.Initialize(0.2 * dx, rho, R, fluid.gravity_coeff * fluid.g);
		fluid_3d.use_surface_tension = true;
		fluid_3d.nden_0 = 3e6;
		fluid_3d.pressure_multiplier = 0.01 * 5e2;
		fluid_3d.curvature_multiplier = 0.2 * 0.1;
		fluid_3d.cohesion_multiplier = 0.2 * 1.;
		fluid_3d.viscosity_multiplier = 0.2 * .5;
		fluid.Initialize(dx, nullptr, &fluid_3d);
	}
	else { fluid.Initialize(dx); }
	//fluid.Initialize(dx);

	for (int i = 0; i < fluid.particles.Size(); i++) fluid.particles.RH(i) = fluid.particles.H(i);

	real V_num = fluid.Compute_Enclosed_Volume();
	fluid.n_pressure_params.Set_Sphere_Baseline(d, R, fluid.Surface_Tension_Coefficient(default_conc), V_num, rho, fluid.thickness);

	fluid.n_pressure_params.capillary_coeff *= 0;

}



template class FluidSPHBubbleDriver<2>;
template class FluidSPHBubbleDriver<3>;
