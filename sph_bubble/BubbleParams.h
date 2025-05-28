//////////////////////////////////////////////////////////////////////////
// Parameter classes for FluidSPHBubble
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __BubbleParams_h__
#define __BubbleParams_h__

#include <iostream>
#include "Common.h"
#include "Kernels.h"

//enum class KernelType { SPIKY, POLY6, CUBIC, GAUSSIAN, QUINTIC };

class OperatorParams {
public:
	std::string calculate_field = "only_fluid";//"only_fluid", "only_boundary", "all"
	KernelType kernel = KernelType::SPIKY;
	OperatorParams() {}
	OperatorParams(std::string _calculate) { calculate_field = _calculate; }
	OperatorParams(std::string _calculate, KernelType _type) { calculate_field = _calculate, kernel = _type; }
};

class VorticityParams {
public:
	KernelType weight_kernel = KernelType::QUINTIC;//its calculation field is "all" automatically
	real force_coeff = 0;
	bool enable_diffuse = false;
	real vorticity_diffuse_coeff = 0;
	VorticityParams() {}
	VorticityParams(KernelType _kernel, real _coeff, bool _diffuse, real _diffuse_coeff) {
		weight_kernel = _kernel;
		force_coeff = _coeff;
		enable_diffuse = _diffuse;
		vorticity_diffuse_coeff = _diffuse_coeff;
	}

	void Initialize_From_File(std::ifstream& fin) {
		//bool expressed in "0" or "1"
		//force_coeff, enable_diffuse, vorticity_diffuse_coeff, seed_vorticity
		if (!fin.is_open()) { std::cerr << "VorticityParams::Initialize_From_File error: ifstream not opened\n"; exit(0); }
		std::string token;
		while (std::getline(fin, token)) {
			if (token == "vorticity") {
				int tmp = 0;
				fin >> force_coeff;
				fin >> tmp; enable_diffuse = (bool)tmp;
				fin >> vorticity_diffuse_coeff;
				return;
			}
		}
		std::cerr << "VorticityParams::Initialize_From_File error: not token vorticity\n"; exit(0);
	}
};

class TangentialPressureParams {
public:
	real laplacian_pressure_coeff = 0.02;
	real divergence_pressure_coeff = 5;
	real height_pressure_coeff = 1;
	real tangential_pressure_coeff = 1e-1;
	real boundary_pressure_coeff = 0.;
	TangentialPressureParams() {}
private:
	real alpha_height = 0;
	real alpha_height_laplacian = 0;
	real alpha_divergence = 0;
public:
	void Calculate_Alpha(real _rho, real _gamma_water, real _thickness, int frame_rate) {
		alpha_height = 1.0 * _rho;
		alpha_height_laplacian = 1 * _rho / (2 * _gamma_water) / (_thickness / (pow((_thickness + 1), 2. / 3.)));
		alpha_divergence = _rho * frame_rate;
	}
	void Set_Baseline(real _rho, real _gamma_water, real _thickness, int frame_rate) {
		Calculate_Alpha(_rho, _gamma_water, _thickness, frame_rate);
		laplacian_pressure_coeff = 5 * 0.0000005 * alpha_height_laplacian;
		divergence_pressure_coeff = 1 * 0.0005 * alpha_divergence;
		height_pressure_coeff = 0. * alpha_height;
		tangential_pressure_coeff = 0.4;
		boundary_pressure_coeff = 0.1 / _thickness;
	}
	void Set_Baseline2(real _rho, real _gamma_water, real _thickness, int frame_rate, real R) {
		Calculate_Alpha(_rho, _gamma_water, _thickness, frame_rate);
		laplacian_pressure_coeff = R / 0.15 * 1 * 20000.;
		divergence_pressure_coeff = R / 0.15 * 50.;
		height_pressure_coeff = R / 0.15 * 1.5 * 10.;
		tangential_pressure_coeff = 0.4;
		boundary_pressure_coeff = 0.03 * 20 * 0.1 / _thickness;
	}
	void Set_Baseline3(real _rho, real _gamma_water, real _thickness, int frame_rate) {
		Calculate_Alpha(_rho, _gamma_water, _thickness, frame_rate);
		laplacian_pressure_coeff = 0.33 / 0.15 * 1 * 20000.;
		divergence_pressure_coeff = 0.6 * 0.33 / 0.15 * 50.;
		height_pressure_coeff = 0.25 * 0.33 / 0.15 * 1.5 * 10.;
		tangential_pressure_coeff = 0.4;
		boundary_pressure_coeff = 0.6 * 0.03 * 20 * 0.1 / _thickness;
	}
	void Set_Weak1(real _rho, real _gamma_water, real _thickness, int frame_rate) {
		Calculate_Alpha(_rho, _gamma_water, _thickness, frame_rate);
		laplacian_pressure_coeff = 2 * 0.0000005 * alpha_height_laplacian;
		divergence_pressure_coeff = 1 * 0.0005 * alpha_divergence;
		height_pressure_coeff = 0. * alpha_height;
		tangential_pressure_coeff = 1.;
		boundary_pressure_coeff = 1.0 / _thickness;
	}
	void Set_Weak2(real _rho, real _gamma_water, real _thickness, int frame_rate) {
		Calculate_Alpha(_rho, _gamma_water, _thickness, frame_rate);
		laplacian_pressure_coeff = 1 * 0.0000005 * alpha_height_laplacian;
		divergence_pressure_coeff = 1 * 0.0005 * alpha_divergence;
		height_pressure_coeff = 0. * alpha_height;
		tangential_pressure_coeff = 0.3;
		boundary_pressure_coeff = .2 / _thickness;
	}
	void Set_Zero(void) {
		laplacian_pressure_coeff = 0;
		divergence_pressure_coeff = 0;
		height_pressure_coeff = 0;
		tangential_pressure_coeff = 0;
		boundary_pressure_coeff = 0;
	}
};

class NormalPressureParams {
public:
	real capillary_coeff = 100000;
	real init_capillary_coeff = 0.;
	std::string air_pressure_mode = "none"; //{"none","ib","pv"}
	real atmo_pressure = 0;//0 for vaccum outside
	//for ib force
	real air_density = 1e-3;//rho_air/rho_water
	real ib_force_coeff = 1.0;
	//for pv force
	real pv_air_force_coeff = 1;
	bool closed = false;
	real pressure_constant = 1;
	real pressure_decay = 0.;
	real decayed_pressure = 1;
	OperatorParams opt_mode;
	NormalPressureParams() {
		opt_mode = OperatorParams("only_fluid", KernelType::SPIKY);
	}
	void Set_Sphere_Baseline(int d, real R, real gamma, real V, real rho, real thickness, real dp_multiplier=1) {
		real kappa = (d == 2) ? 1.0 / R : 2.0 / R; 
		real dp = 4 * gamma * kappa;
		real alpha_capillary = rho * thickness / dp * dp_multiplier;//dp to generate acceleration of dp_multiplier
		capillary_coeff = alpha_capillary;
		air_pressure_mode = "pv";
		atmo_pressure = dp * 10;
		pv_air_force_coeff = 1;
		closed = true;
		pressure_constant = V * (atmo_pressure + dp);
		opt_mode = OperatorParams("only_fluid", KernelType::SPIKY);
	}
	void Set_Sphere_Weak(int d, real R, real gamma, real V, real rho, real thickness) {
		real kappa = (d == 2) ? 1.0 / R : 2.0 / R;
		real dp = 4 * gamma * kappa;
		real alpha_capillary = rho * thickness / dp;//dp to generate acceleration of 1
		capillary_coeff = alpha_capillary;
		air_pressure_mode = "pv";
		atmo_pressure = dp * 10;
		pv_air_force_coeff = 1;
		closed = true;
		pressure_constant = V * (atmo_pressure + dp);
		opt_mode = OperatorParams("only_fluid", KernelType::SPIKY);
	}
	void Set_Hemisphere(int d, real R, real gamma, real rho, real thickness) {
		real kappa = (d == 2) ? 1.0 / R : 2.0 / R;
		real dp = 4 * gamma * kappa;
		real alpha_capillary = rho * thickness / dp;//dp to generate acceleration of 1
		capillary_coeff = alpha_capillary;
		air_pressure_mode = "none";
		closed = false;
		opt_mode = OperatorParams("only_fluid", KernelType::SPIKY);
	}
	void Set_Irregular(int d, real gamma, real V, real rho, real thickness) {
		real nominal_R = pow(V * 3.0 / (4.0 * pi), 1.0 / 3);
		real kappa = 2.0 / nominal_R;
		real dp = 4 * gamma * kappa;
		real alpha_capillary = rho * thickness / dp;//dp to generate acceleration of 1
		capillary_coeff = alpha_capillary * 5e-3;
		air_pressure_mode = "pv";
		atmo_pressure = dp * 5;
		pv_air_force_coeff = 1;
		closed = true;
		pressure_constant = V * (atmo_pressure + dp);
		opt_mode = OperatorParams("only_fluid", KernelType::SPIKY);
	}
	void Set_Circle_Baseline(real R, real gamma, real rho, real thickness) {
		real kappa = 2.0 / R;
		real dp = 4 * gamma * kappa;
		real alpha_capillary = rho * thickness / dp;//dp to generate acceleration of 1
		capillary_coeff = 5 * alpha_capillary;
		air_pressure_mode = "none";
		opt_mode = OperatorParams("all", KernelType::SPIKY);
	}
	void Set_IB(real R, real gamma, real rho, real thickness) {
		real kappa = 2.0 / R;
		real dp = 4 * gamma * kappa;
		real alpha_capillary = rho * thickness / dp;//dp to generate acceleration of 1
		capillary_coeff = alpha_capillary;
		std::cout << "capillary coeff: " << capillary_coeff << "\n";
		air_pressure_mode = "ib";
		closed = false;
		air_density = 1;
		ib_force_coeff = 1.0;
	}
	void Set_Planar(void) {
		air_pressure_mode = "none";
		capillary_coeff = 0;
	}
};

class BoundaryParams {
public:
	std::string particle_geometry_mode = "mirror_compensate";//how to calculate height and surface area for weighting, "traditional" or "surface" or "volumetric"
	bool replenish = false;
	std::string boundary_force_mode = "binary";//boundary force. "adhesion" or "binary" or "none"
	real vis_boundary = 150;
	real grav_boundary = 0.02 * 10000;
	real boundary_force_coeff = 1;
	std::string boundary_mode = "slippery"; //to post-process velocity after calculating F and V
	bool keep_xz_plane = false;
	real replenish_interval = 0.02;
	real replenish_proportion = 0.8;
	real replenish_dx_num = 3.0;
	real keep_velocity_rate = 0.0;
	std::string replenish_criteria = "none";
	bool inherit_rh = false;
	BoundaryParams(){}
	void Set_Circle_Baseline(real vis_coeff) {
		particle_geometry_mode = "mirror_compensate";
		replenish = false;
		boundary_force_mode = "vis";
		vis_boundary = vis_coeff;
		boundary_force_coeff = 1;
		boundary_mode = "analytical";
		keep_xz_plane = true;
	}
	void Set_Blow_Bubble(real vis_coeff) {
		particle_geometry_mode = "mirror_compensate";
		replenish = false;
		boundary_force_mode = "vis";
		vis_boundary = vis_coeff;
		boundary_force_coeff = 1;
		boundary_mode = "analytical";
		keep_xz_plane = false;
	}
	void Set_Pure_Gravity(real grav) {
		particle_geometry_mode = "mirror_compensate";
		replenish = false;
		boundary_force_mode = "binary";
		vis_boundary = 0;
		grav_boundary = grav;
		boundary_force_coeff = 1;
		boundary_mode = "analytical";
	}
	void Set_Pure_Analytical(void) {
		particle_geometry_mode = "mirror_compensate";
		replenish = false;
		boundary_force_mode = "none";
		boundary_mode = "analytical";
	}
	void Set_No_Boundary(void) {
		particle_geometry_mode = "mirror_compensate";
		replenish = false;
		boundary_force_mode = "none";
		boundary_mode = "none";
	}
};

class RenderHeightParams {
public:
	std::string mode = "divergence";//{"divergence","laplacian"}
	bool fix_RH = true;
	bool blend_h = false;
	real blend_constant = 0.;
	RenderHeightParams(){}
	RenderHeightParams(std::string _mode, bool _fix_rh, bool _blend, real _blend_constant) {
		mode = _mode;
		fix_RH = _fix_rh;
		blend_h = _blend;
		blend_constant = _blend_constant;
	}
};

#endif

