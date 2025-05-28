//////////////////////////////////////////////////////////////////////////
// SPH Kernels
// Copyright (c) (2018-),Xiangxin Kong, Mengdi Wang
// Please see simplex/docs/kernels-math-en.md for documentation.
// This file is part of SimpleX,whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __Kernel_h__
#define __Kernel_h__
#include <iostream>
#include <cmath>
#include "Common.h"
#include "Constants.h"

enum class KernelType { POLY6, SPIKY, CUBIC, QUINTIC, GAUSSIAN };
//Note: in mathematical formulations, there're some kernels truncated to 2h or 3h instead of h.
//However we want them to truncate to h, so some are not corresponding well to formulations.
//We can confirm that POLY6, SPIKY, CUBIC, QUINTIC, GAUSSIAN are normalized and numerically validated.
//But we're not quite sure about another kernels' behaviors.
//TODO: fix this.

template<int d> class KernelSPH
{
	Typedef_VectorDii(d);
public:

	real gaussian_trunc = 3;//h=3*sigma, for gaussian kernel

	//Truncated at h. Only non-zero in [0,h)
	real h, h2, h3;//h,h^2,h^3
	real hd1;//1/h
	real h13, h23;//h^(1/3), h^(2/3)

	real coef_poly6, coef_poly6_grad, coef_poly6_lap;
	real coef_spiky, coef_spiky_grad, coef_spiky_lap;
	real coef_vis, coef_vis_grad_1, coef_vis_grad_2, coef_vis_grad_3, coef_vis_lap;
	real coef_cubic_spline;
	real coef_quintic;
	real coef_gaussian, coef_gaussian_grad;

	KernelSPH() { Precompute_Coefs(1.0); }
	KernelSPH(const real _h) { Precompute_Coefs(_h); }
	void Initialize(const real _h) { Precompute_Coefs(_h); }
	void Precompute_Coefs(const real _h);

	real Weight(const real x, const KernelType& kernel_type)const;
	VectorD Grad(const VectorD& vec, const KernelType& kernel_type)const;


	//////////////////////////////////////////////////////////////////////////
	////kernels

	real W_Poly6(const real r) const;
	real W_Spiky(const real r) const;
	real W_Cubic(const real r) const;
	real W_Quintic(const real r) const;
	real W_Gaussian(const real r) const;

	//////////////////////////////////////////////////////////////////////////
	////gradients

	VectorD Grad_Poly6(const VectorD& vr) const;
	VectorD Grad_Spiky(const VectorD& vr) const;
	VectorD Grad_Vis(const VectorD& vr) const;
	VectorD Grad_Cubic(const VectorD& vr) const;
	VectorD Grad_Quintic(const VectorD& vr) const;
	VectorD Grad_Gaussian(const VectorD& vr) const;

	//////////////////////////////////////////////////////////////////////////
	////Laplacian
	real Lap_Poly6(const VectorD& vr) const;
	real Lap_Spiky(const VectorD& vr) const;
	real Lap_Vis(const VectorD& vr) const;
};
#endif
