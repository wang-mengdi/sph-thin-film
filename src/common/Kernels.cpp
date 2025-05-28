#include "Kernels.h"

template<int d>
void KernelSPH<d>::Precompute_Coefs(const real _h)
{
	h = _h;
	h2 = h * h;
	h3 = h * h * h;
	hd1 = (real)1 / h;
	h13 = h * one_third;
	h23 = h * two_thirds;

	if constexpr (d == 3) {
		////poly6
		coef_poly6 = 315.0 / (64.0 * pi * pow(h, 9));
		coef_poly6_grad = -945.0 / (32.0 * pi * pow(h, 9));
		coef_poly6_lap = coef_poly6_grad;
		////spiky
		coef_spiky = 15.0 / (pi * pow(h, 6));
		coef_spiky_grad = -45.0 / (pi * pow(h, 6));
		coef_spiky_lap = -90.0 / (pi * pow(h, 6));
		////vis
		coef_vis = 15.0 / (2.0 * pi * h3);
		coef_vis_grad_1 = coef_vis;
		coef_vis_grad_2 = -1.5 / h3;
		coef_vis_grad_3 = 2.0 / h2;
		coef_vis_lap = 45.0 / (pi * h3);
		////cubic
		// coef_cubic_spline=3.0/(2.0*pi*h3);
		coef_cubic_spline = 8.0 / (pi * h3);
		////interpolation4
		coef_quintic = 81.0 / (359.0 * pi * h3);
		////Gaussian kernel
		coef_gaussian = pow(gaussian_trunc, 3) / (pow(pi, 1.5) * h3);
		coef_gaussian_grad = coef_gaussian * (-2) * pow(gaussian_trunc, 2) / pow(h, 2);
	}
	else if constexpr (d == 2) {
		////poly6
		coef_poly6 = 4.0 / (pi * pow(h, 8));
		coef_poly6_grad = -24.0 / (pi * pow(h, 8));
		coef_poly6_lap = coef_poly6_grad;
		////spiky
		coef_spiky = 10.0 / (pi * pow(h, 5));
		coef_spiky_grad = -30.0 / (pi * pow(h, 5));
		coef_spiky_lap = -60.0 / (pi * pow(h, 5));
		////vis
		coef_vis = 10.0 / (3.0 * pi * h2);
		coef_vis_grad_1 = coef_vis;
		coef_vis_grad_2 = -1.5 / h3;
		coef_vis_grad_3 = 2.0 / h2;
		coef_vis_lap = 20.0 / (pi * h2);
		////cubic
		// coef_cubic_spline=15.0/(7.0*pi*h2);
		coef_cubic_spline = 40.0 / (7.0 * pi * h2);
		////interpolation4
		coef_quintic = 63.0 / (478.0 * pi * h2);
		////Gaussian kernel
		coef_gaussian = pow(gaussian_trunc, 2) / (pi * h2);
		coef_gaussian_grad = coef_gaussian * (-2) * pow(gaussian_trunc, 2) / pow(h, 2);
	}
	else if constexpr (d == 1) {
		////poly6
		coef_poly6 = 35.0 / (32.0 * pow(h, 7));
		coef_poly6_grad = -105.0 / (16.0 * pow(h, 7));
		coef_poly6_lap = coef_poly6_grad;
		////spiky
		coef_spiky = 2.0 / pow(h, 4);
		coef_spiky_grad = -6.0 / pow(h, 4);
		coef_spiky_lap = -12.0 / pow(h, 4);
		////vis
		coef_vis = (real)0;
		coef_vis_grad_1 = coef_vis;
		coef_vis_grad_2 = -1.5 / h3;
		coef_vis_grad_3 = 2.0 / h2;
		coef_vis_lap = (real)0;
		////cubic
		// coef_cubic_spline=1.0/h;
		coef_cubic_spline = 4.0 / (3.0 * h);
		////interpolation4
		coef_quintic = 1.0 / (40 * h);
		////Gaussian kernel
		coef_gaussian = gaussian_trunc / (sqrt(pi) * h);
		coef_gaussian_grad = coef_gaussian * (-2) * pow(gaussian_trunc, 2) / pow(h, 2);
		
	}
}

template<int d>
real KernelSPH<d>::Weight(const real x, const KernelType& kernel_type) const
{
	if (kernel_type == KernelType::POLY6) {
		return W_Poly6(x);
	}
	else if (kernel_type == KernelType::SPIKY) {
		return W_Spiky(x);
	}
	else if (kernel_type == KernelType::CUBIC) {
		return W_Cubic(x);
	}
	else if (kernel_type == KernelType::QUINTIC) {
		return W_Quintic(x);
	}
	else if (kernel_type == KernelType::GAUSSIAN) {
		return W_Gaussian(x);
	}
	else {
		std::cerr << "KernelSPH<d>::Weight error: unknown type\n";
		assert(false);
		return 0.0;
	}
}

template<int d>
Vector<real,d> KernelSPH<d>::Grad(const VectorD& vec, const KernelType& kernel_type) const
{
	if (kernel_type == KernelType::POLY6) {
		return Grad_Poly6(vec);
	}
	else if(kernel_type == KernelType::SPIKY) {
		return Grad_Spiky(vec);
	}
	else if (kernel_type == KernelType::CUBIC) {
		return Grad_Cubic(vec);
	}
	else if (kernel_type == KernelType::QUINTIC) {
		return Grad_Quintic(vec);
	}
	else if (kernel_type == KernelType::GAUSSIAN) {
		return Grad_Gaussian(vec);
	}
	else {
		std::cerr << "Surface_Kernel_Grad error: unknown type\n";
		assert(false);
		return VectorD::Zero();
	}
}

template<int d>
real KernelSPH<d>::W_Poly6(const real r) const
{
	if (r < h) return coef_poly6 * pow(h2 - r * r, 3);
	else return (real)0;
}

template<int d>
real KernelSPH<d>::W_Spiky(const real r) const
{
	if (r < h) return coef_spiky * pow(h - r, 3);
	else return (real)0;
}

template<int d>
real KernelSPH<d>::W_Cubic(const real r) const
{
	//Truncate to [0,h). So h here is 2h in mathematical formula
	real u = r / h;
	if (u >= 0 && u <= 0.5) return coef_cubic_spline * ((real)6 * pow(u, 3) - (real)6 * pow(u, 2) + 1);
	else if (u > 0.5 && u < 1) return coef_cubic_spline * (2 * pow(1 - u, 3));
	else return (real)0;
}

template<int d>
real KernelSPH<d>::W_Quintic(const real r) const
{
	if (r >= 0 && r < h13) {
		return coef_quintic * (pow(3 - 3 * r * hd1, 5) - 6 * pow(2 - 3 * r * hd1, 5) + 15 * pow(1 - 3 * r * hd1, 5));
	}
	else if (r >= h13 && r < h23) {
		return coef_quintic * (pow(3 - 3 * r * hd1, 5) - 6 * pow(2 - 3 * r * hd1, 5));
	}
	else if (r >= h23 && r < h) {
		return coef_quintic * (pow(3 - 3 * r * hd1, 5));
	}
	else { return 0; }
}

template<int d>
real KernelSPH<d>::W_Gaussian(const real r) const
{
	return coef_gaussian * exp(-pow(gaussian_trunc * r / h, 2));
}


template<int d>
Vector<real, d> KernelSPH<d>::Grad_Poly6(const VectorD& vr) const
{
	real r = vr.norm();
	if (r > 0 && r < h)
	{
		return coef_poly6_grad * pow(h2 - r * r, 2) * vr;
	}
	else { return VectorD::Zero(); }
}

template<int d>
Vector<real, d> KernelSPH<d>::Grad_Spiky(const VectorD& vr) const
{
	real r = vr.norm();
	if (r > 0 && r < h) return coef_spiky_grad * pow(h - r, 2) * vr / r;
	else return VectorD::Zero();
}

template<int d>
Vector<real, d> KernelSPH<d>::Grad_Vis(const VectorD& vr) const
{
	real r = vr.norm();
	if (r > 0 && r < h) return coef_vis_grad_1 * (coef_vis_grad_2 * r + coef_vis_grad_3 - h / (2 * r * r * r)) * vr;
	else return VectorD::Zero();
}

template<int d>
Vector<real, d> KernelSPH<d>::Grad_Cubic(const VectorD& vr) const
{
	real r = vr.norm();
	real u = r / h;
	const VectorD gradu = vr / (r * h);
	if (u > 0 && u <= 0.5) return 6.0 * coef_cubic_spline * u * (3.0 * u - 2.0) * gradu;
	else if (u > 0.5 && u < 1) return -6.0 * coef_cubic_spline * pow(1.0 - u, 2) * gradu;
	else return VectorD::Zero();
}

template<int d>
Vector<real, d> KernelSPH<d>::Grad_Quintic(const VectorD& vr) const
{
	real r = vr.norm();
	if (r > 0 && r < h13) {
		return vr / r * coef_quintic * hd1 * (-(real)15 * pow(3 - 3 * r * hd1, 4) + (real)90 * pow(2 - 3 * r * hd1, 4) - (real)225 * pow(1 - 3 * r * hd1, 4));
	}
	else if (r >= h13 && r < h23) {
		return vr / r * coef_quintic * hd1 * (-(real)15 * pow(3 - 3 * r * hd1, 4) + (real)90 * pow(2 - 3 * r * hd1, 4));
	}
	else if (r >= h23 && r < h) {
		return vr / r * coef_quintic * hd1 * (-(real)15 * pow(3 - 3 * r * hd1, 4));
	}
	else { return VectorD::Zero(); }
}

template<int d>
Vector<real, d> KernelSPH<d>::Grad_Gaussian(const VectorD& vr) const
{
	real r = vr.norm();
	return coef_gaussian_grad * exp(-pow(gaussian_trunc * r / h, 2)) * vr;
}


template<int d>
real KernelSPH<d>::Lap_Poly6(const VectorD& vr) const
{
	real r = vr.norm();
	if (r <= h) return coef_poly6_lap * (h2 - r * r) * (3 * h2 - 7 * r * r);
	else return (real)0;
}

template<int d>
real KernelSPH<d>::Lap_Spiky(const VectorD& vr) const
{
	real r = vr.norm();
	if (r > 0 && r <= h) return coef_spiky_lap * (h - r) * (h - 2.0 * r) / r;
	else return (real)0;
}

template<int d>
real KernelSPH<d>::Lap_Vis(const VectorD& vr) const
{
	real r = vr.norm();
	if (r < h) return coef_vis_lap * (h - r) / h3;
	else return (real)0;
}

template class KernelSPH<1>;
template class KernelSPH<2>;
template class KernelSPH<3>;
