#ifndef __SPHBubbleParticles_h__
#define __SPHBubbleParticles_h__
#include "GeometryParticles.h"



template<int d> class SPHBubbleParticles: public GeometryParticles<d>
{Typedef_VectorDii(d);Typedef_MatrixD(d);using Base=GeometryParticles<d>;
public:
	using Base::Size; using Base::Resize; using Base::Add_Element; using Base::Add_Elements; using Base::Join; using Base::Copy_Element_From; using Base::Delete_Elements; using Base::Print_Attributes; using Base::Save_Snapshot; using Base::Load_Snapshot;
	using Base::X;using Base::V;using Base::F;using Base::M;
	using Base::I;using Base::E;using Base::G;
	using MatrixT=Matrix<real,d-1>;					////tangential vector type
	using VectorT=Vector<real,d-1>;					////tangential vector type
	using VectorTi=Vector<int,d-1>;					////tangential vector int

	SPHBubbleParticles() {
		New_Attributes();
		std::cout << "SPHBubbleParticles initialized\n";
		Points<d>::Print_Attributes();
	}

	//////Particle attributes

	////Particle attributes
	Declare_Attribute(real, H, h, Rebind_H);        ////height of current cylinder
	Declare_Attribute(real, KH, kh, Rebind_KH);     ////laplacian of H. "Additional" mean curvature
	Declare_Attribute(real, RH, rh, Rebind_RH);     ////height for rendering
	Declare_Attribute(real, RH_V, rh_v, Rebind_RH_V);////the change rate of rh
	Declare_Attribute(real, Vrt, vrt, Rebind_Vor);//vorticity strength

	////tangential space scalar attributes
	Declare_Attribute(real,Vol,vol,Rebind_Vol);				////volume
	Declare_Attribute(real, Div, div, Rebind_Div);			////divergence
	Declare_Attribute(real,P,p,Rebind_P);					////pressure
	Declare_Attribute(real, P0, p0, Rebind_P0);					////0-degree pressure
	Declare_Attribute(real, P1, p1, Rebind_P1);					////1-degree pressure
	Declare_Attribute(real, P2, p2, Rebind_P2);					////2-degree pressure
	Declare_Attribute(real, PB, pb, Rebind_PB);
	
	//Surfactant Tension Gamma
	Declare_Attribute(real, Gamma, gamma, Rebind_Gamma);
	//Soap Concentration Conc
	Declare_Attribute(real, Conc, conc, Rebind_Conc);

	//Surface area
	Declare_Attribute(real, SA, sa, Rebind_SA);

	//Phase
	Declare_Attribute(int, Phase, phase, Rebind_Phase);


	//Surface area
	Declare_Attribute(VectorD, SN, sn, Rebind_SN);

	////set boundary particles
	Declare_Attribute(int,B,b,Rebind_B); ////1 for boundary particles
	bool Is_Boundary(int idx)const{return (B(idx)==1);}

	Declare_Attribute_Inherent_Func(h,kh,rh,rh_v,vrt,vol,div,p,p0,p1,p2, pb,b,gamma,conc, sa, sn, phase);	////does not need the attributes from the base class


	void Apply_Rotation(int idx, MatrixD R) {
		X(idx) = R * X(idx);
		V(idx) = R * V(idx);
		F(idx) = R * F(idx);
		E(idx) = R * E(idx);
	}
	void Apply_Translation(int idx, VectorD T) {
		X(idx) = X(idx) + T;
	}
};

#endif