#ifndef __LeastSquaresDriver_h__
#define __LeastSquaresDriver_h__
#include "Driver.h"
#include "LeastSquares.h"
#include "Mesh.h"

template<int d> class LeastSquaresDriver : public Driver
{public:
	using Base=Driver;
	using VectorDi=Vector<int,d>;using VectorD=Vector<real,d>;
	using VectorTi=Vector<int,d-1>;using VectorT=Vector<real,d-1>;

	GeometryParticles<d> particles;
	Array<real> data,approx;
	SurfaceMesh<d> mesh;
	Vector2 domain=Vector2(-1.,1.);
	int count=16;

	virtual void Initialize()
	{
		using namespace LeastSquares;
		switch(test){
		case 1:{////LS<1,1>
			if constexpr (d==2){
			LS<1,1> ls;
			int pn=3;
			particles.Resize(pn);
			particles.X(0)=VectorD::Zero();
			particles.X(1)=VectorD::Ones();
			particles.X(2)=VectorD::Ones()*(real)2;
			Update_Data(pn,1);
			ls.Fit(&data[0],pn);
			for (int i=0; i<pn; i++)
			{approx.push_back(ls(data[i*2]));}
			real sqerr=Calculate_LS_Error(&data[0],&approx[0],1);
			std::cout<<"Test LS<1,1>: "<<std::endl;
			std::cout<<"Expect c=["<<0<<","<<1<<"] "<<std::endl;
			std::cout<<"Approx c=["<<ls.c[0]<<","<<ls.c[1]<<"] "<<std::endl;
			std::cout<<"square err="<<sqerr<<std::endl;}
		}break;
		case 2: {////LS<1,2>
			if constexpr (d==2){
			LS<1,2> ls;////f(x)=c0+c1*x+c2*x^2			
			int pn=3;
			particles.Resize(pn);
			particles.X(0)=VectorD(0.,-3.);
			particles.X(1)=VectorD(1.,5.);
			particles.X(2)=VectorD(-1.,1.);
			Update_Data(pn,1);
			ls.Fit(&data[0],pn);
			for (int i=0; i<pn; i++)
			{approx.push_back(ls(data[i*2]));}
			real sqerr=Calculate_LS_Error(&data[0],&approx[0],1);
			std::cout<<"Test LS<1,2>:" <<std::endl;
			std::cout<<"Expect c=["<<-3<<","<<2<<","<<6<<"] "<<std::endl;
			std::cout<<"Approx c=["<<ls.c[0]<<","<<ls.c[1]<<","<<ls.c[2]<< "] "<<std::endl;
			std::cout<<"square err="<<sqerr<<std::endl;
			}
		}break;
		case 3:{////LS<1,3>
			if constexpr (d==2) {////f(x)=x^3+x^2+x+2
			LS<1,3> ls;////f(x)=c0+c1*x+c2*x^2+c3*x^3
			int pn=4;
			particles.Resize(4);
			particles.X(0)=VectorD(0.,2.);
			particles.X(1)=VectorD(1.,5.);
			particles.X(2)=VectorD(-1.,1.);
			particles.X(3)=VectorD(2.,16.);
			Update_Data(pn,1);
			ls.Fit(&data[0],pn);
			for (int i=0; i<pn; i++)
			{approx.push_back(ls(data[i*2]));}
			real sqerr=Calculate_LS_Error(&data[0],&approx[0],1);
			std::cout<<"Test LS<1,3>: "<<std::endl;		
			std::cout<<"Expect c=["<<2<<","<<1<<","<<1<<","<<1<<"] "<<std::endl;
			std::cout<<"Approx c=["<<ls.c[0]<<","<<ls.c[1]<<","<<ls.c[2]<<","<<ls.c[3]<< "] "<<std::endl;
			std::cout<<"square err="<<sqerr<<std::endl;}
		}break;
		case 4: {////LS<2,1>
		}break;
		case 5:{////LS<2,2>
			if constexpr (d==3) {////f(x,y)=c0+c1*x+c2*y+c3*x^2+c4*y^2+c5*xy
				LS<2,2> ls;
				int pn=9;
				particles.Resize(pn);////data from http://www.nealen.com/projects/mls/asapmls.pdf
				particles.X(0)=VectorD(1.,1.,1.);particles.X(1)=VectorD(1.,-1.,-0.5);particles.X(2)=VectorD(-1.,1.,1.);
				particles.X(3)=VectorD(-1,-1,1.);particles.X(4)=VectorD(0.,0.,-1.);particles.X(5)=VectorD(1.,0.,0.);
				particles.X(6)=VectorD(-1.,0.,0.);particles.X(7)=VectorD(0.,1.,0.);particles.X(8)=VectorD(0.,-1.,0);
				Update_Data(pn,2);
				ls.Fit(&data[0],pn);
				for (int i=0; i<pn; i++)
				{ approx.push_back(ls(data[i*3],data[i*3+1]));}
				real sqerr=Calculate_LS_Error(&data[0],&approx[0],2);
				std::cout<<"Test LS<2,2>: "<<std::endl;		
				std::cout<<"Approx c=["<<ls.c[0]<<","<<ls.c[1]<<","<<ls.c[2]<<","<<ls.c[3]<<","<<ls.c[4]<<","<<ls.c[5]<<"]"<<std::endl;
				std::cout<<"square err="<<sqerr<<std::endl;
			}
		}break;
		case 6: {////MLS<2,2>
			if constexpr (d==3) {
				MLS<2,2> mls;
				int pn=3;
				particles.X(0)=VectorD::Zero();
				particles.X(1)=VectorD::Ones();
				particles.X(2)=VectorD::Ones()*(real)2;
				Update_Data(pn,2);
				mls.Fit(&data[0],pn,4.,4.);
				for (int i=0; i<pn; i++)
				{ approx.push_back(mls(data[i*3],data[i*3+1]));}
				//real sqerr=Calculate_MLS_Error(&data[0],&approx[0],2);
				std::cout<<"Test MLS<2,2>: "<<std::endl;		
				std::cout<<"Approx c=["<<mls.c[0]<<","<<mls.c[1]<<","<<mls.c[2]<<","<<mls.c[3]<<","<<mls.c[4]<<","<<mls.c[5]<<"]"<<std::endl;
				//std::cout<<"square err="<<sqerr<<std::endl;
			}
		}break;
		case 7:{////WLS<2,1>
		}break;
		}

		/*
		real dx=(domain[1]-domain[0])/(real)count;
		if constexpr (d==2){
			mesh.Vertices().push_back(Vector2(domain[0],(real)0));
			for(int i=0;i<count;i++){
				mesh.Vertices().push_back(Vector2(domain[0]+dx*(real)i,(real)0));
				mesh.Elements().push_back(Vector2i(i,i+1));}}
		else if constexpr (d==3){
			MeshFunc::Initialize_Herring_Bone_Mesh(count,count,dx,&mesh,0,1);}
			*/
	}

	virtual void Run()
	{
		Update_Least_Squares();
		Write_Output_Files(0);	
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);

		std::cout<<"particles ";
		for(int i=0;i<particles.Size();i++){
			std::cout<<particles.X(i).transpose()<<", ";}

		{std::string file_name=frame_dir+"/particles";
		particles.Write_To_File_3d(file_name);}

		{std::string file_name=frame_dir+(d==2?"/segment_mesh":"/triangle_mesh");
		mesh.Write_To_File_3d(file_name);}

		std::cout<<"Write to frame "<<frame<<std::endl;
	}

	void Test_Least_Squares()
	{
		using namespace LeastSquares;
		/*
		{Array<real> data={0.,1., 1.,5., -1.,9., 5.,-339};////f(x)=-4x^3+6x^2+2x+1
		LS<1,3> ls;////f(x)=c0+c1*x+c2*x^2+c3*x^3
		ls.Fit(&data[0],4);
		real val=ls(0.5);
		real val_2=ls(2.);	
		std::cout<<"Test LS<1,3>: "<<std::endl;
		std::cout<<"Expect c=["<<2<<","<<1<<","<<1<<","<<1<<"] "<<std::endl;
		std::cout<<"Approx c=["<<ls.c[0]<<","<<ls.c[1]<<","<<ls.c[2]<<","<<ls.c[3]<< "] "<<std::endl;
		std::cout<<"val(expect 3): "<<val<<std::endl<<std::endl;
		}*/

		{Array<real> data={0.,0.,0., 1.,1.,1., 2.,2.,2.};
		MLS<2,2> mls;
		mls.Fit(&data[0],3,4.,4.);
		real val=mls(4.,4.);
		std::cout<<"Test MLS<2,2> ";
		std::cout<<"val: "<<val<<std::endl<<std::endl;}
	}

	void Update_Least_Squares()
	{
		/*
		using namespace LeastSquares;
		LS<d-1,3> ls;
		Array<real> data;
		for(int i=0;i<particles.Size();i++){
			for(int j=0;j<d;j++){
				data.push_back(particles.X(i)[j]);}}
		ls.Fit(&data[0],particles.Size());
		*/
		//for(int i=0;i<mesh.Vertices().size();i++){
			//mesh.Vertices()[i][d-1]=ls(V<d-1>(mesh.Vertices()[i]));
			//std::cout<<"vertex: "<<mesh.Vertices()[i].transpose()<<std::endl;}
	}

	void Update_Data(int n, int deg)
	{
		data.clear();
		for(int i=0;i<n;i++){
			for(int j=0;j<=deg;j++){
				data.push_back(particles.X(i)[j]);}}
	}

	real Calculate_LS_Error(const real* data, const real* approx, int dim)
	{
		real sum_sqerr=(real)0;
		for (int i = 0; i < particles.Size(); i++)
		{
			sum_sqerr+=pow((data[(dim+1)*(i+dim)]-approx[i]),2);
		}
		return sum_sqerr;
	}

	real Calculate_WLS_Error(const real* data, const real* approx, int dim, real x, real y)
	{
		real sum_sqerr=(real)0;
		for (int i = 0; i < particles.Size(); i++)
		{
			sum_sqerr+= pow((data[(dim+1)*(i+dim)]-approx[i]),2);
		}
		return sum_sqerr;
	}
};

#endif