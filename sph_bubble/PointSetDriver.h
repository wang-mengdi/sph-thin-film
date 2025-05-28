#ifndef __PointSetDriver_h__
#define __PointSetDriver_h__

#include "PointSetDriver.h"
#include "Driver.h"
#include "PointSet.h"
#include "PointSetFunc.h"
#include "Mesh.h"
#include "AnalyticalFields.h"
#include "LeastSquares.h"
#include "RandomNumber.h"

using namespace AuxFunc;

template<int d> class PointSetDriver : public Driver
{Typedef_VectorDii(d);Typedef_MatrixD(d);using Base=Driver;
public:
	PointSet<d> ps;
	VorticityField<d> field;
	GeometryParticles<d> aux_points;

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		switch(test){
		case 1:
		case 2:
		case 4:{
			for(int i=0;i<ps.points->Size();i++){
				ps.points->V(i)=field.Velocity(ps.points->X(i),time);
				ps.points->X(i)+=ps.points->V(i)*dt;}
			ps.Update();
			//ps.Update_Local_Frame(dt);
			ps.Reinitialize_Local_Frames();	
			ps.Point_Reseeding();
			//ps.Point_Relaxation();
		}break;
		case 3:{
			for(int i=0;i<ps.points->Size();i++){
				ps.points->V(i)=field.Velocity(ps.points->X(i),time);
				ps.points->X(i)+=ps.points->V(i)*dt;}
			ps.Update();
			//ps.Update_Local_Frame(dt);
			ps.Reinitialize_Local_Frames();		

			for(int i=0;i<aux_points.Size();i++){
				aux_points.X(i)=ps.Project_To_Surface(aux_points.X(i));}
		}break;
		}
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);

		{std::string file_name=frame_dir+"/tracker_circles";
		PointSetFunc::Write_Tracker_Circles_To_File<d>(file_name,*ps.points);}

		{std::string file_name=frame_dir+"/segment_mesh";
		PointSetFunc::Write_Local_Frames_To_File(file_name,*ps.points);}
		
		if(aux_points.Size()!=0)
		{std::string file_name=frame_dir+"/tracker_points";
		aux_points.Write_To_File_3d_Fast(file_name);}

		std::cout<<"Write to frame "<<frame<<std::endl;
	}

	virtual void Initialize()
	{
		frame_rate=25;

		switch(test){
		case 1:{	////2D circle and 3D sphere, rigid rotation
			if constexpr (d==2){
				real dx=PointSetFunc::Initialize_Circle_Points(VectorD::Zero(),(real)1,64,*ps.points);
				ps.Initialize(dx,2);}
			else if constexpr (d==3){
				real dx=PointSetFunc::Initialize_Sphere_Points(VectorD::Zero(),(real)1,3,*ps.points);
				ps.Initialize(dx,2);}
			field.test=1;
		}break;
		case 2:{	////vortex motion
			VectorD c=VectorD::Unit(1)*(real).5;real R=(real).2;
			if constexpr (d==2){
				real dx=PointSetFunc::Initialize_Circle_Points(c,R,64,*ps.points);
				ps.Initialize(dx,2);}
			else if constexpr (d==3){
				real dx=PointSetFunc::Initialize_Sphere_Points(VectorD::Zero(),(real)1,3,*ps.points);
				ps.Initialize(dx,2);}
			field.test=2;
		}break;
		case 3:{	////test projection
			real dx;
			if constexpr (d==2){
				dx=PointSetFunc::Initialize_Circle_Points(VectorD::Zero(),(real)1,64,*ps.points);
				ps.Initialize(dx,2);}
			else if constexpr (d==3){
				dx=PointSetFunc::Initialize_Sphere_Points(VectorD::Zero(),(real)1,3,*ps.points);
				ps.Initialize(dx,2);}

			RandomNumber random((real)-1,(real)1);
			for(int i=0;i<ps.points->Size();i++){
				for(int k=0;k<4;k++){
					VectorD perturb=(real)2*dx*random.VectorValue<d>()+ps.points->X(i);
					int s=aux_points.Add_Element();
					aux_points.X(s)=perturb;}}
			field.test=2;
		}break;
		case 4:{	////test relaxation
			real dx;
			if constexpr (d==2){
				dx=PointSetFunc::Initialize_Circle_Points(VectorD::Zero(),(real)1,64,*ps.points);
				ps.Initialize(dx,2);}
			else if constexpr (d==3){
				dx=PointSetFunc::Initialize_Sphere_Points(VectorD::Zero(),(real)1,3,*ps.points);
				ps.Initialize(dx,2);}	
			ps.Update();

			////perturb the initial positions
			RandomNumber random((real)-1,(real)1);
			Array<VectorD> perturbed_x(ps.points->Size());
			for(int i=0;i<ps.points->Size();i++){
				Vector<real,d-1> t=(real)1.*dx*random.VectorValue<d-1>();
				VectorD p;ps.Unproject_To_World(t,ps.points->E(i),p);
				perturbed_x[i]=ps.points->X(i)+p;
				perturbed_x[i]=ps.Project_To_Surface(perturbed_x[i]);}
			for(int i=0;i<ps.points->Size();i++){
				ps.points->X(i)=perturbed_x[i];}
		}break;
		case 5:{	////test string in 2D
			if constexpr (d==2){
				
			}break;
		}break;
		}
	}
};

#endif