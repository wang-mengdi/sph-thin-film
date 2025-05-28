#ifndef __AnimationDriver_h__
#define __AnimationDriver_h__
#include "Driver.h"
#include "AnimationSurface.h"
#include "Mesh.h"
#ifdef USE_TINY_OBJ_LOADER
#include "TinyObjLoader.h"
#include "tiny_obj_loader.h"
#endif

using std::ofstream;
using std::string;
using std::endl;

template<int d> class AnimationDriver : public Driver
{Typedef_VectorDii(d);
public:
	using Base=Driver;
	using Base::current_frame;using Base::last_frame;

	AnimationSurface<d> surface;
	
	Array<std::string> mesh_names={"wolfman"};
	std::string mesh_prefix;

	virtual void Run()
	{
		while(current_frame<last_frame){			
			Load_Obj(current_frame);
			std::cout<<"current_frame="<<current_frame<<std::endl;
			surface.generate_voronoi_mesh();
			surface.generate_noisy_points();
			Write_Output_Files(current_frame);
			current_frame++; // start with frame #0
		}
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);

		////Write Particles

		/*
		{std::string file_name = frame_dir + "/particles";
		// temp: write voronoi centers
		Particles<d> temp_particles;
		temp_particles.Resize(surface.voronoi_centers.size());
		temp_particles.XRef() = surface.voronoi_centers;
		temp_particles.Write_To_File_3d(file_name);}
		//*/
		/*
		{std::string file_name = frame_dir + "/particles";
		// temp: write noisy points
		Particles<d> temp_particles;
		temp_particles.Resize(surface.noisy_points.size());
		temp_particles.XRef() = surface.noisy_points;
		for (int i = 0; i < temp_particles.Size(); i++) {
			temp_particles.M(i) = (real)1;
			temp_particles.F(i) = VectorD::Zero();
			temp_particles.V(i) = VectorD::Zero();
		}
		temp_particles.Write_To_File_3d(file_name);}
		//*/


		//*
		{std::string file_name = frame_dir + "/particles";
		surface.points->Write_To_File_3d(file_name);}

		{std::string file_name=frame_dir+"/tracker_points";
		surface.points->Write_To_File_3d_Fast(file_name);}

		{std::string file_name=frame_dir+"/point_force";
		Write_Segments_To_File_3d_Fast<d,real>(surface.points->XRef(),surface.points->FRef(),file_name);}

		{std::string file_name=frame_dir+"/point_velocity";
		Write_Segments_To_File_3d_Fast<d,real>(surface.points->XRef(),surface.points->VRef(),file_name);}

		{std::string file_name=frame_dir+"/tracker_circles";
		int pn=surface.points->Size();Array<VectorD> normals(pn);
		for(int i=0;i<pn;i++)normals[i]=surface.points->Normal(i);
		Write_Vectors_To_File_3d_Fast<d,real>(surface.points->XRef(),normals,file_name);}
		//*/

		// write the traingle mesh
		{std::string file_name = frame_dir + "/triangle_mesh";
		surface.triangle_mesh->Write_To_File_3d(file_name); }

		std::cout<<"Write to frame "<<frame<<std::endl;

		Write_Obj_Files(frame_dir);
	}
	
	void Write_Obj_Files(const std::string & frame_dir)
	{
		// write noisy points
		string noise_path = frame_dir + "/" + "noisy.obj";
		ofstream noise_file;
		noise_file.open(noise_path);
		write_points(noise_file, surface.noisy_points);
		noise_file.close();

		//write original meshes (points and triangle faces)
		string original_path = frame_dir + "/" + "original.obj";
		ofstream original_file;
		original_file.open(original_path);
		write_points(original_file, surface.triangle_mesh->Vertices());
		write_faces(original_file, surface.triangle_mesh->Elements());
		original_file.close();

		// write voronoi meshes (points and lines)
		string voronoi_path = frame_dir + "/" + "voronoi.obj";
		ofstream voronoi_file;
		voronoi_file.open(voronoi_path);
		write_points(voronoi_file, surface.voronoi_centers);
		write_lines(voronoi_file, surface.voronoi_edges);
		original_file.close();
	}

	void write_points(ofstream & out, Array<VectorD> & points) {
		for (int i = 0; i < points.size(); i++) {
			out << "v";
			for (int j = 0; j < points[i].size(); j++) {
				out << " " << points[i][j];
			}
			out << endl;
		}
	}

	void write_lines(ofstream& out, Array<Vector2i>& edges) {
		for (int i = 0; i < edges.size(); i++) {
			out << "l " << edges[i][0]+1 << " " << edges[i][1]+1 << endl;
		}
	}

	void write_faces(ofstream& out, Array<Vector3i>& elements) {
		for (int i = 0; i < elements.size(); i++) {
			out << "f " << elements[i][0] + 1 << " " << elements[i][1] + 1 << " " << elements[i][2] + 1 << endl;
		}
	}

	virtual void Initialize()
	{
		using namespace Obj;
		switch(test){
		case 1:{	
			if constexpr (d==3){
				mesh_prefix=Path::Data()+"/meshes/"+mesh_names[0]+"/"+mesh_names[0]+"_";
			}
		}break;
		}
	}

	virtual void Load_Obj(int frame)
	{
		if constexpr (d==3){
		auto obj_frame=frame+1;
		std::string surface_mesh_file_name=mesh_prefix+std::to_string(obj_frame) +".obj";
		Array<std::shared_ptr<TriangleMesh<3>>> obj_mesh;
		Obj::Read_From_Obj_File<TriangleMesh<3>>(surface_mesh_file_name, obj_mesh);
		bool read=obj_mesh.size()>0;if(!read)return;
		
		surface.triangle_mesh=obj_mesh[0];

		MeshFunc::Rescale<3>(surface.triangle_mesh->Vertices(),(real)2);
		surface.points->Resize((int)surface.triangle_mesh->Vertices().size());
		surface.points->XRef()= surface.triangle_mesh->Vertices();

		for(int i=0;i<surface.points->Size();i++){
			Initialize_Particle(i);
		}

		}
	}

	virtual void Initialize_Particle(int i)	////Initialize particle with index i
	{
		surface.points->M(i)=(real)1;
		surface.points->F(i)=VectorD::Zero();
		surface.points->V(i)=VectorD::Zero();
	}
};

#endif

