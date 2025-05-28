#include "BoundaryCondition.h"
#include "Particles.h"

template<int d> void BoundaryConditionMacGrid<d>::Enforce_Boundary_Conditions(FaceField<real, d>& v) 
{
	for (auto p : psi_N_values) {
		int axis=p.first[0];
		int face_index=p.first[1];
		real value=p.second;
		v.face_fields[axis].array[face_index]=value;}
}

template<int d> void BoundaryConditionMacGrid<d>::Write_Psi_D_To_File_3d(std::string file_name) 
{
	Particles<d> particles;
	for (auto p : psi_D_values) {
		VectorDi cell=mac_grid.grid.Cell_Coord(p.first);
		VectorD pos=mac_grid.grid.Center(cell);
		int i=particles.Add_Element(); 
		particles.X(i)=pos;}
	particles.Write_To_File_3d(file_name);
}

template<int d> void BoundaryConditionMacGrid<d>::Write_Psi_N_To_File_3d(std::string file_name) 
{
	Particles<d> particles;
	for (auto p : psi_N_values) {
		int axis=p.first[0];
		VectorDi face=mac_grid.Face_Coord(axis,p.first[1]);
		VectorD pos=mac_grid.Face_Center(axis, face);
		int i=particles.Add_Element(); particles.X(i)=pos;}
	particles.Write_To_File_3d(file_name);
}

template class BoundaryConditionMacGrid<2>;
template class BoundaryConditionMacGrid<3>;

template<int d> void BoundaryConditionMacGridViscosity<d>::Enforce_Boundary_Conditions(FaceField<real,d>& v)
{
	for(auto p:psi_D_values){
		int axis=p.first[0];
		int face_index=p.first[1];
		real value=p.second;
		v.face_fields[axis].array[face_index]=value;}
}

template class BoundaryConditionMacGridViscosity<2>;
template class BoundaryConditionMacGridViscosity<3>;