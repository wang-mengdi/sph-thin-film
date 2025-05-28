#include "EdgeField.h"
#include "File.h"

template<typename T, int d>
inline void EdgeField<T, d>::Resize(const VectorDi & cell_counts)
{
	if (grid.cell_counts == cell_counts) return;
	grid.Initialize(cell_counts);
	for (int i = 0; i < d; i++)
	{
		VectorDi edge_counts = cell_counts + VectorDi::Ones() - VectorDi::Unit(d);
		edge_fields[i].Resize(edge_counts);
	}
}

template<typename T, int d>
inline void EdgeField<T, d>::Fill(const T & value)
{
	for (int i = 0; i < d; i++)
		edge_fields[i].Fill(value);
}

template<typename T, int d>
inline T & EdgeField<T, d>::operator()(const int axis, const VectorDi & coord)
{
	return edge_fields[axis](coord);
}

template<typename T, int d>
inline const T & EdgeField<T, d>::operator()(const int axis, const VectorDi & coord) const
{
	return edge_fields[axis](coord);
}

template<typename T, int d>
inline void EdgeField<T, d>::Write_Binary(std::ostream & output) const
{
	File::Write_Binary(output, grid.cell_counts);
	for(int i=0;i<d;i++) File::Write_Binary_Array(output, &edge_fields[i].array[0], (int)edge_fields[i].array.size());
}

template<typename T, int d>
inline void EdgeField<T, d>::Write_Binary(const std::string & file_name) const
{
	std::ofstream output(file_name, std::ios::binary);
	if (!output) { std::cerr << "FaceField<T, d>::Write_Binary error: cannot open file " << file_name << "\n "; exit(0); }
	Write_Binary(output);
	output.close();
}

template<typename T, int d>
inline void EdgeField<T, d>::Read_Binary(std::istream & input)
{
	VectorDi cell_counts; File::Read_Binary(input, cell_counts); Resize(cell_counts);
	for (int i = 0; i < d; i++)File::Read_Binary_Array(input, &edge_fields[i].array[0], (int)edge_fields[i].array.size());
}

template<typename T, int d>
inline void EdgeField<T, d>::Read_Binary(const std::string & file_name)
{
	std::ifstream input(file_name, std::ios::binary);
	if (!input) { std::cerr << "FaceField<T, d>::Read_Binary error: cannot open file " << file_name << "\n "; exit(0); }
	Read_Binary(input);
	input.close();

}

template<typename T, int d>
inline void EdgeField<T, d>::Write_To_File_3d(const std::string & file_name) const
{
	if constexpr (d == 3) {
		Field<Vector<T, d>, d> v; Edge_To_Cell_Conversion(*this, v);
		File::Write_Binary_To_File(file_name, v);
	}
	else {
		Field<Vector<T, d>, d> v; Edge_To_Cell_Conversion(*this, v);
		Field<Vector<T, 3>, 3> v3; VF_Dim_Conversion<T, d, 3>(v, v3);
		File::Write_Binary_To_File(file_name, v3);
	}
}

template<typename T, int d>
void EdgeField<T, d>::Edge_To_Cell_Conversion(const EdgeField<T, d>& edge_field, Field<Vector<T, d>, d>& cell_field)
{
	const Grid<d>& grid = edge_field.grid;
	cell_field.Resize(grid.cell_counts);
	if constexpr (d == 2)
	{
		for (int j = 0; j < grid.cell_counts[1]; j++) for (int i = 0; i < grid.cell_counts[0]; i++)
		{
			T vx(0), vy(0);
			vx = (T)0.5*(edge_field.edge_fields[0](Vector2i(i, j)) + edge_field.edge_fields[0](Vector2i(i, j + 1)));
			vy = (T)0.5*(edge_field.edge_fields[1](Vector2i(i, j)) + edge_field.edge_fields[1](Vector2i(i + 1, j)));
			cell_field(i, j) = Vector<T,d>(vx, vy);
		}
	}
	else
	{
		for (int k = 0; k < grid.cell_counts[2]; k++) for (int j = 0; j < grid.cell_counts[1]; j++) for (int i = 0; i < grid.cell_counts[0]; i++)
		{
			T vx(0), vy(0), vz(0);
			vx = (T)0.25*(
				edge_field.edge_fields[0](Vector3i(i, j, k)) + edge_field.edge_fields[0](Vector3i(i, j + 1, k)) + \
				edge_field.edge_fields[0](Vector3i(i, j, k + 1)) + edge_field.edge_fields[0](Vector3i(i, j + 1, k + 1))
				);
			vy = (T)0.25*(
				edge_field.edge_fields[1](Vector3i(i, j, k)) + edge_field.edge_fields[1](Vector3i(i + 1, j, k)) + \
				edge_field.edge_fields[1](Vector3i(i, j, k + 1)) + edge_field.edge_fields[1](Vector3i(i + 1, j, k + 1))
				);
			vz = (T)0.25*(
				edge_field.edge_fields[2](Vector3i(i, j, k)) + edge_field.edge_fields[2](Vector3i(i + 1, j, k)) + \
				edge_field.edge_fields[2](Vector3i(i, j + 1, k)) + edge_field.edge_fields[2](Vector3i(i + 1, j + 1, k))
				);
			cell_field(i, j, k) = Vector<T, d>(vx, vy, vz);
		}
	}
}

template class EdgeField<int, 2>;
template class EdgeField<int, 3>;
template class EdgeField<float, 2>;
template class EdgeField<float, 3>;
template class EdgeField<double, 2>;
template class EdgeField<double, 3>;
