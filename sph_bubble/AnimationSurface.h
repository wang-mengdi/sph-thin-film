#ifndef __AnimationSurface_h__
#define __AnimationSurface_h__
#include <cstdlib>
#include "PointSet.h"

template<class T> void print(const T& t) { std::cout << t << std::endl; }
template<class T, class... Args> void print(const T& t, const Args&... r) { std::cout << t << "\t"; print(r...); }
#define show(a) print(#a, a);

template<int d> class AnimationSurface : public PointSet<d>
{
	Typedef_VectorDii(d); 
	Typedef_MatrixD(d);
	using Base=PointSet<d>;
public:
	std::shared_ptr<TriangleMesh<d>> triangle_mesh;
	Array<VectorD> noisy_points;
	Array<VectorD> voronoi_centers;
	Array<Vector2i> voronoi_edges;

	void generate_noisy_points(float noise = 0.025f, float drop_rate = 0.05f, float noise_rate = 0.5f) {
		noisy_points.clear();
		for (int i = 0; i < triangle_mesh->Vertices().size(); i++) {
			auto rand_drop_point = (double)rand() / RAND_MAX;
			if (rand_drop_point <= drop_rate) {
				// drop this point
				// show(rand_drop_point);
				continue;
			}
			auto temp_point = triangle_mesh->Vertices()[i];
			auto rand_add_noisy = (double)rand() / RAND_MAX;
			if (rand_add_noisy <= noise_rate) {
				// add noise [-noise, +noise]
				for (int i = 0; i < temp_point.size(); i++) {
					temp_point[i] += ((double)rand() / RAND_MAX - 0.5) * 2 * noise;
				}
			}
			// show(temp_point);
			noisy_points.push_back(temp_point);
		}

	}

	void generate_voronoi_mesh() {
		// foreach elements, calculate center and assgin to dual_points
		voronoi_centers.clear();
		voronoi_centers.resize(triangle_mesh->Elements().size());
		show(voronoi_centers.size());
		for (int i = 0; i < triangle_mesh->Elements().size(); i++) {
			auto e = triangle_mesh->Elements()[i];
			voronoi_centers[i] = calculate_center(triangle_mesh->Vertices()[e[0]], triangle_mesh->Vertices()[e[1]], triangle_mesh->Vertices()[e[2]]);
		}

		// for i in 0:n, for j in i+1:n, check 3 edges to find adjacent faces and assign to dual_edges
		// edges in obj: l idx1 idx2
		voronoi_edges.clear();
		for (int i = 0; i < triangle_mesh->Elements().size(); i++) {
			const auto & face1 = triangle_mesh->Elements()[i];
			for (int j = i + 1; j < triangle_mesh->Elements().size(); j++) {
				const auto & face2 = triangle_mesh->Elements()[j];
				for (int k = 0; k < 3; k++) {
					int p1 = face1(k), p2 = face1((k + 1) % 3);
					for (int kk = 0; kk < 3; kk++) {
						if (p1 == face2(kk) && p2 == face2((kk + 3 - 1) % 3) ||
							p2 == face2(kk) && p1 == face2((kk + 3 - 1) % 3) 
							) {
							// face1 and face2 has a same edge
							voronoi_edges.emplace_back(i, j);
						}
					}
				}
			}
		}
		show(voronoi_edges.size());
	}


private:
	// center of the circumscribed circle
	// https://en.wikipedia.org/wiki/Circumscribed_circle#Barycentric_coordinates
	VectorD calculate_center(const VectorD & a, const VectorD & b, const VectorD & c) {
		VectorD bc = c-b;
		VectorD ca = a-c;
		VectorD ab = b-a;

		double aa = bc.dot(bc);
		double bb = ca.dot(ca);
		double cc = ab.dot(ab);

		double Barycentric_1 = aa * (bb + cc - aa);
		double Barycentric_2 = bb * (cc + aa - bb);
		double Barycentric_3 = cc * (aa + bb - cc);
		double Barycentric_sum = Barycentric_1 + Barycentric_2 + Barycentric_3;

		VectorD circCenter = (a * Barycentric_1 + b * Barycentric_2 + c * Barycentric_3 ) / Barycentric_sum;

		return circCenter;
	}



};
#endif