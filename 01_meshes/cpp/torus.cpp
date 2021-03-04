#include <set>
#include <vector>
#include <cmath>
#include <gmsh.h>

using std::vector;
using gmsh::model::geo::addPoint;
using gmsh::model::geo::addCircleArc;
using gmsh::model::geo::addLine;
using gmsh::model::geo::addCurveLoop;
using gmsh::model::geo::addPlaneSurface;
using gmsh::model::geo::addSurfaceFilling;
using gmsh::model::geo::addSurfaceLoop;
using gmsh::model::geo::addVolume;

struct point {double x, y, z;};

void rotateZ(point &X, const point &O, double angle) {
	double c = cos(angle), s = sin(angle);
	double x = O.x + c*(X.x - O.x) - s*(X.y - O.y);
	double y = O.y + s*(X.x - O.x) + c*(X.y - O.y);
	X.x = x, X.y = y;
}

void rotateY(point &X, const point &O, double angle) {
	double c = cos(angle), s = sin(angle);
	double z = O.z + c*(X.z - O.z) - s*(X.x - O.x);
	double x = O.x + s*(X.z - O.z) + c*(X.x - O.x);
	X.z = z, X.x = x;
}

int createTorus(double R, double r, double lc, int N, int M) {
	vector<vector<int>> points(N);
	vector<int> Os(N);
	for (int k = 0; k < N; k++)
		points[k].resize(M);

	const point O = {.x = 0, .y = 0, .z = 0};
	int Op = addPoint(0, 0, 0, lc);
	for (int i = 0; i < N; i++) {
		point O1 = {.x = R, .y = 0, .z = 0};
		for (int j = 0; j < M; j++) {
			point X = {.x = R + r, .y = 0, .z = 0};
			rotateZ(X, O1, 2*M_PI/M * j);
			rotateY(X, O, 2*M_PI/N * i);
			points[i][j] = addPoint(X.x, X.y, X.z, lc);
		}
		rotateY(O1, O, 2*M_PI/N * i);
		Os[i] = addPoint(O1.x, O1.y, O1.z, lc);
	}

	vector<vector<int>> edges(N), lines(N);
	for (int k = 0; k < N; k++) {
		edges[k].resize(M);
		lines[k].resize(M);
	}
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++) {
			int ni = (i + 1) % N, nj = (j + 1) % M;
			edges[i][j] = addCircleArc(points[i][j], Os[i], points[i][nj]);
			lines[i][j] = addCircleArc(points[i][j], Op, points[ni][j]);
		}
	vector<int> surfaces;
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++) {
			int ni = (i + 1) % N, nj = (j + 1) % M;
			surfaces.push_back(addSurfaceFilling({addCurveLoop({edges[i][j], lines[i][nj], -edges[ni][j], -lines[i][j]})}));
		}
	return addSurfaceLoop(surfaces);
}

int main(int argc, char **argv)
{
	gmsh::initialize();

	gmsh::model::add("torus");

	addVolume({createTorus(1.0, 0.4, 0.2, 12, 4), createTorus(1.0, 0.2, 0.2, 12, 4)});

	gmsh::model::geo::synchronize();

	gmsh::model::mesh::generate(3);

	gmsh::write("torus.msh");

	std::set<std::string> args(argv, argv + argc);
	if(!args.count("-nopopup")) gmsh::fltk::run();

	gmsh::finalize();

	return 0;
}
