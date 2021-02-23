#include <set>
#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <gmsh.h>

using std::vector;
using std::pair;

int main(int argc, char **argv)
{
	gmsh::initialize();

	gmsh::model::add("heart");

	try {
		gmsh::merge("heart_pendant.stl");
	} catch(...) {
		gmsh::logger::write("Could not load STL mesh: bye!");
		gmsh::finalize();
		return 1;
	}


  // Angle between two triangles above which an edge is considered as sharp:
	double angle = 0;
	bool forceParametrizablePatches = true;
	bool includeBoundary = false;
	double curveAngle = 180;

	gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary,
										forceParametrizablePatches,
										curveAngle * M_PI / 180.);
	gmsh::model::mesh::createGeometry();


	vector<pair<int, int>> s;
	gmsh::model::getEntities(s, 2);
	vector<int> sl;
	int n = 0;
	for(auto surf : s)
		sl.push_back(surf.second), n++;
	int l = gmsh::model::geo::addSurfaceLoop(sl);
	gmsh::model::geo::addVolume({l});
	std::cout << n << '\n';

	gmsh::model::geo::synchronize();

	int f = gmsh::model::mesh::field::add("MathEval");
	gmsh::model::mesh::field::setString(f, "F", "0.1");
	gmsh::model::mesh::field::setAsBackgroundMesh(f);

	gmsh::model::mesh::generate(3);

	std::set<std::string> args(argv, argv + argc);
	if(!args.count("-nopopup"))
		gmsh::fltk::run();

	gmsh::write("heart.msh");
	gmsh::finalize();
	return 0;
}
