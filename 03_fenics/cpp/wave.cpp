#include <string>
#include <iostream>
#include <dolfin.h>
#include "wave.h"

using namespace dolfin;
using std::make_shared;

class DirichletBoundary : public SubDomain {
	bool inside(const Array<double> &x, bool on_boundary) const {
		return x[0] < DOLFIN_EPS || x[0] > 10000.0 - DOLFIN_EPS;
	}
};

class BoundaryField : public Expression {
	void eval(Array<double> &values, const Array<double> &r) const {
		double x = r[0], y = r[1];
		if (x < DOLFIN_EPS)
			values[0] = 10000.0;
		else
			values[0] = 0.0;
	}
};

class dU0 : public Expression {
	void eval(Array<double> &values, const Array<double> &x) const {
		values[0] = 0.0;
		values[1] = 0.0;
	}
};

int main(int argc, char *argv[]) {
	if (argc < 2) {
		std::cerr << "Usage: wave k\n";
		return 1;
	}
	auto mesh = make_shared<Mesh>(UnitSquareMesh::create({{256, 256}}, CellType::Type::triangle));
	auto V = make_shared<wave::FunctionSpace>(mesh);

	auto u0 = make_shared<BoundaryField>();
	auto boundary = make_shared<DirichletBoundary>();
	DirichletBC bc(V, u0, boundary);

	wave::BilinearForm a(V, V);
	a.k = make_shared<Constant>(std::stoi(argv[1]));
	wave::LinearForm L(V);
	L.du0 = make_shared<dU0>();

	Function u(V);
	solve(a == L, u, bc);

	File file("wave.pvd");
	file << u;

	return 0;
}
