#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>

using std::vector;
using std::string;
using std::pair;
using std::size_t;

class CalcNode {
friend class CalcMesh;

protected:
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;

public:
	CalcNode() : x(0.0), y(0.0), z(0.0), vx(0.0), vy(0.0), vz(0.0) {}

	CalcNode(double x, double y, double z, double vx, double vy, double vz)
		: x(x), y(y), z(z), vx(vx), vy(vy), vz(vz) {}

	void move(double tau) {
		x += vx * tau;
		y += vy * tau;
		z += vz * tau;
	}
};

class Element {
friend class CalcMesh;

protected:
	unsigned long nodesIds[4];
};

class CalcMesh {
protected:
	vector<CalcNode> nodes;
	vector<Element> elements;
	string name;
	double t = 0;
	vector<double> (*velField)(double x, double y, double z, double t);

public:
	CalcMesh(const vector<double> &nodesCoords,
	const vector<size_t> &tetrsPoints,
	vector<double> velField(double x, double y, double z, double t),
	const string &name) {
		this->name = name;
		this->velField = velField;
		nodes.resize(nodesCoords.size() / 3);
		for(unsigned int i = 0; i < nodesCoords.size() / 3; i++) {
			double X = nodesCoords[i*3];
			double Y = nodesCoords[i*3 + 1];
			double Z = nodesCoords[i*3 + 2];
			vector<double> v = velField(X, Y, Z, 0);
			nodes[i] = CalcNode(X, Y, Z, v[0], v[1], v[2]);
		}

		elements.resize(tetrsPoints.size() / 4);
		for(unsigned int i = 0; i < tetrsPoints.size() / 4; i++) {
			elements[i].nodesIds[0] = tetrsPoints[i*4] - 1;
			elements[i].nodesIds[1] = tetrsPoints[i*4 + 1] - 1;
			elements[i].nodesIds[2] = tetrsPoints[i*4 + 2] - 1;
			elements[i].nodesIds[3] = tetrsPoints[i*4 + 3] - 1;
		}
	}

	void doTimeStep(double tau) {
		for(unsigned int i = 0; i < nodes.size(); i++) {
			nodes[i].move(tau);
			vector<double> v = velField(nodes[i].x, nodes[i].y, nodes[i].z, t);
			nodes[i].vx = v[0];
			nodes[i].vy = v[1];
			nodes[i].vz = v[2];
		}
		t += tau;
	}

	void snapshot(unsigned int snap_number) {
		// Сетка в терминах VTK
		vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
		// Точки сетки в терминах VTK
		vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

		// Векторное поле на точках сетки
		auto vel = vtkSmartPointer<vtkDoubleArray>::New();
		vel->SetName("velocity");
		vel->SetNumberOfComponents(3);

		// Обходим все точки нашей расчётной сетки
		for(unsigned int i = 0; i < nodes.size(); i++) {
			// Вставляем новую точку в сетку VTK-снапшота
			dumpPoints->InsertNextPoint(nodes[i].x, nodes[i].y, nodes[i].z);

			// Добавляем значение векторного поля в этой точке
			double _vel[3] = {nodes[i].vx, nodes[i].vy, nodes[i].vz};
			vel->InsertNextTuple(_vel);
		}

		// Грузим точки в сетку
		unstructuredGrid->SetPoints(dumpPoints);

		// Присоединяем векторное и скалярное поля к точкам
		unstructuredGrid->GetPointData()->AddArray(vel);

		// А теперь пишем, как наши точки объединены в тетраэдры
		for(unsigned int i = 0; i < elements.size(); i++) {
			auto tetra = vtkSmartPointer<vtkTetra>::New();
			tetra->GetPointIds()->SetId(0, elements[i].nodesIds[0]);
			tetra->GetPointIds()->SetId(1, elements[i].nodesIds[1]);
			tetra->GetPointIds()->SetId(2, elements[i].nodesIds[2]);
			tetra->GetPointIds()->SetId(3, elements[i].nodesIds[3]);
			unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
		}

		string fileName = name + "-" + std::to_string(snap_number) + ".vtu";
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		writer->SetFileName(fileName.c_str());
		writer->SetInputData(unstructuredGrid);
		writer->Write();
	}
};

vector<double> field(double x, double y, double z, double t) {
	return {1.9*sin(y)*sin(t), 0, 0};
}

int main() {
	// Шаг точек по пространству
	double h = 4.0;
	// Шаг по времени
	double tau = 0.01;

	const unsigned int GMSH_TETR_CODE = 4;

	// Теперь придётся немного упороться:
	// (а) построением сетки средствами gmsh,
	// (б) извлечением данных этой сетки в свой код.
	gmsh::initialize();
	gmsh::model::add("heart");

	// Считаем STL
	try {
		gmsh::merge("../../stl/heart_pendant.stl");
	} catch(...) {
		gmsh::logger::write("Could not load STL mesh: bye!");
		gmsh::finalize();
		return -1;
	}

	// Восстановим геометрию
	double angle = 40;
	bool forceParametrizablePatches = false;
	bool includeBoundary = true;
	double curveAngle = 180;
	gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary, forceParametrizablePatches, curveAngle * M_PI / 180.);
	gmsh::model::mesh::createGeometry();

	// Зададим объём по считанной поверхности
	vector<pair<int, int> > s;
	gmsh::model::getEntities(s, 2);
	vector<int> sl;
	for(auto surf : s) sl.push_back(surf.second);
	int l = gmsh::model::geo::addSurfaceLoop(sl);
	gmsh::model::geo::addVolume({l});

	gmsh::model::geo::synchronize();

	// Зададим мелкость желаемой сетки
	int f = gmsh::model::mesh::field::add("MathEval");
	gmsh::model::mesh::field::setString(f, "F", "4");
	gmsh::model::mesh::field::setAsBackgroundMesh(f);

	// Построим сетку
	gmsh::model::mesh::generate(3);

	// Теперь извлечём из gmsh данные об узлах сетки
	vector<double> nodesCoord;
	vector<size_t> nodeTags;
	vector<double> parametricCoord;
	gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

	// И данные об элементах сетки тоже извлечём, нам среди них нужны только тетраэдры, которыми залит объём
	vector<size_t> *tetrsNodesTags = nullptr;
	vector<int> elementTypes;
	vector<vector<size_t>> elementTags;
	vector<vector<size_t>> elementNodeTags;
	gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);
	for(unsigned int i = 0; i < elementTypes.size(); i++) {
		if(elementTypes[i] != GMSH_TETR_CODE)
			continue;
		tetrsNodesTags = &elementNodeTags[i];
	}

	if(tetrsNodesTags == nullptr) {
		cout << "Can not find tetra data. Exiting." << endl;
		gmsh::finalize();
		return -2;
	}

	cout << "The model has " <<  nodeTags.size() << " nodes and " << tetrsNodesTags->size() / 4 << " tetrs." << endl;

	// На всякий случай проверим, что номера узлов идут подряд и без пробелов
	for(int i = 0; i < nodeTags.size(); ++i) {
		// Индексация в gmsh начинается с 1, а не с нуля. Ну штош, значит так.
		assert(i == nodeTags[i] - 1);
	}
	// И ещё проверим, что в тетраэдрах что-то похожее на правду лежит.
	assert(tetrsNodesTags->size() % 4 == 0);

	// TODO: неплохо бы полноценно данные сетки проверять, да

	CalcMesh mesh(nodesCoord, *tetrsNodesTags, field, "heart_pendant");

	gmsh::finalize();

	for (int n = 0; n < 700; n++) {
		mesh.snapshot(n);
		mesh.doTimeStep(0.2);
	}

	return 0;
}
