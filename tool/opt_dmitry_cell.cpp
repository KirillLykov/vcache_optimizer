#include <iostream>
#include <string>
#include <fstream>
#include <array>
#include <iomanip>
//#include  <initializer_list>
#include <vcache_optimizer/vcache_optimizer.hpp>
#include <vcache_optimizer/stdout_callback.hpp>


struct mytriangle
{
	unsigned int indices[3];

	explicit mytriangle() {}

	explicit mytriangle(unsigned int const i1, unsigned int const i2, unsigned int const i3)
	{
		indices[0] = i1;
		indices[1] = i2;
		indices[2] = i3;
	}

	unsigned int operator[](int const index) const
	{
		return indices[index];
	}
};


struct myvtx
{
	float value[3];

	myvtx() {}
	myvtx(float value[3]): value{value[0], value[1], value[2]} {}
	//myvtx(std::initializer_list<float> lazy): value{lazy[0], lazy[1], lazy[2]} {}
};


struct mymesh
{
  typedef std::vector < myvtx > vertices_t;
  typedef std::vector < std::pair<float, float> > edges_t;
	typedef std::vector < mytriangle > triangles_t;
	typedef std::vector < std::array<float, 4> > dih_t;

	vertices_t vertices;
	edges_t edges;
  triangles_t triangles;
  dih_t dihs;

	void custom_operation(const std::vector< std::pair< bool, unsigned int > >& permutation)
	{
    for (unsigned int edge_index = 0; edge_index < edges.size(); ++edge_index)
    {
      std::pair<float, float>& e = edges[edge_index];
      e.first = permutation[e.first].second;
      e.second = permutation[e.second].second;
    }

    for (unsigned int dih_index = 0; dih_index < dihs.size(); ++dih_index)
    {
      auto& d = dihs[dih_index];
      for (size_t i = 0; i < 4; ++i)
        d[i] = permutation[d[i]].second;
    }
	}
};


namespace vcache_optimizer
{
  template < >
  struct mesh_traits < mymesh >
  {
    typedef std::string submesh_id_t;
    typedef myvtx vertex_t;
    typedef mytriangle triangle_t;
    typedef unsigned int vertex_index_t;
    typedef unsigned int triangle_index_t;
  };
}

mytriangle create_new_triangle(mymesh &, unsigned int const i1, unsigned int const i2, unsigned int const i3)
{
	return mytriangle(i1, i2, i3);
}

std::size_t get_num_triangles(mymesh const &m, std::string const &)
{
	return m.triangles.size();
}

std::size_t get_num_vertices(mymesh const &m, std::string const &)
{
	return m.vertices.size();
}

mytriangle get_triangle(mymesh const &m, std::string const &, unsigned int const &index)
{
	return m.triangles[index];
}

myvtx get_vertex(mymesh const &m, std::string const &, unsigned int const &index)
{
	return m.vertices[index];
}

void set_triangle(mymesh &m, std::string const &, unsigned int const &index, mytriangle const &new_triangle)
{
	m.triangles[index] = new_triangle;
}

void set_vertex(mymesh &m, std::string const &, unsigned int const &index, myvtx const &new_vertex)
{
	m.vertices[index] = new_vertex;
}

void readCellMesh(const std::string& inputFileName, mymesh& mesh)
{
  std::ifstream in(inputFileName.c_str());
  std::string line;

  if (in.good())
  {
    std::cout << "Reading file " << inputFileName << std::endl;
  }
  else
  {
    std::cout << inputFileName << ": no such file" << std::endl;
    exit(1);
  }

  size_t nparticles;
  size_t ntriang;
  size_t nbonds;
  size_t ndihedrals;
  in >> nparticles >> nbonds >> ntriang >> ndihedrals;
  size_t tmp1, tmp2, aid;
  // Atoms section
  for (size_t i = 0; i < nparticles; ++i)
  {
    myvtx pos;
    in >> tmp1 >> tmp2 >> aid >> pos.value[0] >> pos.value[1] >> pos.value[2];
    assert(tmp1 == i + 1);
    mesh.vertices.push_back(pos);

    if (aid != 1) break;
  }

  std::string buf;
  std::getline(in, buf); // skip empty line

  for (size_t i=0; i<nbonds; i++)
  {
    std::pair<float, float> edge;
    in >> tmp1 >> tmp2 >> edge.first >> edge.second;
    --edge.first; --edge.second;
    assert(edge.first < nparticles && edge.second < nparticles);
    mesh.edges.push_back(edge);
  }

  std::getline(in, buf); // skip empty line

  // Angles section --> triangles
  for (size_t i=0; i<ntriang; i++)
  {
    mytriangle tr;
    in >> tmp1 >> tmp2 >> tr.indices[0] >> tr.indices[1] >> tr.indices[2];

    --tr.indices[0]; --tr.indices[1]; --tr.indices[2];
    assert(tr.indices[0] < nparticles && tr.indices[1] < nparticles && tr.indices[2] < nparticles);
    mesh.triangles.push_back(tr);
  }

  std::getline(in, buf); // skip empty line

  // Dihedrals
  for (size_t i = 0; i < ndihedrals; ++i)
  {
    std::array<float, 4> dih;
    in >> tmp1 >> tmp2 >> dih[0] >> dih[1] >> dih[2] >> dih[3];

    --dih[0]; --dih[1]; --dih[2]; --dih[3];
    assert(dih[0] < nparticles && dih[1] < nparticles && dih[2] < nparticles && dih[3] < nparticles);
    mesh.dihs.push_back(dih);
  }

  assert(mesh.vertices.size() == nparticles && mesh.edges.size() == nbonds &&
      mesh.triangles.size() == ntriang  && mesh.dihs.size() == ndihedrals);
  in.close();
}

void writeCellMesh(const std::string& outputFileName, const mymesh& mesh)
{
  std::ofstream out(outputFileName.c_str());
  std::string line;

  size_t nparticles = mesh.vertices.size();
  size_t nbonds = mesh.edges.size();
  size_t ntriang = mesh.triangles.size();
  size_t ndihedrals = mesh.dihs.size();
  out << nparticles << "\n" << nbonds << "\n" << ntriang << "\n" << ndihedrals<< "\n\n";

  // Atoms section
  for (size_t i = 0; i < nparticles; ++i)
  {
    myvtx pos;
    out << std::setprecision(7) << i + 1 << " 1 1 " << mesh.vertices[i].value[0] << " " << mesh.vertices[i].value[1] << " " << mesh.vertices[i].value[2] << std::endl;
  }

  out << std::endl;

  for (size_t i=0; i<nbonds; i++)
  {
    out << i + 1 << " 1 "<< mesh.edges[i].first + 1 << " " << mesh.edges[i].second + 1 << std::endl;
  }

  out << std::endl;

  // Angles section --> triangles
  for (size_t i=0; i<ntriang; i++)
  {
    mytriangle tr;
    out << i + 1 << " 1 " << mesh.triangles[i].indices[0] + 1 << " " << mesh.triangles[i].indices[1] + 1
        << " " << mesh.triangles[i].indices[2] + 1 << std::endl;
  }

  out << std::endl;

  // Dihedrals
  for (size_t i = 0; i < ndihedrals; ++i)
  {
    out << i + 1 << " 1 " << mesh.dihs[i][0] + 1 << " " << mesh.dihs[i][1] + 1
        << " " << mesh.dihs[i][2] + 1 << " " << mesh.dihs[i][3] + 1 << std::endl;
  }

  out.close();

  std::cout << "The reordered mesh was saved to " << outputFileName << std::endl;
}

int main(int argc, char* argv[])
{
  if (argc != 3) {
    std::cout << "usage: ./opt-cell <input-mesh-file-name> <output-mesh-file-name>\n";
    exit(-1);
  }

  std::string inputFileName = argv[1];
  std::string outputFileName = argv[2];

	vcache_optimizer::vcache_optimizer < mymesh > optimizer;

	mymesh mesh;
	readCellMesh(inputFileName, mesh);

	vcache_optimizer::stdout_callback cb; 
	optimizer(mesh, "", cb);

	writeCellMesh(outputFileName, mesh);

	return 0;
}

