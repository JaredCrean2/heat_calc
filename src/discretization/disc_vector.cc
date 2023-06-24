
#include "discretization/disc_vector.h"
#include "mesh/dirichlet_update_map.h"
#include "mesh/mesh.h"
#include "parallel_exchange.h"

DiscVector::DiscVector(DiscPtr disc) : 
  m_disc(disc),
  m_vec(boost::extents[m_disc->getMesh()->getNumDofs()]),
  m_array(m_disc->getMesh()->getNumVolumeGroups())
{
  for (int i=0; i < m_disc->getMesh()->getNumVolumeGroups(); ++i)
  {
    auto& vol_group = m_disc->getMesh()->getElements(i);
    m_array[i].resize(
      boost::extents[vol_group.getNumElems()][vol_group.getNumSolPtsPerElement()]);
  }
}

DiscVector& DiscVector::operator=(const DiscVector& other)
{
  assertAlways(getDisc() == other.getDisc(), "Cannot copy DiscVectors that were created on different Discretizations");
  if (this == &other)
    return *this;

  std::copy(&(other.m_vec[0]), &(other.m_vec[0]) + other.m_vec.shape()[0], &(m_vec[0]));

  for (size_t i=0; i < m_array.size(); ++i)
  {
    auto& array_in = other.m_array[i];
    auto& array_out = m_array[i];
    for (int el=0; el < array_in.shape()[0]; ++el)
      for (int j=0; j < array_in.shape()[1]; ++j)
        array_out[el][j] = array_in[el][j];
  }

  m_is_array_modified = other.m_is_array_modified;
  m_is_vec_modified   = other.m_is_vec_modified;

  return *this;
}



void DiscVector::markVectorModified()
{ 
  if (m_is_array_modified)
    throw std::runtime_error("cannot mark vector modified when array is already modified");

  m_is_vec_modified = true;
}


void DiscVector::markArrayModified()
{ 
  if (m_is_vec_modified)
    throw std::runtime_error("cannot mark array modified when vector is already modified");

  m_is_array_modified = true;
}


// overwrites array with vector.  Does not overwrite all values because
// Dirichlet values are not present in the vector
void DiscVector::syncVectorToArray()
{
  //std::cout << "\nEntered syncVectorToArray" << std::endl;
  if (!m_is_vec_modified)
  {
    std::cerr << "Warning: called syncVectorToArray when vector has not been modified, array will not be updated" << std::endl;
    return;
  }
  
  auto dof_numbering = m_disc->getDofNumbering();

  for (int i=0; i < m_disc->getNumVolDiscs(); ++i)
  {
    VolDiscPtr vol_disc = m_disc->getVolDisc(i);
    const auto& dof_nums = dof_numbering->getDofs(i);
    auto& array_i = m_array[i];

    for (int el=0; el < vol_disc->getNumElems(); ++el)
      for (int node=0; node < vol_disc->getNumSolPtsPerElement(); ++node)
      {
        auto dof = dof_nums[el][node];
        if (dof_numbering->isDofActive(dof))
        {
          //std::cout << "dof " << dof << " -> " << el << ", " << node << ", value = " << m_vec[dof] << std::endl;
          array_i[el][node] = m_vec[dof];
        }
      }
  }

  m_is_vec_modified   = false;
  m_is_array_modified = false;
}


ArrayType<Real, 1>& DiscVector::getVector()
{
  if (m_is_array_modified)
    throw std::runtime_error("cannot get vector when array is already modified");

  return m_vec;
}

const ArrayType<Real, 1>& DiscVector::getVector() const
{
  if (m_is_array_modified)
    throw std::runtime_error("cannot get vector when array is already modified");

  return m_vec;
}


ArrayType<Real, 2>& DiscVector::getArray(const Index idx)
{
  if (m_is_vec_modified)
    throw std::runtime_error("cannot get array when vector is already modified");

  return m_array.at(idx);
}

const ArrayType<Real, 2>& DiscVector::getArray(const Index idx) const
{
  if (m_is_vec_modified)
    throw std::runtime_error("cannot get array when vector is already modified");

  return m_array.at(idx);
}

void DiscVector::updateDependentDirichletValues()
{
  using ArrayNode = Mesh::DirichletUpdateMap::ArrayNode;

  auto dirichlet_update_map = getDisc()->getDofNumbering()->getDirichletUpdateMap();

  std::vector<ArrayNode> dest_nodes;

  for (int i=0; i < dirichlet_update_map->getNumLocalSrcNodes(); ++i)
  {
    ArrayNode src_node = dirichlet_update_map->getLocalNodes(i, dest_nodes);

    double val_src = getArray(src_node.block)[src_node.el][src_node.localnode];
    for (auto& dest_node : dest_nodes)
      getArray(dest_node.block)[dest_node.el][dest_node.localnode] = val_src;
  }

  MPI_Comm comm = MPI_COMM_WORLD;
  ParallelExchange<double> exchanger(comm, 1103);
  for (int i=0; i < commSize(comm); ++i)
  {
    exchanger.getRecvBuf(dirichlet_update_map->getRecvCount(i));

    for (const ArrayNode& send_node : dirichlet_update_map->getSendDofs(i))
      exchanger.getSendBuf(i).push_back(getArray(send_node.block)[send_node.el][send_node.localnode]);
  }

  auto unpacker = [&](int rank, const std::vector<double>& buf)
  {
    for (size_t i=0; i < buf.size(); ++i)
    {
      dirichlet_update_map->getRecvDofs(rank, i, dest_nodes);
      for (const ArrayNode& dest_node : dest_nodes)
        getArray(dest_node.block)[dest_node.el][dest_node.localnode] = buf[i];
    }
  };

  exchanger.finishCommunication(unpacker);
}


void DiscVector::set(const Real val)
{
  for (auto& val_v : m_vec)
    val_v = val;

  for (int i=0; i < m_disc->getNumVolDiscs(); ++i)
  {
    auto& array_i = m_array[i];
    for (unsigned int i=0; i < array_i.shape()[0]; ++i)
      for (unsigned int j=0; j < array_i.shape()[1]; ++j)
        array_i[i][j] = val;
  }

  m_is_vec_modified   = false;
  m_is_array_modified = false;
}


DiscVectorPtr makeDiscVector(DiscPtr disc)
{
  return std::make_shared<DiscVector>(disc);
}

void printArray(DiscVectorPtr vec)
{
  Real val_sum = 0;
  auto disc = vec->getDisc();
  for (int i=0; i < disc->getNumVolDiscs(); ++i)
  {
    std::cout << "Block " << i << std::endl;
    auto vol_disc  = disc->getVolDisc(i);
    ArrayType<Real, 2> sol_coords(boost::extents[vol_disc->getNumSolPtsPerElement()][3]);

    auto& arr_i = vec->getArray(i);

    for (unsigned int el=0; el < arr_i.shape()[0]; ++el)
    {
      vol_disc->getVolumeSolCoords(el, sol_coords);
      Real el_sum = 0;
      for (unsigned int j=0; j < arr_i.shape()[1]; ++j)
      {
        //std::cout << "arr(" << el << ", " << j << ") = " << arr_i[el][j] << std::endl;

        std::string coord_str = std::string("coords = ") + std::to_string(sol_coords[j][0]) + ", " + 
                                std::to_string(sol_coords[j][1]) + ", " + std::to_string(sol_coords[j][2]) + ", ";
        int nspaces = std::max(40 - coord_str.size(), size_t(1));
        std::string space_str(std::max(nspaces, 1), ' ');
        std::cout << "el " << el << " " << coord_str << space_str << "arr(" << el << ", " << j << ") = " << arr_i[el][j] << std::endl;

        val_sum += arr_i[el][j];
        el_sum += arr_i[el][j];
      }

      std::cout << "el_sum = " << el_sum << std::endl;
    }
  }

  std::cout << "val_sum = " << val_sum << std::endl;
}


void printVector(DiscVectorPtr vec)
{
  std::cout << "printing vector form" << std::endl;
  vec->syncArrayToVector();

  auto disc = vec->getDisc();
  for (int i=0; i < disc->getNumVolDiscs(); ++i)
  {
    std::cout << "Block " << i << std::endl;
    auto vol_disc  = disc->getVolDisc(i);
    auto& dofs     = disc->getDofNumbering()->getDofs(i);
    auto& vec_vals = vec->getVector();
    std::cout << "dofs.shape = " << dofs.shape()[0] << ", " << dofs.shape()[1] << std::endl;
    std::cout << "vec_vals.shape = " << vec_vals.shape()[0] << std::endl;
    ArrayType<Real, 2> sol_coords(boost::extents[vol_disc->getNumSolPtsPerElement()][3]);

    for (int el=0; el < vol_disc->getNumElems(); ++el)
    {
      vol_disc->getVolumeSolCoords(el, sol_coords);
      for (unsigned int j=0; j < vol_disc->getNumSolPtsPerElement(); ++j)
      {
        //std::cout << "el = " << el << ", j = " << j << std::endl;
        auto dof = dofs[el][j];
        if (disc->getDofNumbering()->isDofActive(dof))
        {
          //std::cout << "dof = " << dof << std::endl;
          //std::cout << "arr(" << el << ", " << j << ") = " << vec_vals[dof] << std::endl;
          std::string coord_str = std::string("coords = ") + std::to_string(sol_coords[j][0]) + ", " + 
                                  std::to_string(sol_coords[j][1]) + ", " + std::to_string(sol_coords[j][2]) + ", ";
          int nspaces = std::max(40 - coord_str.size(), size_t(1));
          std::string space_str(nspaces, ' ');
          std::cout << "el " << el << " " << coord_str << space_str << "arr(" << el << ", " << j << ") = " << vec_vals[dof] << std::endl;

         // std::cout << "coords = " << sol_coords[j][0] << ", " << sol_coords[j][1] << ", " << sol_coords[j][2]
         //           << ",  arr(" << el << ", " << j << ") = " << vec_vals[dof] << std::endl;
        }
      }
    }
  }
}

void copyToVector(const ArrayType<Real, 1>& vec_in, DiscVectorPtr vec_out)
{
  auto& vec_out_vec = vec_out->getVector();
  for (int i=0; i < vec_out_vec.shape()[0]; ++i)
    vec_out_vec[i] = vec_in[i];

  vec_out->markVectorModified();

  vec_out->syncVectorToArray();
}

void copyFromVector(const DiscVectorPtr vec_in, ArrayType<Real, 1>& vec_out)
{
  if (!vec_in->isVectorCurrent())
    vec_in->syncArrayToVector();

  const auto& vec_in_vec = vec_in->getVector();
  for (int i=0; i < vec_in_vec.shape()[0]; ++i)
    vec_out[i] = vec_in_vec[i];
}