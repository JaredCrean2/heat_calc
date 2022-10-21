#ifndef GMI_DATA_STRUCTURE_H
#define GMI_DATA_STRUCTURE_H

#include <vector>
#include <array>
#include <deque>
#include <cassert>
#include <map>
#include <memory>
#include <algorithm>
#include "error_handling.h"
#include "gmi.h"

namespace mesh_gmi {

class GMIEntity
{
  public:
    GMIEntity(int dim, int tag, int index) :
      m_dim(dim),
      m_tag(tag),
      m_index(index)
    {}

    int getDim() const { return m_dim; }

    int getTag() const {return m_tag; }

    int getIndex() const { return m_index; }

    void addUpward(const GMIEntity& entity) 
    {
      assertAlways(entity.getDim() == getDim() + 1, "Upward adjacency must have my dimension + 1");
      auto idx = entity.getIndex();
      auto it = std::find(m_upward_indices.begin(), m_upward_indices.end(), idx);
      if (it != m_upward_indices.end())
        throw std::runtime_error("duplicate upward adjacencies are not allowed");

      m_upward_indices.push_back(idx); 
    }

    int countUpward() const { return m_upward_indices.size(); }

    int getUpwardIndex(int idx) const { return m_upward_indices[idx]; }

    void addDownward(const GMIEntity& entity) 
    { 
      auto idx = entity.getIndex();
      auto it = std::find(m_downward_indices.begin(), m_downward_indices.end(), idx);
      if (it != m_downward_indices.end())
        throw std::runtime_error("duplicate downward adjacencies are not allowed");

      m_downward_indices.push_back(idx); 
    }

    int countDownward() const { return m_downward_indices.size(); }

    int getDownwardIndex(int idx) const { return m_downward_indices[idx]; }

  private:
    int m_dim;
    int m_tag;
    int m_index;
    std::vector<int> m_upward_indices;
    std::vector<int> m_downward_indices;
};

namespace impl {

struct GMIIterator
{
  explicit GMIIterator(int dim, int idx=0) :
    dim(dim),
    idx(idx)
  {}

  int dim;
  int idx;
};
}

class GMITopo
{
  private:
    using StorageType = std::deque<GMIEntity>;

  public:
    GMIEntity& createEntity(int dim, int tag)
    {
      assert(dim >= 0 && dim <= 3);
      assert(m_tag_to_index[dim].count(tag) == 0);

      int index = m_entities[dim].size();
      m_entities[dim].emplace_back(dim, tag, index);
      m_tag_to_index[dim][tag] = index;
      if (m_entity_counts)
        m_entity_counts[dim]++;

      return m_entities[dim].back();
    }

    bool hasEntityByIndex(int dim, int idx)
    {
      assert(dim >=0 && dim <= 3);
      return idx >= 0 && size_t(idx) < m_entities[dim].size();
    }

    bool hasEntityByTag(int dim, int tag)
    {
      assert(dim >= 0 && dim <= 3);
      return m_tag_to_index[dim].count(tag) > 0;
    }

    GMIEntity& getEntityByIndex(int dim, int idx)
    {
      assert(dim >=0 && dim <= 3);
      assert(idx >= 0 && unsigned(idx) < m_entities[dim].size());
      return m_entities[dim][idx];
    }

    const GMIEntity& getEntityByIndex(int dim, int idx) const
    {
      assert(dim >=0 && dim <= 3);
      assert(idx >= 0 && unsigned(idx) < m_entities[dim].size());
      return m_entities[dim][idx];
    }

    GMIEntity& getEntityByTag(int dim, int tag)
    {
      return getEntityByIndex(dim, m_tag_to_index[dim].at(tag));
    }

    const GMIEntity& getEntityByTag(int dim, int tag) const
    {
      return getEntityByIndex(dim, m_tag_to_index[dim].at(tag));
    }

    int getNumEntities(int dim) const { return m_entities[dim].size(); }

    //-------------------------------------------------------------------------
    // GMI Interface
    impl::GMIIterator* gmiBegin(int dim) { return new impl::GMIIterator(dim, 0); }

    GMIEntity* gmiNext(impl::GMIIterator* it)
    {
      if (it->idx == m_entities[it->dim].size())
        return nullptr;
      else
        return &(m_entities[it->dim][it->idx++]);
    }

    void gmiEnd(impl::GMIIterator* it) { delete it; }

    GMIEntity* gmiFind(int dim, int tag)
    {
      if (m_tag_to_index[dim].count(tag) == 0)
        throw std::runtime_error("unable to find GMIEntity of dim " + std::to_string(dim) + ", tag " + std::to_string(tag));
      
      return &(getEntityByTag(dim, tag));
    }

    void setEntityCountsArray(int* entity_counts)
    {
      m_entity_counts = entity_counts;
      for (int i=0; i < 4; ++i)
        m_entity_counts[i] = m_entities[i].size();
    }
  

  private:
    std::array<StorageType, 4> m_entities;  // index as m_entities[dim][tag]
    std::array<std::map<int, int>, 4> m_tag_to_index;
    int* m_entity_counts = nullptr;
};

void verify(const GMITopo& topo);

struct gmi_model* createGMITopo(std::shared_ptr<GMITopo> topo);

struct gmi_model* createGMITopo();


}  // namespace




#endif 