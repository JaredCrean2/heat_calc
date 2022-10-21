#include "mesh/gmiDataStructure.h"
#include <stdexcept>
#include <iostream>


namespace mesh_gmi {

bool hasDownward(const GMITopo& topo, const GMIEntity& entity, int idx_down)
{
  for (int i=0; i < entity.countDownward(); ++i)
    if (entity.getDownwardIndex(i) == idx_down)
      return true;

  return false;
}

bool hasUpward(const GMITopo& topo, const GMIEntity& entity, int idx_up)
{
  for (int i=0; i < entity.countUpward(); ++i)
    if (entity.getUpwardIndex(i) == idx_up)
      return true;

  return false;
}


void verify(const GMITopo& topo)
{
  for (int dim=0; dim <= 3; ++dim)
  {
    for (int i=0; i < topo.getNumEntities(dim); ++i)
    {
      auto e = topo.getEntityByIndex(dim, i);

      // check downward
      if (dim > 0)
      {
        for (int j=0; j < e.countDownward(); ++j)
        {
          int idx_j = e.getDownwardIndex(j);
          auto e_j = topo.getEntityByIndex(dim-1, idx_j);
          if (!hasUpward(topo, e_j, i))
            throw std::runtime_error(std::string("GMIEntity with dimension ") + std::to_string(dim) 
                                     + ", index " + std::to_string(i) 
                                     + " has downward adjacency with index " + std::to_string(idx_j) +
                                     ", but it does not have the corresponding upward adjacency");
        }
      }

      if (dim < 3)
      {
        for (int j=0; j < e.countUpward(); ++j)
        {
          int idx_j = e.getUpwardIndex(j);
          auto e_j = topo.getEntityByIndex(dim+1, idx_j);
          if (!hasDownward(topo, e_j, i))
            throw std::runtime_error(std::string("GMIEntity with dimension ") + std::to_string(dim) 
                                    + ", index " + std::to_string(i) 
                                    + " has upward adjacency with index " + std::to_string(idx_j) +
                                    ", but it does not have the corresponding downward adjacency");
        }
      }
    }
  }
}


//-----------------------------------------------------------------------------
// GMI iterface
struct GMITopoModel
{
  struct gmi_model_ops const* ops;
  int n[4];
  std::shared_ptr<GMITopo> topo;
};

//TODO: use GMIEntity->getIndex() instead of getTag()


struct gmi_iter* gmi_begin(struct gmi_model* m, int dim) 
{ 
  GMITopoModel* model = reinterpret_cast<GMITopoModel*>(m);
  return reinterpret_cast<struct gmi_iter*>(model->topo->gmiBegin(dim));
}

struct gmi_ent* gmi_next(struct gmi_model* m, struct gmi_iter* i)
{ 
  GMITopoModel* model = reinterpret_cast<GMITopoModel*>(m);
  impl::GMIIterator* it = reinterpret_cast<impl::GMIIterator*>(i);
  return reinterpret_cast<struct gmi_ent*>(model->topo->gmiNext(it));
}

void gmi_end(struct gmi_model* m, struct gmi_iter* i)
{
  GMITopoModel* model = reinterpret_cast<GMITopoModel*>(m);
  impl::GMIIterator* it = reinterpret_cast<impl::GMIIterator*>(i);
  model->topo->gmiEnd(it);
}

int gmi_dim(struct gmi_model* m, struct gmi_ent* e)
{
  GMIEntity* ent = reinterpret_cast<GMIEntity*>(e);
  return ent->getDim();
}

int gmi_tag(struct gmi_model* m, struct gmi_ent* e)
{
  GMIEntity* ent = reinterpret_cast<GMIEntity*>(e);
  return ent->getTag();
}

struct gmi_ent* gmi_find(struct gmi_model* m, int dim, int tag)
{
  GMITopoModel* model = reinterpret_cast<GMITopoModel*>(m);
  return reinterpret_cast<struct gmi_ent*>(model->topo->gmiFind(dim, tag));
}

struct gmi_set* gmi_adjacent(struct gmi_model* m, struct gmi_ent* e, int dim)
{
  GMITopoModel* model = reinterpret_cast<GMITopoModel*>(m);
  GMIEntity* ent = reinterpret_cast<GMIEntity*>(e);

  int e_dim = ent->getDim();
  if (dim == e_dim + 1)
  {
    struct gmi_set* upward = gmi_make_set(ent->countUpward());
    upward->n = ent->countUpward();
    for (int i=0; i < ent->countUpward(); ++i)
    {
      int idx = ent->getUpwardIndex(i);
      upward->e[i] = reinterpret_cast<struct gmi_ent*>(&(model->topo->getEntityByIndex(dim, idx)));
      //upward->e[i] = mesh_gmi::gmi_find(m, dim, tag);
    }

    return upward;
  } else if (dim == e_dim - 1)
  {
    struct gmi_set* downward = gmi_make_set(ent->countDownward());
    downward->n = ent->countDownward();
    for (int i=0; i < ent->countDownward(); ++i)
    {
      int idx = ent->getDownwardIndex(i);
      downward->e[i] = reinterpret_cast<struct gmi_ent*>(&(model->topo->getEntityByIndex(dim, idx)));
      //downward->e[i] = mesh_gmi::gmi_find(m, dim, tag);
    }

    return downward;
  } else
    throw std::runtime_error("only one level of geometric adjacencies supported");
}

void gmi_destroy(struct gmi_model* m)
{
  GMITopoModel* model = reinterpret_cast<GMITopoModel*>(m);
  delete model;
}

struct gmi_model_ops mesh_gmi_model_ops
{
  .begin    = &gmi_begin,
  .next     = &gmi_next,
  .end      = &gmi_end,
  .dim      = &gmi_dim,
  .tag      = &gmi_tag,
  .find     = &gmi_find,
  .adjacent = &gmi_adjacent,
  .destroy  = &gmi_destroy
};

// Note: the gmi library is attempting to create a C++ 
//       class with virtual member functions without the language support
//       necessary to do that is a sane way.
//       As a result, this struct must be of the form:
//       struct Foo
//       {
//          members of gmi_model in the exact same order they appear in gmi.h
//          any additional members
//       }
//       Note: this is actually how all gmi models are constructed, even inside
//       Pumi.  See, for example, how gmi_null.c::create() callocs a struct of
//       size gmi_base, then returns the address of gmi_base::model, which is 
//       a gmi_model.  This only works because gmi_base::model is the first data
//       member of a C struct, so its address is the same as the enclosing struct.
//       All of this is to say: gmi_model* objects are not really gmi_model* objects
//       in all of gmi, they are just some related class.  So what we are doing here
//       is no worse than what is intended.

struct gmi_model* createGMITopo(std::shared_ptr<GMITopo> topo) 
{ 
  auto model = new GMITopoModel; 
  model->ops = &mesh_gmi_model_ops;
  topo->setEntityCountsArray(model->n);
  model->topo = topo;

  return reinterpret_cast<struct gmi_model*>(model);
}

struct gmi_model* createGMITopo() 
{
  auto topo = std::make_shared<GMITopo>();
  return createGMITopo(topo);
}

}