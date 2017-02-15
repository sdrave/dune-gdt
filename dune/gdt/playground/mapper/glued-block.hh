// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_MAPPER_GLUED_BLOCK_HH
#define DUNE_GDT_MAPPER_GLUED_BLOCK_HH

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/type_utils.hh>

#if HAVE_DUNE_GRID_MULTISCALE
# include <dune/grid/multiscale/glued.hh>
#endif

#include <dune/gdt/spaces/interface.hh>

#include "../../mapper/interface.hh"

namespace Dune {
namespace GDT {
namespace Mapper {

#if HAVE_DUNE_GRID_MULTISCALE


template< class MacroGridType, class LocalGridType, class LocalSpaceType >
class GluedBlock;


namespace internal {


template< class MacroGridType, class LocalGridType, class LocalSpaceType >
class GluedBlockTraits
{
  static_assert(std::is_base_of< SpaceInterface< typename LocalSpaceType::Traits,
                                                 LocalSpaceType::dimDomain,
                                                 LocalSpaceType::dimRange,
                                                 LocalSpaceType::dimRangeCols >,
                                 LocalSpaceType >::value,
                "LocalSpaceType has to be derived from SpaceInterface!");
public:
  typedef GluedBlock< MacroGridType, LocalGridType, LocalSpaceType > derived_type;
  typedef typename LocalSpaceType::EntityType EntityType;
  typedef typename LocalSpaceType::MapperType::BackendType BackendType;
}; // class GluedBlockTraits


} // namespace internal


template< class MacroGridType, class LocalGridType, class LocalSpaceImp >
class GluedBlock
  : public MapperInterface< internal::GluedBlockTraits< MacroGridType, LocalGridType, LocalSpaceImp > >
{
  typedef MapperInterface< internal::GluedBlockTraits< MacroGridType, LocalGridType, LocalSpaceImp > > BaseType;
public:
  typedef internal::GluedBlockTraits< MacroGridType, LocalGridType, LocalSpaceImp > Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType  EntityType;
  typedef LocalSpaceImp                LocalSpaceType;

  typedef grid::Multiscale::Glued< MacroGridType, LocalGridType > GluedGridType;

private:
  template< class L, class E >
  class Compute
  {
    static_assert(AlwaysFalse< L >::value, "Not implemented for this kind of entity (only codim 0)!");
  };

  template< class L >
  class Compute< L, typename GluedGridType::MicroEntityType >
  {
    typedef typename GluedGridType::MicroEntityType Comdim0EntityType;

  public:
    static size_t numDofs(GluedGridType& glued_grid,
                          const std::vector< std::shared_ptr< const L > >& local_spaces,
                          const Comdim0EntityType& entity)
    {
      const auto local_entity_ptr = glued_grid.global_to_local_entity(entity);
      const auto& local_entity = *local_entity_ptr;
      const auto local_block = find_block_of(glued_grid, entity);
      const auto& local_space = *local_spaces[local_block];
      return local_space.mapper().numDofs(local_entity);
    }

    static void globalIndices(GluedGridType& glued_grid,
                              const std::vector< std::shared_ptr< const L > >& local_spaces,
                              const std::vector< size_t >& global_start_indices,
                              const Comdim0EntityType& entity,
                              Dune::DynamicVector< size_t >& ret)
    {
      const size_t block = find_block_of(glued_grid, entity);
      local_spaces[block]->mapper().globalIndices(entity, ret);
      const auto local_entity_ptr = glued_grid.global_to_local_entity(entity);
      const auto& local_entity = *local_entity_ptr;
      const size_t num_dofs = local_spaces[block]->mapper().numDofs(local_entity);
      assert(ret.size() >= num_dofs);
      for (size_t ii = 0; ii < num_dofs; ++ii)
        ret[ii] += global_start_indices[block];
    }

    static size_t mapToGlobal(GluedGridType& glued_grid,
                              const std::vector< std::shared_ptr< const L > >& local_spaces,
                              const std::vector< size_t >& global_start_indices,
                              const Comdim0EntityType& entity,
                              const size_t& localIndex)
    {
      const size_t block = find_block_of(glued_grid, entity);
      const auto local_entity_ptr = glued_grid.global_to_local_entity(entity);
      const auto& local_entity = *local_entity_ptr;
      const size_t block_local_index = local_spaces[block]->mapper().mapToGlobal(local_entity,
                                                                                 localIndex);
      return global_start_indices[block] + block_local_index;
    }

  private:
    template< class EntityType >
    static size_t find_block_of(GluedGridType& glued_grid, const EntityType& entity)
    {
      const auto global_entity_index = glued_grid.global_grid_view().indexSet().index(entity);
      return glued_grid.global_to_local_indices()[global_entity_index].first;
    }
  }; // class Compute< ..., EntityType >

public:
  GluedBlock(GluedGridType& glued_grid,
             const std::vector< std::shared_ptr< const LocalSpaceType > > local_spaces)
    : glued_grid_(glued_grid)
    , local_spaces_(local_spaces)
    , num_blocks_(local_spaces_.size())
    , size_(0)
    , max_num_dofs_(0)
  {
    if (local_spaces_.size() != glued_grid_.num_subdomains())
      DUNE_THROW(Stuff::Exceptions::shapes_do_not_match,
                            "You have to provide a local space for each subdomain of the multiscale grid!\n"
                            << "  Size of the given multiscale grid: " << glued_grid_.num_subdomains() << "\n"
                            << "  Number of local spaces given: " << local_spaces_.size());
    for (size_t bb = 0; bb < num_blocks_; ++bb) {
      max_num_dofs_ = std::max(max_num_dofs_, local_spaces_[bb]->mapper().maxNumDofs());
      global_start_indices_.push_back(size_);
      size_ += local_spaces_[bb]->mapper().size();
    }
  } // GluedBlock(...)

  size_t numBlocks() const
  {
    return num_blocks_;
  }

  size_t localSize(const size_t block) const
  {
    assert(block < num_blocks_);
    return local_spaces_[block]->mapper().size();
  }

  size_t mapToGlobal(const size_t block, const size_t localIndex) const
  {
    assert(block < num_blocks_);
    return global_start_indices_[block] + localIndex;
  }

  const BackendType& backend() const
  {
    DUNE_THROW(NotImplemented, "");
    return local_spaces_[0]->mapper().backend();
  }

  size_t size() const
  {
    return size_;
  }

  size_t maxNumDofs() const
  {
    return max_num_dofs_;
  }

  size_t numDofs(const EntityType& entity) const
  {
    DSC::TimedLogger().get("gdt.mapper.glued-block.numDofs").warn()
        << "Requiring access to global micro grid!" << std::endl;
    return Compute< LocalSpaceType, EntityType >::numDofs(glued_grid_, local_spaces_, entity);
  }

  void globalIndices(const EntityType& entity, Dune::DynamicVector< size_t >& ret) const
  {
    DSC::TimedLogger().get("gdt.mapper.glued-block.globalIndices").warn()
        << "Requiring access to global micro grid!" << std::endl;
    Compute< LocalSpaceType, EntityType >::globalIndices(glued_grid_, local_spaces_, global_start_indices_, entity, ret);
  }

  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    DSC::TimedLogger().get("gdt.mapper.glued-block.mapToGlobal").warn()
        << "Requiring access to global micro grid!" << std::endl;
    return Compute< LocalSpaceType, EntityType >::mapToGlobal(glued_grid_,
                                                              local_spaces_,
                                                              global_start_indices_,
                                                              entity,
                                                              localIndex);
  } // ... mapToGlobal(...)

private:
  GluedGridType& glued_grid_;
  std::vector< std::shared_ptr< const LocalSpaceType > > local_spaces_;
  size_t num_blocks_;
  size_t size_;
  size_t max_num_dofs_;
  std::vector< size_t > global_start_indices_;
}; // class GluedBlock


#else // HAVE_DUNE_GRID_MULTISCALE


template< class MacroGridType, class LocalGridType, class LocalSpaceType >
class GluedBlock
{
  static_assert(AlwaysFalse< MacroGridType >::value, "You are missing dune-grid-multiscale!");
};


#endif // HAVE_DUNE_GRID_MULTISCALE

} // namespace Mapper
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_MAPPER_GLUED_BLOCK_HH
