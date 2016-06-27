// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_GLUED_BLOCK_HH
#define DUNE_GDT_SPACES_GLUED_BLOCK_HH

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/common/timedlogging.hh>

#if HAVE_DUNE_GRID_MULTISCALE
# include <dune/grid/multiscale/glued.hh>
#endif

#include <dune/gdt/playground/mapper/glued-block.hh>

#include "../../spaces/interface.hh"

namespace Dune {
namespace GDT {
namespace Spaces {

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
//  typedef grid::Multiscale::Default< typename LocalSpaceType::GridViewType::Grid > MsGridType;
public:
  typedef GluedBlock< MacroGridType, LocalGridType, LocalSpaceType > derived_type;
  static const int                                     polOrder = LocalSpaceType::polOrder;
  typedef typename LocalSpaceType::BackendType         BackendType;
  typedef Mapper::GluedBlock<  MacroGridType, LocalGridType, LocalSpaceType > MapperType;
  typedef typename LocalSpaceType::BaseFunctionSetType BaseFunctionSetType;
  typedef typename LocalSpaceType::CommunicatorType    CommunicatorType;
private:
  typedef grid::Multiscale::Glued< MacroGridType, LocalGridType > GluedGridType;
public:
  typedef typename GluedGridType::MicroGridViewType    GridViewType;
  typedef typename LocalSpaceType::RangeFieldType      RangeFieldType;

  static const Stuff::Grid::ChoosePartView part_view_type  = LocalSpaceType::part_view_type;
  static const bool                        needs_grid_view = LocalSpaceType::needs_grid_view;
}; // class GluedBlockTraits


} // namespace internal


template< class MacroGridType, class LocalGridType, class LocalSpaceImp >
class GluedBlock
  : public SpaceInterface< internal::GluedBlockTraits< MacroGridType, LocalGridType, LocalSpaceImp >,
                           LocalSpaceImp::dimDomain,
                           LocalSpaceImp::dimRange,
                           LocalSpaceImp::dimRangeCols >
{
  typedef SpaceInterface< internal::GluedBlockTraits< MacroGridType, LocalGridType, LocalSpaceImp >,
                          LocalSpaceImp::dimDomain,
                          LocalSpaceImp::dimRange,
                          LocalSpaceImp::dimRangeCols >             BaseType;
  typedef GluedBlock< MacroGridType, LocalGridType, LocalSpaceImp > ThisType;
public:
  typedef internal::GluedBlockTraits< MacroGridType, LocalGridType, LocalSpaceImp > Traits;
  typedef typename Traits::BackendType         BackendType;
  typedef typename Traits::MapperType          MapperType;
  typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
  typedef LocalSpaceImp                        LocalSpaceType;

  typedef typename BaseType::PatternType      PatternType;
  typedef typename BaseType::GridViewType     GridViewType;
  typedef typename BaseType::EntityType       EntityType;
  typedef typename BaseType::CommunicatorType CommunicatorType;

  typedef grid::Multiscale::Glued< MacroGridType, LocalGridType > GluedGridType;

  GluedBlock(GluedGridType& glued_grid,
             const std::vector< std::shared_ptr< const LocalSpaceType > >& local_spaces)
    : glued_grid_(glued_grid)
    , local_spaces_(local_spaces)
    , mapper_(std::make_shared< MapperType >(glued_grid, local_spaces_))
    , global_micro_grid_(nullptr)
  {
    if (local_spaces_.size() != glued_grid_.num_subdomains())
      DUNE_THROW(Stuff::Exceptions::shapes_do_not_match,
                 "You have to provide a local space for each subdomain of the multiscale grid!\n"
                 << "  Size of the given multiscale grid: " << glued_grid_.num_subdomains() << "\n"
                 << "  Number of local spaces given: " << local_spaces_.size());
  } // GluedBlock(...)

  GluedBlock(const ThisType& other) = default;
  GluedBlock(ThisType&& source) = default;
  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  const GluedGridType& glued_grid() const
  {
    return glued_grid_;
  }

  const std::vector< std::shared_ptr< const LocalSpaceType > >& local_spaces() const
  {
    return local_spaces_;
  }

  const GridViewType& grid_view() const
  {
    DSC::TimedLogger().get("gdt.spaces.glued-block.grid_view").warn()
        << "Requiring access to global micro grid!" << std::endl;
    prepare_global_micro_grid();
    return *global_micro_grid_;
  } // ... grid_view(...)

  const BackendType& backend() const
  {
    DUNE_THROW(NotImplemented, "");
    return local_spaces_[0]->backend();
  }

  const MapperType& mapper() const
  {
    return *mapper_;
  }

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    DUNE_THROW(NotImplemented, "Does not work reliably yet!");
    DSC::TimedLogger().get("gdt.spaces.glued-block.base_function_set").warn()
        << "Requiring access to global micro grid!" << std::endl;
    return local_spaces_[find_block_of(entity)]->base_function_set(glued_grid_.global_to_local_entity(entity));
  }

  template< class ConstraintsType >
  void local_constraints(const EntityType& /*entity*/, ConstraintsType& /*ret*/) const
  {
    DUNE_THROW(NotImplemented, "I am not sure yet how to implement this!");
  }

  template< class G, class S, size_t d, size_t r, size_t rC >
  PatternType compute_pattern(const GridView< G >& /*local_grid_view*/,
                              const SpaceInterface< S, d, r, rC >& /*ansatz_space*/) const
  {
    DUNE_THROW(NotImplemented, "I am not sure yet how to implement this!");
    return PatternType();
  }

  CommunicatorType& communicator() const
  {
    DUNE_THROW(NotImplemented, "I am not sure yet how to implement this!");
    return local_spaces_[0]->communicator();
  }

private:
  void prepare_global_micro_grid() const
  {
    if (!global_micro_grid_)
      global_micro_grid_ = std::make_shared<GridViewType>(glued_grid_.global_grid_view());
  }

  template< class EntityType >
  size_t find_block_of(const EntityType& entity) const
  {
    prepare_global_micro_grid();
    const auto global_entity_index = glued_grid_.global_grid_view().indexSet().index(entity);
    return glued_grid_.global_to_local_indices()[global_entity_index].first;
  }

  GluedGridType& glued_grid_;
  const std::vector< std::shared_ptr< const LocalSpaceType > > local_spaces_;
  const std::shared_ptr< const MapperType > mapper_;
  mutable std::shared_ptr<GridViewType> global_micro_grid_;
}; // class GluedBlock


#else // HAVE_DUNE_GRID_MULTISCALE


template< class MacroGridType, class LocalGridType, class LocalSpaceType >
class GluedBlock
{
  static_assert(Dune::AlwaysFalse< MacroGridType >::value, "You are missing dune-grid-multiscale!");
};


#endif // HAVE_DUNE_GRID_MULTISCALE

} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_GLUED_BLOCK_HH
