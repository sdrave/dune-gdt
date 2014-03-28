// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCT_INTERFACES_HH
#define DUNE_GDT_PRODUCT_INTERFACES_HH

#include <type_traits>

#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container/interfaces.hh>

#include <dune/gdt/space/interface.hh>
#include <dune/gdt/discretefunction/default.hh>

namespace Dune {
namespace GDT {


template <class Traits>
class ProductInterface : protected Stuff::CRTPInterface<ProductInterface<Traits>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType FieldType;

  const GridViewType& grid_view() const
  {
    CHECK_CRTP(this->as_imp(*this).grid_view());
    return this->as_imp(*this).grid_view();
  }

  template <class RangeType, class SourceType>
  FieldType apply2(const RangeType& range, const SourceType& source) const
  {
    CHECK_CRTP(this->as_imp(*this).apply2(range, source));
    return this->as_imp(*this).apply2(range, source);
  }
}; // class ProductInterface


template <class Traits>
class LocalizableProductInterface : protected Stuff::CRTPInterface<LocalizableProductInterface<Traits>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::SourceType SourceType;
  typedef typename Traits::FieldType FieldType;

  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;

private:
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, SourceType>::value,
                "SourceType has to be derived from Stuff::IsLocalizableFunction!");
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, RangeType>::value,
                "RangeType has to be derived from Stuff::IsLocalizableFunction!");
  static_assert(std::is_same<typename SourceType::EntityType, EntityType>::value,
                "The EntityType of SourceType and GridViewType have to match!");
  static_assert(std::is_same<typename RangeType::EntityType, EntityType>::value,
                "The EntityType of RangeType and GridViewType have to match!");
  static_assert(std::is_same<typename SourceType::DomainFieldType, DomainFieldType>::value,
                "The DomainFieldType of SourceType and GridViewType have to match!");
  static_assert(std::is_same<typename RangeType::DomainFieldType, DomainFieldType>::value,
                "The DomainFieldType of RangeType and GridViewType have to match!");
  static_assert(SourceType::dimDomain == dimDomain, "The dimDomain of SourceType and GridViewType have to match!");
  static_assert(RangeType::dimDomain == dimDomain, "The dimDomain of RangeType and GridViewType have to match!");

public:
  const GridViewType& grid_view() const
  {
    CHECK_CRTP(this->as_imp(*this).grid_view());
    return this->as_imp(*this).grid_view();
  }

  const RangeType& range() const
  {
    CHECK_CRTP(this->as_imp(*this).range());
    return this->as_imp(*this).range();
  }

  const SourceType& source() const
  {
    CHECK_CRTP(this->as_imp(*this).source());
    return this->as_imp(*this).source();
  }

  FieldType apply2()
  {
    CHECK_CRTP(this->as_imp(*this).apply2());
    return this->as_imp(*this).apply2();
  }
}; // class LocalizableProductInterface


template <class Traits>
class AssemblableProductInterface : protected Stuff::CRTPInterface<AssemblableProductInterface<Traits>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::RangeSpaceType RangeSpaceType;
  typedef typename Traits::SourceSpaceType SourceSpaceType;
  typedef typename Traits::MatrixType MatrixType;

  typedef typename MatrixType::ScalarType FieldType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;

  typedef Stuff::LA::SparsityPatternDefault PatternType;

private:
  static_assert(std::is_base_of<SpaceInterface<typename RangeSpaceType::Traits>, RangeSpaceType>::value,
                "RangeSpaceType has to be derived from SpaceInterface!");
  static_assert(std::is_base_of<SpaceInterface<typename SourceSpaceType::Traits>, SourceSpaceType>::value,
                "SourceSpaceType has to be derived from SpaceInterface!");
  static_assert(std::is_same<typename RangeSpaceType::GridViewType, GridViewType>::value,
                "The GridViewType of RangeSpaceType and GridViewType have to match!");
  static_assert(std::is_same<typename SourceSpaceType::GridViewType, GridViewType>::value,
                "The GridViewType of SourceSpaceType and GridViewType have to match!");
  static_assert(std::is_base_of<Stuff::LA::MatrixInterface<typename MatrixType::Traits>, MatrixType>::value,
                "MatrixType has to be derived from Stuff::LA::MatrixInterface!");

public:
  static PatternType pattern(const RangeSpaceType& range_space)
  {
    return pattern(range_space, range_space);
  }

  static PatternType pattern(const RangeSpaceType& range_space, const SourceSpaceType& source_space)
  {
    return pattern(range_space, source_space, *(range_space.grid_view()));
  }

  static PatternType pattern(const RangeSpaceType& range_space, const SourceSpaceType& source_space,
                             const GridViewType& grid_view)
  {
    return derived_type::pattern(range_space, source_space, grid_view);
  }

  const GridViewType& grid_view() const
  {
    CHECK_CRTP(this->as_imp(*this).grid_view());
    return this->as_imp(*this).grid_view();
  }

  const RangeSpaceType& range_space() const
  {
    CHECK_CRTP(this->as_imp(*this).range_space());
    return this->as_imp(*this).range_space();
  }

  const SourceSpaceType& source_space() const
  {
    CHECK_CRTP(this->as_imp(*this).source_space());
    return this->as_imp(*this).source_space();
  }

  void assemble()
  {
    CHECK_AND_CALL_CRTP(this->as_imp(*this).assemble());
  }

  MatrixType& matrix()
  {
    CHECK_CRTP(this->as_imp(*this).matrix());
    return this->as_imp(*this).matrix();
  }

  const MatrixType& matrix() const
  {
    CHECK_CRTP(this->as_imp(*this).matrix());
    return this->as_imp(*this).matrix();
  }

  template <class R, class S>
  FieldType apply2(const Stuff::LA::VectorInterface<R>& range, const Stuff::LA::VectorInterface<S>& source)
  {
    typedef typename R::derived_type RangeType;
    typedef typename S::derived_type SourceType;
    assert(range.size() == matrix().rows());
    assert(source.size() == matrix().cols());
    assemble();
    auto tmp = range.copy();
    matrix().mv(static_cast<const SourceType&>(source), tmp);
    return static_cast<const RangeType&>(range).dot(tmp);
  } // ... apply2(...)

  template <class R, class S>
  FieldType apply2(const ConstDiscreteFunction<RangeSpaceType, R>& range,
                   const ConstDiscreteFunction<SourceSpaceType, S>& source)
  {
    return apply2(range.vector(), source.vector());
  }
}; // class AssemblableProductInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCT_INTERFACES_HH