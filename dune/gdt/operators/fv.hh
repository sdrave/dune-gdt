// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_OPERATORS_FV_HH
#define DUNE_GDT_OPERATORS_FV_HH

#include <memory>
#include <type_traits>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/grid/walker/apply-on.hh>
#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/interfaces.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/functionals/base.hh>
#include <dune/gdt/local/fluxes/interfaces.hh>
#include <dune/gdt/local/fluxes/godunov.hh>
#include <dune/gdt/local/fluxes/laxfriedrichs.hh>
#include <dune/gdt/local/fluxes/kinetic.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/fv.hh>
#include <dune/gdt/local/operators/fv.hh>
#include <dune/gdt/local/operators/integrals.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/spaces/interface.hh>

#include <dune/gdt/playground/spaces/dg/dune-pdelab-wrapper.hh>

#include "interfaces.hh"
#include "base.hh"

namespace Dune {
namespace GDT {


enum class NumericalFluxes
{
  godunov,
  godunov_with_reconstruction,
  laxfriedrichs,
  laxfriedrichs_with_reconstruction,
  local_laxfriedrichs,
  local_laxfriedrichs_with_reconstruction,
  kinetic
};

// forwards
template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          SlopeLimiters slope_limiter = SlopeLimiters::minmod>
class AdvectionLaxFriedrichsOperator;

template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter = SlopeLimiters::minmod>
class AdvectionGodunovOperator;

template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter = SlopeLimiters::minmod>
class AdvectionKineticOperator;

template <class RHSEvaluationImp>
class AdvectionRHSOperator;


namespace internal {


// TODO: add static assert once type of BoundaryValueFunctionImp is decided
template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter_type>
class AdvectionTraitsBase
{
  static_assert(is_analytical_flux<AnalyticalFluxImp>::value,
                "AnalyticalFluxImp has to be derived from AnalyticalFluxInterface!");
  //  static_assert(Stuff::is_???< BoundaryValueFunctionImp >::value,
  //                "BoundaryValueFunctionImp has to be derived from ???!");
public:
  static const SlopeLimiters slope_limiter = slope_limiter_type;
  typedef AnalyticalFluxImp AnalyticalFluxType;
  typedef BoundaryValueFunctionImp BoundaryValueFunctionType;
  static const size_t dimDomain = AnalyticalFluxType::dimDomain;
  static const size_t dimRange = AnalyticalFluxType::dimRange;
  static const size_t dimRangeCols = AnalyticalFluxType::dimRangeCols;
  typedef typename AnalyticalFluxType::DomainFieldType FieldType;
  typedef typename AnalyticalFluxType::FluxJacobianRangeType JacobianType;
}; // class AdvectionTraitsBase


template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          SlopeLimiters slope_limiter_type>
class AdvectionLaxFriedrichsOperatorTraits
    : public AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter_type>
{
  static_assert(XT::Functions::is_localizable_function<LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be derived from XT::Functions::LocalizableFunctionInterface!");

  typedef AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter_type> BaseType;

public:
  using BaseType::slope_limiter;
  typedef LocalizableFunctionImp LocalizableFunctionType;
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueFunctionType;
  using BaseType::dimDomain;
  typedef typename Dune::GDT::LaxFriedrichsLocalNumericalCouplingFlux<AnalyticalFluxType,
                                                                      LocalizableFunctionType,
                                                                      dimDomain>
      NumericalCouplingFluxType;
  typedef typename Dune::GDT::
      LaxFriedrichsLocalDirichletNumericalBoundaryFlux<AnalyticalFluxType,
                                                       typename BoundaryValueFunctionType::NonparametricType,
                                                       LocalizableFunctionType,
                                                       dimDomain>
          NumericalBoundaryFluxType;
  typedef AdvectionLaxFriedrichsOperator<AnalyticalFluxImp,
                                         BoundaryValueFunctionImp,
                                         LocalizableFunctionImp,
                                         slope_limiter>
      derived_type;
}; // class AdvectionLaxFriedrichsOperatorTraits

template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter_type>
class AdvectionGodunovOperatorTraits
    : public AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter_type>
{
  typedef AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter_type> BaseType;

public:
  using BaseType::slope_limiter;
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueFunctionType;
  using BaseType::dimDomain;
  typedef typename Dune::GDT::GodunovLocalNumericalCouplingFlux<AnalyticalFluxType, dimDomain>
      NumericalCouplingFluxType;
  typedef typename Dune::GDT::GodunovLocalNumericalBoundaryFlux<AnalyticalFluxType,
                                                                typename BoundaryValueFunctionType::NonparametricType,
                                                                dimDomain>
      NumericalBoundaryFluxType;
  typedef AdvectionGodunovOperator<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter> derived_type;
}; // class AdvectionGodunovOperatorTraits

template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter_type>
class AdvectionKineticOperatorTraits
    : public AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter_type>
{
  typedef AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter_type> BaseType;

public:
  using BaseType::slope_limiter;
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueFunctionType;
  using BaseType::dimDomain;
  typedef typename Dune::GDT::KineticLocalNumericalCouplingFlux<AnalyticalFluxType, dimDomain>
      NumericalCouplingFluxType;
  typedef typename Dune::GDT::KineticLocalNumericalBoundaryFlux<AnalyticalFluxType,
                                                                typename BoundaryValueFunctionType::NonparametricType,
                                                                dimDomain>
      NumericalBoundaryFluxType;
  typedef AdvectionKineticOperator<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter> derived_type;
}; // class AdvectionKineticOperatorTraits


template <class RHSEvaluationImp>
class AdvectionRHSOperatorTraits
{
public:
  typedef AdvectionRHSOperator<RHSEvaluationImp> derived_type;
  typedef RHSEvaluationImp RHSEvaluationType;
  typedef typename RHSEvaluationImp::DomainFieldType FieldType;
  typedef typename RHSEvaluationImp::RhsJacobianRangeType JacobianType;
}; // class AdvectionRHSOperatorTraits


} // namespace internal


template <class AnalyticalFluxImp,
          class NumericalCouplingFluxImp,
          class NumericalBoundaryFluxImp,
          class BoundaryValueFunctionImp,
          class SourceImp,
          class RangeImp>
class AdvectionLocalizableDefault
    : public Dune::GDT::LocalizableOperatorBase<typename RangeImp::SpaceType::GridViewType, SourceImp, RangeImp>
{
  typedef Dune::GDT::LocalizableOperatorBase<typename RangeImp::SpaceType::GridViewType, SourceImp, RangeImp> BaseType;

  static_assert(is_analytical_flux<AnalyticalFluxImp>::value,
                "AnalyticalFluxImp has to be derived from AnalyticalFluxInterface!");
  static_assert(is_local_numerical_coupling_flux<NumericalCouplingFluxImp>::value,
                "NumericalCouplingFluxImp has to be derived from LocalNumericalCouplingFluxInterface!");
  static_assert(is_local_numerical_boundary_flux<NumericalBoundaryFluxImp>::value,
                "NumericalBoundaryFluxImp has to be derived from LocalNumericalBoundaryFluxInterface!");
  //  static_assert(std::is_base_of< ???, BoundaryValueFunctionImp >::value,
  //                "BoundaryValueFunctionImp has to be derived from ???!");
  static_assert(is_discrete_function<SourceImp>::value, "SourceImp has to be derived from DiscreteFunction!");
  static_assert(is_discrete_function<RangeImp>::value, "RangeImp has to be derived from DiscreteFunction!");

public:
  typedef AnalyticalFluxImp AnalyticalFluxType;
  typedef NumericalCouplingFluxImp NumericalCouplingFluxType;
  typedef NumericalBoundaryFluxImp NumericalBoundaryFluxType;
  typedef BoundaryValueFunctionImp BoundaryValueFunctionType;
  typedef SourceImp SourceType;
  typedef RangeImp RangeType;
  typedef typename SourceType::RangeFieldType RangeFieldType;
  typedef typename RangeType::SpaceType::GridViewType GridViewType;
  static const size_t dimDomain = GridViewType::dimension;
  typedef typename Dune::GDT::LocalCouplingFvOperator<NumericalCouplingFluxType> LocalCouplingOperatorType;
  typedef typename Dune::GDT::LocalBoundaryFvOperator<NumericalBoundaryFluxType> LocalBoundaryOperatorType;

  template <class... LocalOperatorArgTypes>
  AdvectionLocalizableDefault(const AnalyticalFluxType& analytical_flux,
                              const std::shared_ptr<BoundaryValueFunctionType>& boundary_values,
                              const SourceType& source,
                              RangeType& range,
                              LocalOperatorArgTypes&&... local_operator_args)
    : BaseType(rng.space().grid_view(), src, rng)
    , local_operator_(analytical_flux, std::forward<LocalOperatorArgTypes>(local_operator_args)...)
    , local_boundary_operator_(
          analytical_flux, boundary_values, std::forward<LocalOperatorArgTypes>(local_operator_args)...)
  {
    this->append(local_operator_, new XT::Grid::ApplyOn::InnerIntersectionsPrimally<GridViewType>());
    this->append(local_operator_, new XT::Grid::ApplyOn::PeriodicIntersectionsPrimally<GridViewType>());
    this->append(local_boundary_operator_, new XT::Grid::ApplyOn::NonPeriodicBoundaryIntersections<GridViewType>());
  }

private:
  const LocalCouplingOperatorType local_operator_;
  const LocalBoundaryOperatorType local_boundary_operator_;
}; // class AdvectionLocalizableDefault


template <class SourceImp, class RangeImp, class BoundaryValueFunctionImp, class MatrixImp, SlopeLimiters slope_limiter>
class LinearReconstructionLocalizable
    : public Dune::GDT::LocalizableOperatorBase<typename RangeImp::SpaceType::GridViewType, SourceImp, RangeImp>
{
  typedef Dune::GDT::LocalizableOperatorBase<typename RangeImp::SpaceType::GridViewType, SourceImp, RangeImp> BaseType;
  typedef LinearReconstructionLocalizable<SourceImp, RangeImp, BoundaryValueFunctionImp, MatrixImp, slope_limiter>
      ThisType;

public:
  typedef SourceImp SourceType;
  typedef RangeImp RangeType;
  typedef BoundaryValueFunctionImp BoundaryValueFunctionType;
  typedef MatrixImp MatrixType;
  typedef typename SourceType::RangeFieldType RangeFieldType;
  typedef typename RangeType::SpaceType::GridViewType GridViewType;
  static const size_t dimDomain = GridViewType::dimension;
  typedef typename Dune::GDT::LocalReconstructionFvOperator<MatrixType, BoundaryValueFunctionType, slope_limiter>
      LocalOperatorType;

  LinearReconstructionLocalizable(const SourceType& src,
                                  RangeType& rng,
                                  const MatrixType& eigenvectors,
                                  const MatrixType& eigenvectors_inverse,
                                  const std::shared_ptr<BoundaryValueFunctionType>& boundary_values)
    : BaseType(range.space().grid_view(), source, range)
    , local_operator_(eigenvectors, eigenvectors_inverse, boundary_values)
    , source_(src)
    , range_(rng)
  {
    this->append(local_operator_);
  }

private:
  const LocalOperatorType local_operator_;
  const SourceType& source_;
  RangeType& range_;
}; // class LinearReconstructionLocalizable


// TODO: remove eigen dependency of GodunovLocalNumericalCouplingFlux/GodunovLocalNumericalBoundaryFlux
#if HAVE_EIGEN

namespace internal {


template <size_t domainDim, size_t rangeDim, class MatrixType, class EigenMatrixType, class AnalyticalFluxType>
struct EigenvectorInitializer
{
  static void initialize(const AnalyticalFluxType& /*analytical_flux*/,
                         const bool /*flux_is_linear*/,
                         const bool use_linear_reconstruction,
                         std::shared_ptr<MatrixType>& /*eigenvectors*/,
                         std::shared_ptr<MatrixType>& /*eigenvectors_inverse*/)
  {
    if (use_linear_reconstruction) {
      DUNE_THROW(Dune::NotImplemented, "Linear reconstruction is only implemented in 1D!");
    }
  }
}; // struct EigenvectorInitializer<...>

template <size_t rangeDim, class MatrixType, class EigenMatrixType, class AnalyticalFluxType>
struct EigenvectorInitializer<1, rangeDim, MatrixType, EigenMatrixType, AnalyticalFluxType>
{
  static void initialize(const AnalyticalFluxType& analytical_flux,
                         const bool flux_is_linear,
                         const bool use_linear_reconstruction,
                         std::shared_ptr<MatrixType>& eigenvectors,
                         std::shared_ptr<MatrixType>& eigenvectors_inverse)
  {
    if (use_linear_reconstruction) {
      assert(flux_is_linear && "Linear reconstruction is only implemented for linear analytical fluxes!");
      // calculate matrix of eigenvectors of A, where A is the jacobian of the linear analytical flux, i.e. u_t + A*u_x
      // = 0. As the analytical flux is linear, the jacobian A is constant, so it is enough to evaluate at 0.
      ::Eigen::EigenSolver<typename EigenMatrixType::BackendType> eigen_solver(
          Dune::XT::Common::from_string<EigenMatrixType>(
              Dune::XT::Common::to_string(
                  analytical_flux.jacobian(typename AnalyticalFluxType::RangeType(0),
                                           typename AnalyticalFluxType::EntityType{},
                                           typename AnalyticalFluxType::EntityType::Geometry::LocalCoordinate(0))))
              .backend());
      assert(eigen_solver.info() == ::Eigen::Success);
      const auto eigen_eigenvectors = eigen_solver.eigenvectors();
#ifndef NDEBUG
      for (size_t ii = 0; ii < rangeDim; ++ii)
        for (size_t jj = 0; jj < rangeDim; ++jj)
          assert(eigen_eigenvectors(ii, jj).imag() < 1e-15);
#endif
      eigenvectors = std::make_shared<MatrixType>(Dune::XT::Common::from_string<MatrixType>(
          Dune::XT::Common::to_string(EigenMatrixType(eigen_eigenvectors.real()))));
      eigenvectors_inverse = std::make_shared<MatrixType>(Dune::XT::Common::from_string<MatrixType>(
          Dune::XT::Common::to_string(EigenMatrixType(eigen_eigenvectors.inverse().real()))));
    }
  }
}; // struct EigenvectorInitializer<1, ...>

template <class NumericalCouplingFluxType,
          class NumericalBoundaryFluxType,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols,
          SlopeLimiters slope_limiter>
struct AdvectionOperatorApplier
{
  template <class AnalyticalFluxType,
            class BoundaryValueFunctionType,
            class MatrixType,
            class SourceType,
            class RangeType,
            class... LocalOperatorArgTypes>
  static void apply(const AnalyticalFluxType& analytical_flux,
                    const BoundaryValueFunctionType& boundary_values,
                    const SourceType& source,
                    RangeType& range,
                    const XT::Common::Parameter param,
                    const bool use_linear_reconstruction,
                    const std::shared_ptr<MatrixType>& eigenvectors,
                    const std::shared_ptr<MatrixType>& eigenvectors_inverse,
                    LocalOperatorArgTypes&&... local_operator_args)
  {
    const std::shared_ptr<typename BoundaryValueFunctionType::NonparametricType> current_boundary_values =
        boundary_values.with_parameter(param);
    if (use_linear_reconstruction) {
      typedef DunePdelabDgProductSpaceWrapper<typename SourceType::SpaceType::GridViewType,
                                              1, // polOrder
                                              RangeFieldType,
                                              dimRange,
                                              dimRangeCols>
          DGSpaceType;
      typedef DiscreteFunction<DGSpaceType, typename SourceType::VectorType> ReconstructedDiscreteFunctionType;
      const auto dg_space = Dune::XT::Common::make_unique<DGSpaceType>(range.space().grid_view());
      const auto reconstruction =
          Dune::XT::Common::make_unique<ReconstructedDiscreteFunctionType>(*dg_space, "reconstructed");
      LinearReconstructionLocalizable<SourceType,
                                      ReconstructedDiscreteFunctionType,
                                      typename BoundaryValueFunctionType::NonparametricType,
                                      MatrixType,
                                      slope_limiter>
          reconstruction_operator(
              source, *reconstruction, *eigenvectors, *eigenvectors_inverse, current_boundary_values);
      reconstruction_operator.apply();
      AdvectionLocalizableDefault<AnalyticalFluxType,
                                  NumericalCouplingFluxType,
                                  NumericalBoundaryFluxType,
                                  typename BoundaryValueFunctionType::NonparametricType,
                                  ReconstructedDiscreteFunctionType,
                                  RangeType>
          localizable_operator(analytical_flux,
                               current_boundary_values,
                               *reconstruction,
                               range,
                               std::forward<LocalOperatorArgTypes>(local_operator_args)...);
      localizable_operator.apply();
    } else {
      AdvectionLocalizableDefault<AnalyticalFluxType,
                                  NumericalCouplingFluxType,
                                  NumericalBoundaryFluxType,
                                  typename BoundaryValueFunctionType::NonparametricType,
                                  SourceType,
                                  RangeType>
          localizable_operator(analytical_flux,
                               current_boundary_values,
                               source,
                               range,
                               std::forward<LocalOperatorArgTypes>(local_operator_args)...);
      localizable_operator.apply(true);
    }
  }
}; // struct AdvectionOperatorApplier

template <class AnalyticalFluxType,
          class BoundaryValueFunctionType,
          class SourceType,
          class RangeType,
          class MatrixType,
          class Tuple,
          class FunctionType,
          size_t... Is>
constexpr auto apply_impl(const AnalyticalFluxType& analytical_flux,
                          const BoundaryValueFunctionType& boundary_values,
                          const SourceType& source,
                          RangeType& range,
                          const XT::Common::Parameter& param,
                          const bool use_linear_reconstruction,
                          std::shared_ptr<MatrixType> eigenvectors,
                          std::shared_ptr<MatrixType> eigenvectors_inverse,
                          Tuple t,
                          FunctionType f,
                          std::index_sequence<Is...>)
{
  return f(analytical_flux,
           boundary_values,
           source,
           range,
           param,
           use_linear_reconstruction,
           eigenvectors,
           eigenvectors_inverse,
           get<Is>(t)...);
}

template <class AnalyticalFluxType,
          class BoundaryValueFunctionType,
          class SourceType,
          class RangeType,
          class MatrixType,
          class Tuple,
          class FunctionType>
constexpr auto apply(const AnalyticalFluxType& analytical_flux,
                     const BoundaryValueFunctionType& boundary_values,
                     const SourceType& source,
                     RangeType& range,
                     const XT::Common::Parameter& param,
                     const bool use_linear_reconstruction,
                     std::shared_ptr<MatrixType> eigenvectors,
                     std::shared_ptr<MatrixType> eigenvectors_inverse,
                     Tuple t,
                     FunctionType f)
{
  return apply_impl(analytical_flux,
                    boundary_values,
                    source,
                    range,
                    param,
                    use_linear_reconstruction,
                    eigenvectors,
                    eigenvectors_inverse,
                    t,
                    f,
                    std::make_index_sequence<std::tuple_size<Tuple>{}>{});
}


} // namespace internal

template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          SlopeLimiters slope_lim>
class AdvectionLaxFriedrichsOperator
    : public Dune::GDT::OperatorInterface<internal::AdvectionLaxFriedrichsOperatorTraits<AnalyticalFluxImp,
                                                                                         BoundaryValueFunctionImp,
                                                                                         LocalizableFunctionImp,
                                                                                         slope_lim>>
{
public:
  typedef internal::AdvectionLaxFriedrichsOperatorTraits<AnalyticalFluxImp,
                                                         BoundaryValueFunctionImp,
                                                         LocalizableFunctionImp,
                                                         slope_lim>
      Traits;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::BoundaryValueFunctionType BoundaryValueFunctionType;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  static const size_t dimRangeCols = Traits::dimRangeCols;
  static const SlopeLimiters slope_limiter = Traits::slope_limiter;
  typedef typename AnalyticalFluxType::RangeFieldType RangeFieldType;
  typedef typename Traits::NumericalCouplingFluxType NumericalCouplingFluxType;
  typedef typename Traits::NumericalBoundaryFluxType NumericalBoundaryFluxType;

protected:
  typedef typename Dune::XT::LA::EigenDenseMatrix<RangeFieldType> EigenMatrixType;
  typedef typename Dune::XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

public:
  AdvectionLaxFriedrichsOperator(const AnalyticalFluxType& analytical_flux,
                                 const BoundaryValueFunctionType& boundary_values,
                                 const LocalizableFunctionType& dx,
                                 const bool flux_is_linear = false,
                                 const bool use_linear_reconstruction = false,
                                 const bool use_local_laxfriedrichs_flux = false,
                                 const bool entity_geometries_equal = false)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , dx_(dx)
    , flux_is_linear_(flux_is_linear)
    , use_linear_reconstruction_(use_linear_reconstruction)
    , use_local_laxfriedrichs_flux_(use_local_laxfriedrichs_flux)
    , entity_geometries_equal_(entity_geometries_equal)
  {
    internal::EigenvectorInitializer<dimDomain, dimRange, MatrixType, EigenMatrixType, AnalyticalFluxType>::initialize(
        analytical_flux_, flux_is_linear, use_linear_reconstruction, eigenvectors_, eigenvectors_inverse_);
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter param = {}) const
  {
    internal::AdvectionOperatorApplier<NumericalCouplingFluxType,
                                       NumericalBoundaryFluxType,
                                       RangeFieldType,
                                       dimRange,
                                       dimRangeCols,
                                       slope_limiter>::apply(analytical_flux_,
                                                             boundary_values_,
                                                             source,
                                                             range,
                                                             param,
                                                             use_linear_reconstruction_,
                                                             eigenvectors_,
                                                             eigenvectors_inverse_,
                                                             dx_,
                                                             flux_is_linear_,
                                                             use_local_laxfriedrichs_flux_,
                                                             entity_geometries_equal_);
  }

private:
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueFunctionType& boundary_values_;
  const LocalizableFunctionType& dx_;
  const bool flux_is_linear_;
  const bool use_linear_reconstruction_;
  const bool use_local_laxfriedrichs_flux_;
  const bool entity_geometries_equal_;
  std::shared_ptr<MatrixType> eigenvectors_;
  std::shared_ptr<MatrixType> eigenvectors_inverse_;
}; // class AdvectionLaxFriedrichsOperator


template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_lim>
class AdvectionKineticOperator
    : public Dune::GDT::OperatorInterface<internal::AdvectionKineticOperatorTraits<AnalyticalFluxImp,
                                                                                   BoundaryValueFunctionImp,
                                                                                   slope_lim>>
{
  typedef typename internal::AdvectionKineticOperatorTraits<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_lim>
      Traits;

public:
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::BoundaryValueFunctionType BoundaryValueFunctionType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  static const size_t dimRangeCols = Traits::dimRangeCols;
  static const SlopeLimiters slope_limiter = Traits::slope_limiter;
  typedef typename AnalyticalFluxType::RangeFieldType RangeFieldType;
  typedef typename Traits::NumericalCouplingFluxType NumericalCouplingFluxType;
  typedef typename Traits::NumericalBoundaryFluxType NumericalBoundaryFluxType;

protected:
  typedef typename Dune::XT::LA::EigenDenseMatrix<RangeFieldType> EigenMatrixType;
  typedef typename Dune::XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

public:
  AdvectionKineticOperator(const AnalyticalFluxType& analytical_flux, const BoundaryValueFunctionType& boundary_values)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter param) const
  {
    internal::AdvectionOperatorApplier<NumericalCouplingFluxType,
                                       NumericalBoundaryFluxType,
                                       RangeFieldType,
                                       dimRange,
                                       dimRangeCols,
                                       slope_limiter>::apply(analytical_flux_,
                                                             boundary_values_,
                                                             source,
                                                             range,
                                                             param,
                                                             false,
                                                             eigenvectors_,
                                                             eigenvectors_inverse_,
                                                             param);
  }

  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueFunctionType& boundary_values_;
  std::shared_ptr<MatrixType> eigenvectors_;
  std::shared_ptr<MatrixType> eigenvectors_inverse_;
}; // class AdvectionKineticOperator


// TODO: 0 boundary by default, so no need to specify boundary conditions for periodic grid views
template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_lim>
class AdvectionGodunovOperator
    : public Dune::GDT::OperatorInterface<internal::AdvectionGodunovOperatorTraits<AnalyticalFluxImp,
                                                                                   BoundaryValueFunctionImp,
                                                                                   slope_lim>>
{
public:
  typedef internal::AdvectionGodunovOperatorTraits<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_lim> Traits;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::BoundaryValueFunctionType BoundaryValueFunctionType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  static const size_t dimRangeCols = Traits::dimRangeCols;
  static const SlopeLimiters slope_limiter = Traits::slope_limiter;
  typedef typename AnalyticalFluxType::RangeFieldType RangeFieldType;
  typedef typename Traits::NumericalCouplingFluxType NumericalCouplingFluxType;
  typedef typename Traits::NumericalBoundaryFluxType NumericalBoundaryFluxType;

protected:
  typedef typename Dune::XT::LA::EigenDenseMatrix<RangeFieldType> EigenMatrixType;
  typedef typename Dune::XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

public:
  AdvectionGodunovOperator(const AnalyticalFluxType& analytical_flux,
                           const BoundaryValueFunctionType& boundary_values,
                           const bool flux_is_linear = false,
                           const bool use_linear_reconstruction = false)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , flux_is_linear_(flux_is_linear)
    , use_linear_reconstruction_(use_linear_reconstruction)
  {
    internal::EigenvectorInitializer<dimDomain, dimRange, MatrixType, EigenMatrixType, AnalyticalFluxType>::initialize(
        analytical_flux_, flux_is_linear, use_linear_reconstruction, eigenvectors_, eigenvectors_inverse_);
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter param = {}) const
  {
    internal::AdvectionOperatorApplier<NumericalCouplingFluxType,
                                       NumericalBoundaryFluxType,
                                       RangeFieldType,
                                       dimRange,
                                       dimRangeCols,
                                       slope_limiter>::apply(analytical_flux_,
                                                             boundary_values_,
                                                             source,
                                                             range,
                                                             param,
                                                             use_linear_reconstruction_,
                                                             eigenvectors_,
                                                             eigenvectors_inverse_,
                                                             flux_is_linear_);
  }

private:
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueFunctionType& boundary_values_;
  const bool flux_is_linear_;
  const bool use_linear_reconstruction_;
  std::shared_ptr<MatrixType> eigenvectors_;
  std::shared_ptr<MatrixType> eigenvectors_inverse_;
}; // class AdvectionGodunovOperator

#else // HAVE_EIGEN

template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          SlopeLimiters slope_limiter>
class AdvectionLaxFriedrichsOperator
{
  static_assert(AlwaysFalse<AnalyticalFluxImp>::value, "You are missing eigen!");
};

template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter>
class AdvectionKineticOperator
{
  static_assert(AlwaysFalse<AnalyticalFluxImp>::value, "You are missing eigen!");
};

template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter>
class AdvectionGodunovOperator
{
  static_assert(AlwaysFalse<AnalyticalFluxImp>::value, "You are missing eigen!");
};

#endif // HAVE_EIGEN

template <class GridViewType, class MatrixType, class DiscreteFunctionType>
class MatrixSolveFunctor : public XT::Grid::Functor::Codim0<GridViewType>
{
  typedef typename XT::Grid::Functor::Codim0<GridViewType> BaseType;

public:
  using typename BaseType::EntityType;

  MatrixSolveFunctor(const std::vector<MatrixType>& matrices,
                     const DiscreteFunctionType& rhs,
                     DiscreteFunctionType& solution)
    : matrices_(matrices)
    , rhs_(rhs)
    , solution_(solution)
  {
  }

  virtual void apply_local(const EntityType& entity)
  {
    // get mapper
    const auto& mapper = rhs_.space().mapper();

    // copy rhs to DynamicVector
    DynamicVector<typename MatrixType::value_type> local_solution(mapper.numDofs(entity), 0.);
    DynamicVector<typename MatrixType::value_type> local_rhs(local_solution.size(), 0.);
    const auto& rhs_vector = rhs_.vector();
    auto& solution_vector = solution_.vector();
    const auto global_indices = mapper.globalIndices(entity);
    for (size_t ii = 0; ii < local_rhs.size(); ++ii)
      local_rhs[ii] = rhs_vector.get_entry(global_indices[ii]);
    // solve
    matrices_[rhs_.space().grid_view().indexSet().index(entity)].solve(local_solution, local_rhs);
    // write solution
    for (size_t ii = 0; ii < local_rhs.size(); ++ii)
      solution_vector.set_entry(global_indices[ii], local_solution[ii]);
  }

private:
  const std::vector<MatrixType>& matrices_;
  const DiscreteFunctionType& rhs_;
  DiscreteFunctionType& solution_;
};

template <class GridViewType, class MatrixType, class DiscreteFunctionType>
class MatrixApplyFunctor : public XT::Grid::Functor::Codim0<GridViewType>
{
  typedef typename XT::Grid::Functor::Codim0<GridViewType> BaseType;

public:
  using typename BaseType::EntityType;

  MatrixApplyFunctor(const std::vector<MatrixType>& matrices,
                     const DiscreteFunctionType& vector,
                     DiscreteFunctionType& result)
    : matrices_(matrices)
    , vector_(vector)
    , result_(result)
  {
  }

  virtual void apply_local(const EntityType& entity)
  {
    // get mapper
    const auto& mapper = vector_.space().mapper();

    // copy rhs to DynamicVector
    DynamicVector<typename MatrixType::value_type> local_vector(mapper.numDofs(entity), 0.);
    DynamicVector<typename MatrixType::value_type> local_result(local_vector.size(), 0.);
    const auto& vector_vector = vector_.vector();
    auto& result_vector = result_.vector();
    const auto global_indices = mapper.globalIndices(entity);
    for (size_t ii = 0; ii < local_vector.size(); ++ii)
      local_vector[ii] = vector_vector.get_entry(global_indices[ii]);
    // solve
    matrices_[vector_.space().grid_view().indexSet().index(entity)].mv(local_vector, local_result);

    // write solution
    for (size_t ii = 0; ii < local_vector.size(); ++ii)
      result_vector.set_entry(global_indices[ii], local_result[ii]);
  }

private:
  const std::vector<MatrixType>& matrices_;
  const DiscreteFunctionType& vector_;
  DiscreteFunctionType& result_;
};

template <class RHSEvaluationImp>
class AdvectionRHSOperator : public Dune::GDT::OperatorInterface<internal::AdvectionRHSOperatorTraits<RHSEvaluationImp>>
{
  static_assert(is_rhs_evaluation<RHSEvaluationImp>::value, "RHSEvaluationImp has to be derived from RHSInterface!");

public:
  typedef internal::AdvectionRHSOperatorTraits<RHSEvaluationImp> Traits;
  typedef typename Traits::RHSEvaluationType RHSEvaluationType;
  using BaseType = typename Dune::GDT::OperatorInterface<internal::AdvectionRHSOperatorTraits<RHSEvaluationImp>>;

  AdvectionRHSOperator(const RHSEvaluationType& rhs_evaluation)
    : rhs_evaluation_(rhs_evaluation)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& /*param*/) const
  {
    std::fill(range.vector().begin(), range.vector().end(), 0);
    LocalVolumeIntegralFunctional<LocalFvRhsIntegrand<RHSEvaluationType, SourceType>> local_functional(rhs_evaluation_,
                                                                                                       source);
    VectorFunctionalBase<typename RangeType::VectorType,
                         typename RangeType::SpaceType,
                         typename RangeType::SpaceType::GridViewType,
                         typename RangeType::DomainFieldType>
        functional_assembler(range.vector(), range.space());
    functional_assembler.append(local_functional);
    functional_assembler.assemble(true);
  }

  // assembles jacobian (jacobian is assumed to be zero initially)
  template <class SourceType, class MatrixTraits>
  void assemble_jacobian(XT::LA::MatrixInterface<MatrixTraits, typename SourceType::RangeFieldType>& jac,
                         const SourceType& source,
                         const XT::Common::Parameter& /*param*/ = {}) const
  {
    typedef LocalVolumeIntegralOperator<LocalFvRhsJacobianIntegrand<RHSEvaluationType, SourceType>> LocalOperatorType;
    LocalOperatorType local_operator(rhs_evaluation_, source);
    LocalVolumeTwoFormAssembler<LocalOperatorType> local_assembler(local_operator);
    SystemAssembler<typename SourceType::SpaceType> assembler(source.space());
    assembler.append(local_assembler, jac);
    assembler.assemble(true);
  }

  // assembles jacobian (jacobian is assumed to be zero initially)
  template <class SourceType, class MatrixType>
  void assemble_newton_matrix(std::vector<MatrixType>& newton_matrices,
                              const SourceType& source,
                              const XT::Common::Parameter& param) const
  {
    typedef LocalVolumeIntegralOperator<LocalFvRhsNewtonIntegrand<RHSEvaluationType, SourceType>> LocalOperatorType;
    LocalOperatorType local_operator(rhs_evaluation_, source, param);
    LocalVolumeTwoFormAssembler<LocalOperatorType> local_assembler(local_operator);
    SystemAssembler<typename SourceType::SpaceType> assembler(source.space());
    assembler.append(local_assembler, newton_matrices);
    assembler.assemble(false);
  }

  // solves with local jacobian on each entity
  template <class SourceType, class RangeType, class MatrixType>
  void solve(const std::vector<MatrixType>& newton_matrices,
             const SourceType& rhs,
             RangeType& solution,
             const XT::Common::Parameter& /*param*/ = {}) const
  {
    MatrixSolveFunctor<typename SourceType::SpaceType::GridViewType, MatrixType, SourceType> functor(
        newton_matrices, rhs, solution);
    SystemAssembler<typename SourceType::SpaceType> assembler(rhs.space());
    assembler.append(functor);
    assembler.assemble(false);
  }

  // applies local jacobian on each entity
  template <class SourceType, class RangeType, class MatrixType>
  void mv(const std::vector<MatrixType>& newton_matrices,
          const SourceType& vector,
          RangeType& result,
          const XT::Common::Parameter& /*param*/ = {}) const
  {
    MatrixApplyFunctor<typename SourceType::SpaceType::GridViewType, MatrixType, SourceType> functor(
        newton_matrices, vector, result);
    SystemAssembler<typename SourceType::SpaceType> assembler(vector.space());
    assembler.append(functor);
    assembler.assemble(false);
  }

private:
  const RHSEvaluationType& rhs_evaluation_;
}; // class AdvectionRHSOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_HH
