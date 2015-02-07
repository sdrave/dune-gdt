// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_EVALUATION_INTERFACE_HH
#define DUNE_GDT_EVALUATION_INTERFACE_HH

#include <memory>

#include <dune/common/dynmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/functions/interfaces.hh>

namespace Dune {
namespace GDT {
namespace LocalEvaluation {


/**
 *  \brief  Interface for local evaluations that depend on a codim 0 entity.
 *  \tparam numArguments  The number of local bases.
 *  \note   All evaluations have to be copyable!
 */
template< class Traits, size_t numArguments >
class Codim0Interface
{
  static_assert(AlwaysFalse< Traits >::value, "There is no interface for this numArguments!");
};


/**
 *  \brief  Interface for unary codim 0 evaluations.
 */
template< class Traits >
class Codim0Interface< Traits, 1 >
  : public Stuff::CRTPInterface< Codim0Interface< Traits, 1 >, Traits >
{
public:
  typedef typename Traits::derived_type           derived_type;
  typedef typename Traits::EntityType             EntityType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::DomainFieldType        DomainFieldType;
  static const size_t                             dimDomain = Traits::dimDomain;

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp().localFunctions(entity));
    return this->as_imp().localFunctions(entity);
  }

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam R   RangeFieldType
   *  \tparam r   dimRange of the testBase
   *  \tparam rC  dimRangeRows of the testBase
   */
  template< class R, size_t r, size_t rC >
  size_t order(const LocalfunctionTupleType& localFunctions_in,
               const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, r, rC >& testBase)
  const
  {
    CHECK_CRTP(this->as_imp().order(localFunctions_in, testBase));
    return this->as_imp().order(localFunctions_in, testBase);
  }

  /**
   *  \brief  Computes a unary codim 0 evaluation.
   *  \tparam R   RangeFieldType
   *  \tparam r   dimRange of the testBase
   *  \tparam rC  dimRangeRows of the testBase
   *  \attention ret is assumed to be zero!
   */
  template< class R, size_t r, size_t rC >
  void evaluate(const LocalfunctionTupleType& localFunctions_in,
                const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, r, rC >& testBase,
                const Dune::FieldVector< DomainFieldType, dimDomain >& localPoint,
                Dune::DynamicVector< R >& ret) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().evaluate(localFunctions_in, testBase, localPoint, ret));
  }
}; // class Codim0Interface< Traits, 1 >


/**
 *  \brief  Interface for binary codim 0 evaluations.
 **/
template< class Traits >
class Codim0Interface< Traits, 2 >
  : public Stuff::CRTPInterface< Codim0Interface< Traits, 2 >, Traits >
{
public:
  typedef typename Traits::derived_type           derived_type;
  typedef typename Traits::EntityType             EntityType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::DomainFieldType        DomainFieldType;
  static const size_t                             dimDomain = Traits::dimDomain;

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp().localFunctions(entity));
    return this->as_imp().localFunctions(entity);
  }

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam R       RangeFieldType
   *  \tparam r{T,A}  dimRange of the {testBase,ansatzBase}
   *  \tparam rC{T,A} dimRangeRows of the {testBase,ansatzBase}
   */
  template< class R, size_t rT, size_t rCT, size_t rA, size_t rCA >
  size_t order(const LocalfunctionTupleType& localFunctions,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& testBase,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatzBase) const
  {
    CHECK_CRTP(this->as_imp().order(localFunctions, testBase, ansatzBase));
    return this->as_imp().order(localFunctions, testBase, ansatzBase);
  }

  /**
   *  \brief  Computes a binary codim 0 evaluation.
   *  \tparam R         RangeFieldType
   *  \tparam r{L,T,A}  dimRange of the {localFunction,testBase,ansatzBase}
   *  \tparam rC{L,T,A} dimRangeRows of the {localFunction,testBase,ansatzBase}
   *  \attention ret is assumed to be zero!
   */
  template< class R, size_t rT, size_t rCT, size_t rA, size_t rCA >
  void evaluate(const LocalfunctionTupleType& localFunctions,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& testBase,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatzBase,
                const Dune::FieldVector< DomainFieldType, dimDomain >& localPoint,
                Dune::DynamicMatrix< R >& ret) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().evaluate(localFunctions, testBase, ansatzBase, localPoint, ret));
  }
}; // class Codim0Interface< Traits, 2 >


/**
 *  \brief  Interface for local evaluations that depend on an intersection.
 *  \tparam numArguments  The number of local bases.
 *  \note   All evaluations have to be copyable!
 */
template< class Traits, size_t numArguments >
class Codim1Interface
{
  static_assert(AlwaysFalse< Traits >::value, "There is no interface for this numArguments!");
};


/**
 *  \brief  Interface for unary codim 1 evaluations.
 */
template< class Traits >
class Codim1Interface< Traits, 1 >
    : public Stuff::CRTPInterface< Codim1Interface< Traits, 1 >, Traits >
{
public:
  typedef typename Traits::derived_type           derived_type;
  typedef typename Traits::EntityType             EntityType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::DomainFieldType        DomainFieldType;
  static const size_t                             dimDomain = Traits::dimDomain;

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp().localFunctions(entity));
    return this->as_imp().localFunctions(entity);
  }

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam R   RangeFieldType
   *  \tparam r   dimRange of the testBase
   *  \tparam rC  dimRangeRows of the testBase
   */
  template< class R, size_t r, size_t rC >
  size_t order(const LocalfunctionTupleType& localFunctions,
               const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, r, rC >& testBase)
  const
  {
    CHECK_CRTP(this->as_imp().order(localFunctions, testBase));
    return this->as_imp().order(localFunctions, testBase);
  }

  /**
   *  \brief  Computes a binary codim 1 evaluation.
   *  \tparam IntersectionType    A model of Dune::Intersection< ... >
   *  \tparam R                   RangeFieldType
   *  \tparam r                   dimRange of the testBase
   *  \tparam rC                  dimRangeRows of the testBase
   *  \attention ret is assumed to be zero!
   */
  template< class IntersectionType, class R, size_t r, size_t rC >
  void evaluate(const LocalfunctionTupleType& localFunctions,
                const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, r, rC >& testBase,
                const IntersectionType& intersection,
                const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& localPoint,
                Dune::DynamicVector< R >& ret) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().evaluate(localFunctions, testBase, intersection, localPoint, ret));
  }
}; // class Codim1Interface< Traits, 1 >


/**
 *  \brief  Interface for binary codim 1 evaluations.
 */
template< class Traits >
class Codim1Interface< Traits, 2 >
  : public Stuff::CRTPInterface< Codim1Interface< Traits, 2 >, Traits >
{
public:
  typedef typename Traits::derived_type           derived_type;
  typedef typename Traits::EntityType             EntityType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::DomainFieldType        DomainFieldType;
  static const size_t                             dimDomain = Traits::dimDomain;

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp().localFunctions(entity));
    return this->as_imp().localFunctions(entity);
  }

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam R                   RangeFieldType
   *  \tparam r{T,A}              dimRange of the {testBase,ansatzBase}
   *  \tparam rC{T,A}             dimRangeRows of the {testBase,ansatzBase}
   */
  template< class R, size_t rT, size_t rCT, size_t rA, size_t rCA >
  size_t order(const LocalfunctionTupleType& localFunctions,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& testBase,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatzBase) const
  {
    CHECK_CRTP(this->as_imp().order(localFunctions, testBase, ansatzBase));
    return this->as_imp().order(localFunctions, testBase, ansatzBase);
  }

  /**
   *  \brief  Computes a binary codim 1 evaluation.
   *  \tparam IntersectionType    A model of Dune::Intersection< ... >
   *  \tparam R                   RangeFieldType
   *  \tparam r{T,A}              dimRange of the {testBase*,ansatzBase*}
   *  \tparam rC{T,A}             dimRangeRows of the {testBase*,ansatzBase*}
   *  \attention ret is assumed to be zero!
   */
  template< class IntersectionType, class R, size_t rT, size_t rCT, size_t rA, size_t rCA >
  void evaluate(const LocalfunctionTupleType& localFunctions,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& testBase,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatzBase,
                const IntersectionType& intersection,
                const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& localPoint,
                Dune::DynamicMatrix< R >& ret) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().evaluate(localFunctions, testBase, ansatzBase, intersection, localPoint,
                                                     ret));
  }
}; // class Codim1Interface< Traits, 2 >


/**
 *  \brief  Interface for quaternary codim 1 evaluations.
 */
template< class Traits >
class Codim1Interface< Traits, 4 >
  : public Stuff::CRTPInterface< Codim1Interface< Traits, 4 >, Traits >
{
public:
  typedef typename Traits::derived_type           derived_type;
  typedef typename Traits::EntityType             EntityType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::DomainFieldType        DomainFieldType;
  static const size_t                             dimDomain = Traits::dimDomain;

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp().localFunctions(entity));
    return this->as_imp().localFunctions(entity);
  }

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam R                     RangeFieldType
   *  \tparam r{T,A}                dimRange of the {testBase*,ansatzBase*}
   *  \tparam rC{T,A}               dimRangeRows of the {testBase*,ansatzBase*}
   */
  template< class R, size_t rT, size_t rCT, size_t rA, size_t rCA >
  size_t order(const LocalfunctionTupleType localFunctionsEntity,
               const LocalfunctionTupleType localFunctionsNeighbor,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& testBaseEntity,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatzBaseEntity,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& testBaseNeighbor,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatzBaseNeighbor) const
  {
    CHECK_CRTP(this->as_imp().order(localFunctionsEntity, localFunctionsNeighbor,
                                         testBaseEntity, ansatzBaseEntity,
                                         testBaseNeighbor, ansatzBaseNeighbor));
    return this->as_imp().order(localFunctionsEntity, localFunctionsNeighbor,
                                     testBaseEntity, ansatzBaseEntity,
                                     testBaseNeighbor, ansatzBaseNeighbor);
  }

  /**
   *  \brief  Computes a quaternary codim 1 evaluation.
   *  \tparam IntersectionType      A model of Dune::Intersection< ... >
   *  \tparam R                     RangeFieldType
   *  \tparam r{T,A}                dimRange of the {testBase*,ansatzBase*}
   *  \tparam rC{T,A}               dimRangeRows of the {testBase*,ansatzBase*}
   *  \attention entityEntityRet, entityEntityRet, entityEntityRet and neighborEntityRet are assumed to be zero!
   */
  template< class IntersectionType, class R, size_t rT, size_t rCT, size_t rA, size_t rCA >
  void evaluate(const LocalfunctionTupleType& localFunctionsEntity,
                const LocalfunctionTupleType& localFunctionsNeighbor,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& testBaseEntity,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatzBaseEntity,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& testBaseNeighbor,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatzBaseNeighbor,
                const IntersectionType& intersection,
                const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& localPoint,
                Dune::DynamicMatrix< R >& entityEntityRet,
                Dune::DynamicMatrix< R >& neighborNeighborRet,
                Dune::DynamicMatrix< R >& entityNeighborRet,
                Dune::DynamicMatrix< R >& neighborEntityRet) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().evaluate(localFunctionsEntity,
                                                     localFunctionsNeighbor,
                                                     testBaseEntity, ansatzBaseEntity,
                                                     testBaseNeighbor, ansatzBaseNeighbor,
                                                     intersection,
                                                     localPoint,
                                                     entityEntityRet, neighborNeighborRet,
                                                     entityNeighborRet, neighborEntityRet));
  }
}; // class Codim1Interface< Traits, 4 >


} // namespace LocalEvaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EVALUATION_INTERFACE_HH
