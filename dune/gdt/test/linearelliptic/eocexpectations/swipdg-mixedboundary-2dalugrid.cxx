#include <config.h>

#if HAVE_ALUGRID

#include "swipdg-mixedboundary-2dalugrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


// polorder 1, conforming

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluConform2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 1>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L2") {
    if (test_case.num_refinements() == 1)
      return {3.89e-02, 9.55e-03};
    else
      return {4.02e-02, 1.12e-02, 2.83e-03, 6.33e-04};
  } else if (type == "H1_semi" || type == "energy") {
    if (test_case.num_refinements() == 1)
      return {2.57e-01, 1.18e-01};
    else
      return {2.69e-01, 1.39e-01, 6.87e-02, 3.08e-02};
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 2, conforming

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluConform2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 2>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L2") {
    if (test_case.num_refinements() == 1)
      return {3.59e-03, 6.36e-04};
    else
      return {3.58e-03, 6.25e-04, 1.21e-04, 2.68e-05};
  } else if (type == "H1_semi" || type == "energy") {
    if (test_case.num_refinements() == 1)
      return {4.70e-02, 1.59e-02};
    else
      return {4.81e-02, 1.79e-02, 7.19e-03, 2.85e-03};
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 1, noncoforming

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluSimplex2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 1>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L2") {
    if (test_case.num_refinements() == 1)
      return {6.46e-02, 1.75e-02};
    else
      return {6.84e-02, 2.17e-02, 5.78e-03, 1.27e-03};
  } else if (type == "H1_semi" || type == "energy") {
    if (test_case.num_refinements() == 1)
      return {3.12e-01, 1.47e-01};
    else
      return {3.28e-01, 1.76e-01, 8.72e-02, 3.89e-02};
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 2, noncoforming

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluSimplex2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 2>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L2") {
    if (test_case.num_refinements() == 1)
      return {9.67e-03, 1.50e-03};
    else
      return {9.67e-03, 1.52e-03, 2.75e-04, 5.63e-05};
  } else if (type == "H1_semi" || type == "energy") {
    if (test_case.num_refinements() == 1)
      return {9.25e-02, 2.97e-02};
    else
      return {9.33e-02, 3.19e-02, 1.19e-02, 4.66e-03};
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
