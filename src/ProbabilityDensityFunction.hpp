#ifndef QAHWA_PROBABILITY_DENSITY_FUNCTION
#define QAHWA_PROBABILITY_DENSITY_FUNCTION
#include "Histogram.hpp"
#include "HistogramMaker.hpp"
#include "UserDefinedFunctions.hpp"

namespace qahwa
{

using ProbabilityDensityFunction = Histogram<double>;
using PotentialOfMeanForce       = Histogram<double>;

ProbabilityDensityFunction
make_pdf(const Histogram<std::size_t>& hist);

ProbabilityDensityFunction
standardize(const ProbabilityDensityFunction& pdf);

}

#endif /* QAHWA_PROBABILITY_DENSITY_FUNCTION */
