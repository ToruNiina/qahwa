#ifndef QAHWA_PROBABILITY_DENSITY_FUNCTION
#define QAHWA_PROBABILITY_DENSITY_FUNCTION
#include "Histogram.hpp"

namespace qahwa
{

using ProbabilityDensityFunction = Histogram<double>;

void standardize(ProbabilityDensityFunction& pdf);

}

#endif /* QAHWA_PROBABILITY_DENSITY_FUNCTION */
