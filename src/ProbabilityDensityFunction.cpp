#include "ProbabilityDensityFunction.hpp"

namespace qahwa
{

void standardize(ProbabilityDensityFunction& pdf)
{
    double sum = 0e0;
    for(auto iter = pdf.cbegin(); iter != pdf.cend(); ++iter)
        sum += *iter;
    const double standardize_constant = sum;

    for(auto iter = pdf.begin(); iter != pdf.cend(); ++iter)
        *iter /= standardize_constant;

    return;
}

}
