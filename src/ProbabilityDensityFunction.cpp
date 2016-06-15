#include "ProbabilityDensityFunction.hpp"

namespace qahwa
{

ProbabilityDensityFunction make_pdf(const Histogram<std::size_t>& hist)
{
    std::size_t sum = 0;
    for(auto iter = hist.cbegin(); iter != hist.cend(); ++iter)
        sum += *iter;

    const double standardize_constant = static_cast<double>(sum);

    ProbabilityDensityFunction retval(hist.bins(), hist.range_begin(), hist.range_end());

    auto riter = retval.begin();
    for(auto iter = hist.cbegin(); iter != hist.cend(); ++iter)
    {
        *riter = *iter / standardize_constant;
        ++riter;
    }
    assert(riter == retval.end());

    return retval;
}

ProbabilityDensityFunction standardize(const ProbabilityDensityFunction& pdf)
{
    double sum = 0e0;
    for(auto iter = pdf.cbegin(); iter != pdf.cend(); ++iter)
        sum += *iter;

    const double standardize_constant = sum;

    ProbabilityDensityFunction retval(pdf.bins(), pdf.range_begin(), pdf.range_end());

    auto riter = retval.begin();
    for(auto iter = pdf.cbegin(); iter != pdf.cend(); ++iter)
    {
        *riter = *iter / standardize_constant;
        ++riter;
    }
    assert(riter == retval.end());

    return retval;
}

}
