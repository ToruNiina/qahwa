#include "UserDefinedFunctions.hpp"
#include "HistogramMaker.hpp"
#include "CmdParser.hpp"
#include "InputFileReader.hpp"
#include <fstream>

using namespace qahwa;

int main(int argc, char *argv[])
{
    try{

    CommandLineParser cmdline(argc, argv);
    cmdline.parse();

    switch(cmdline.mode())
    {
        case Mode::HISTOGRAM_DCD:
        {
            DCDReader dcdreader(cmdline.argv_at(2));
            dcdreader.read();
            auto trajectory = dcdreader.data().traj();
            HistogramMaker histmaker(100);
            ReactionCoordinate<double> rctcrd;
            histmaker.set_range(histmaker.find_range(trajectory, rctcrd));
            auto histogram = histmaker.make_histogram(trajectory, rctcrd);
            std::ofstream ofs("qahwa_histogram.dat");
            if(!ofs.good())
                throw std::runtime_error("file open error : qahwa_histogram.dat");
            ofs << histogram;
            break;
        }
        case Mode::HISTOGRAM:
        {
            throw std::runtime_error("not implemented yet");
        }
        case Mode::UNBIASED:
        {
            throw std::runtime_error("not implemented yet");
        }
        case Mode::WHAM:
        {
            throw std::runtime_error("not implemented yet");
        }  
        default:
            throw std::logic_error("invalid enum");
    }

    }
    catch(std::exception& excpt)
    {
        std::cerr << "Error: " << excpt.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
