#include "HistogramMaker.hpp"
#include "CmdParser.hpp"
#include "InputFileReader.hpp"

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
            throw std::runtime_error("not implemented yet");
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
