/*
 eGADA enhanced Genome Alteration Detection Algorithm

 eGADA is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 any later version.

 eGADA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with eGADA.  If not, see <http://www.gnu.org/licenses/>.

 Author:
 Yu S. Huang polyactis@gmail.com

 */
#include "BaseGADA.h"

#include <boost/program_options.hpp>  //for program options
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>

using namespace std;
using namespace boost;
namespace po = boost::program_options;

class eGADA
{
   public:
    int argc;
    char **argv;  // Somehow, can't use "char* argv[]" even though
                  // they are same.
    // otherwise "argv=_argv;" will not work for this error: incompatible types
    // in assignment of ‘char**’ to ‘char* [0]’
    string programName;
    boost::format usageDoc;
    boost::format examplesDoc;

    long M;
    long i;
    double *tn;

    int report;
    int reportIntervalDuringBE;  // how often to report progress during backward
                                 // elimination, default is 100K

    long debug;  // verbosity... set equal to 1 to see messages of SBLandBE(). 0
                 // to not see them
    double T;  // Backward elimination threshold
    // double T2=5.0; //Segment collapse to base Non-Alteration level threshold
    double BaseAmp;  // Base-level
    double a;  // SBL hyperprior parameter
    double sigma2;  // Variance observed, if negative value, it will be
                    // estimated by the mean of the differences
    // I would recommend to be estimated on all the chromosomes and as a trimmed
    // mean.
    long MinSegLen;  // Length in number of probes for a CNA segment to be
                     // called significan.
    long SelectClassifySegments;  // Classify segment into altered state (1),
                                  // otherwise 0
    long SelectEstimateBaseAmp;  // Estimate Neutral hybridization amplitude.
    double convergenceDelta;  // 1E-10 or 1E-8 seems to work well for this
                              // parameter. -- => ++ conv time
    // 1E8 better than 1E10 seems to work well for this parameter. -- => -- conv
    // time
    long maxNoOfIterations;      //=50000, //10000 is enough usually
    double convergenceMaxAlpha;  // 1E8 Maximum number of iterations to reach
                                 // convergence...
    double convergenceB;         // a number related to convergence = 1E-20

    string inputFname;
    string outputFname;

    eGADA(int _argc, char *_argv[]);  // the commandline version
    eGADA() : debug(0)
    {
        initParameters();
    }

    eGADA(long _debug) : debug(_debug)
    {
        initParameters();
    }

    virtual ~eGADA()
    {
        // cleanupMemory();	//comment it out to avoid repetitive free, in
        // python module cleanupMemory is called by the end of run().
    }

    void initParameters();
//#ifdef GADABIN

#ifndef GADABIN  // 2009-11-21 boost python module code included under if macro
                 // GADABIN (eGADA standalone) is not defined.
    void readInIntensity(boost::python::list intensity_list);
    boost::python::list run(boost::python::list intensity_list, double aAlpha,
                            double TBackElim, long MinSegLen);
#endif
#ifndef GADALib
    // stream causes error " note: synthesized method ... required
    // here" because stream is noncopyable.
    po::options_description optionDescription;  //("Allowed options")
    po::positional_options_description positionOptionDescription;
    po::variables_map optionVariableMap;

    std::ifstream inputFile;
    std::ofstream outputFile;
    boost::iostreams::filtering_streambuf<boost::iostreams::input>
        inputFilterStreamBuffer;
    boost::iostreams::filtering_streambuf<boost::iostreams::output>
        outputFilterStreamBuffer;

    // stream causes error " note: synthesized method ... required
    // here" because stream is noncopyable.
    virtual void constructOptionDescriptionStructure();
    virtual void parseCommandlineOptions();
    virtual void openOutputFile();
    virtual void openOneInputFile(
        string &inputFname,
        boost::iostreams::filtering_streambuf<boost::iostreams::input> &
            inputFilterStreamBuffer);
    virtual void closeFiles();
    void commandlineRun();
#endif
    void cleanupMemory()
    {
        free(tn);
        // free(SegState);	// SegState is not always allocated
        // with extra memory
    }
};

eGADA::eGADA(int _argc, char *_argv[]) : argc(_argc), argv(_argv)
{
    debug = 0;
    report = 0;
    programName = _argv[0];
    // strcpy(_argv[0], programName.c_str());	// does no work

    cerr << "program name is " << programName << "." << endl;

    usageDoc = boost::format("%1% -i INPUTFNAME -o OUTPUTFNAME [OPTIONS]\n") %
               programName;
    examplesDoc = boost::format(
                      "%1% -i /tmp/input.tsv.gz -o /tmp/output.gz -M 1000 "
                      "--convergenceDelta 0.01 \n") %
                  programName;
}

void eGADA::initParameters()
{
    T = 5.0;  // Backward elimination threshold
    // double T2=5.0; //Segment collapse to base Non-Alteration level threshold
    a = 0.2;  // SBL hyperprior parameter
    sigma2 = -1;  // Variance observed, if negative value, it will be estimated
                  // by the mean of the differences
    // I would recommend to be estimated on all the chromosomes and as a trimmed
    // mean.
    MinSegLen = 0;  // Length in number of probes for a CNA segment to be called
                    // significan.
    SelectClassifySegments =
        0;  // Classify segment into altered state (1), otherwise 0
    SelectEstimateBaseAmp = 1;  // Estimate Neutral hybridization amplitude.

    convergenceDelta = 1E-8;  // or 1E-8 seems to work well for this parameter.
                              // -- => ++ conv time
    // 1E8 better than 1E10 seems to work well for this parameter. -- => -- conv
    // time
    maxNoOfIterations = 50000;  // 10000 is enough usually
    convergenceMaxAlpha =
        1E8;  // 1E8 Maximum number of iterations to reach convergence...
    convergenceB = 1E-20;  // a number related to convergence = 1E-20
}

#ifndef GADABIN  // 2009-11-21 boost python module code included under if macro
                 // GADABIN (eGADA standalone) is not defined.

/*
 * Read intensity from a python list.
 */
void eGADA::readInIntensity(boost::python::list intensity_list)
{
#if defined(DEBUG)
    cerr << boost::format("# Start reading ... \n");
#endif
    long no_of_probes =
        boost::python::extract<long>(intensity_list.attr("__len__")());
    M = no_of_probes;
    tn = (double *)calloc(no_of_probes, sizeof(double));

    for (long i = 0; i < no_of_probes; i++)
    {
        tn[i] = boost::python::extract<double>(intensity_list[i]);
    }

#if defined(DEBUG)
    cerr << boost::format("# M=%1% probes in input file\n") % M;
#endif
}

boost::python::list eGADA::run(boost::python::list intensity_list, double aAlpha,
                              double TBackElim, long MinSegLen)
{
    readInIntensity(intensity_list);
    // sigma2 has to be set here. sigma2 set in eGADA::initParameters()
    // refers to the eGADA::sigma2;
    // here i suspect is the global one. Without the sentence below, the results
    // seem to be weird.
    // It used to be fine without it.
    sigma2 = -1;  // Variance observed, if negative value, it will be estimated
                  // by the mean of the differences
    BaseGADA baseGADA =
        BaseGADA(tn, M, sigma2, BaseAmp, aAlpha, TBackElim, MinSegLen, debug,
                 convergenceDelta, maxNoOfIterations, convergenceMaxAlpha,
                 convergenceB, reportIntervalDuringBE);
    baseGADA.SBLandBE();
    // K = SBLandBE(tn, M, &sigma2, aAlpha, TBackElim, MinSegLen, &Iext, &Wext,
    // debug, delta, numEMsteps, noOfBreakpointsAfterSBL, convergenceDelta,
    // maxNoOfIterations, convergenceMaxAlpha, convergenceB);
    if (debug)
    {
        std::cerr
            << boost::format(
                   " %1% breakpoints after SBL, %2% breakpoints after BE.\n") %
                   baseGADA.noOfBreakpointsAfterSBL % baseGADA.K;
    }
    if (debug)
    {
        // std::cerr<< boost::format("Backward elimination (T=%2%) and remove
        // segments that are shorter than %1% ... ") % MinSegLen % T;
        std::cerr << boost::format(
                         "Starting IextToSegLen() & IextWextToSegAmp() ... ");
    }

// BEwTandMinLen(Wext, Iext, &K, sigma2, T, MinSegLen, debug);

#if defined(DEBUG)
    cerr << boost::format("# Overall mean %1%.\n") % baseGADA.Wext[0];
    cerr << boost::format("# Sigma^2=%1%.\n") % baseGADA.sigma2;
    cerr << boost::format("# Found %1% breakpoints after SBL\n") %
                baseGADA.noOfBreakpointsAfterSBL;
#endif

// BEwTandMinLen(Wext, Iext, &K, sigma2, TBackElim, MinSegLen, debug);

#if defined(DEBUG)
    cerr << boost::format("# Kept %1% breakpoints after BE\n") % baseGADA.K;
#endif
    baseGADA.IextToSegLen();
    baseGADA.IextWextToSegAmp();

    if (debug)
    {
        cerr << boost::format("# Making segments\n");
    }

    boost::python::list return_ls;
    // fprintf(fout,"Start\tStop\tLength\tAmpl\n");
    for (i = 0; i < baseGADA.K + 1; i++)
    {
        boost::python::list d_row;
        d_row.append(baseGADA.Iext[i] + 1);
        d_row.append(baseGADA.Iext[i + 1]);
        d_row.append(baseGADA.SegLen[i]);
        d_row.append(baseGADA.SegAmp[i]);
        // cerr<< boost::format("%1% \t %2% \t %3% %4% \n")%Iext[i] % Iext[i+1]
        // % SegLen[i] % SegAmp[i];
        return_ls.append(d_row);
    }
    if (debug)
    {
        cerr << boost::format(" return %1% segments to python.\n") %
                    (baseGADA.K + 1);
    }
    // cleanupMemory();	//release tn, SegLen, SegAmp, and others
    return return_ls;
}

BOOST_PYTHON_MODULE(eGADA)
{
    using namespace boost::python;
    class_<eGADA>("eGADA").def(init<long>()).def("run", &eGADA::run);
}
#endif

#ifndef GADALib
// stream causes error " note: synthesized method ... required here"
// because stream is noncopyable.
// boost smart pointer (or cxx11) solves the problem but too much work

void eGADA::constructOptionDescriptionStructure()
{
    optionDescription.add_options()("help,h", "produce help message")(
        "TBackElim,T", po::value<double>(&T)->default_value(5.0),
        " is the backward elimination critical value for a breakpoint. i.e. "
        "minimum (mean1-mean2)/stddev difference between two adjacent "
        "segments.")("aAlpha,a", po::value<double>(&a)->default_value(0.5),
                     "is the SBL hyper prior parameter for a breakpoint. It is "
                     "the  shape parameter of the Gamma distribution. Higher "
                     "(lower) value means less (more) breakpoints expected a "
                     "priori.")("MinSegLen,M",
                                po::value<long>(&MinSegLen)->default_value(0),
                                "is the minimum size in number of probes for a "
                                "segment to be deemed significant.")(
        "BaseAmp", po::value<double>(&BaseAmp)->default_value(0.0),
        "Mean amplitude associated to the Neutral state. If not "
        "provided, and c option is used, then it is estimated as the "
        "median value of all probes hybridization values after running "
        "the algorithm. We recomend to estimate this on chromosomes "
        "that are known to have a Neutral state on most areas. In some "
        "cases this value may be known if we have been applied some "
        "normalization, preprocessing or using another sample as ref.")(
        "sigma2,s", po::value<double>(&sigma2)->default_value(-1),
        "Variance observed, if negative value, it will be estimated by the "
        "mean of the differences. "
        "I would recommend to be estimated on all the chromosomes and as a "
        "trimmed mean.")(
        "SelectClassifySegments,c",
        po::value<long>(&SelectClassifySegments)->default_value(0),
        "Classify segment into altered state (1), otherwise 0")(
        "SelectEstimateBaseAmp",
        po::value<long>(&SelectEstimateBaseAmp)->default_value(1),
        "toggle this to estimate BaseAmp from data, rather than "
        "user-supplied.")(
        "convergenceDelta",
        po::value<double>(&convergenceDelta)->default_value(1E-8),
        "a delta number controlling convergence")(
        "maxNoOfIterations",
        po::value<long>(&maxNoOfIterations)->default_value(50000),
        "maximum number of iterations for EM convergence algorithm to run "
        "before being stopped")(
        "convergenceMaxAlpha",
        po::value<double>(&convergenceMaxAlpha)->default_value(1E8),
        "one convergence related number, not sure what it does.")(
        "convergenceB", po::value<double>(&convergenceB)->default_value(1E-20),
        "one convergence related number, not sure what it does")(
        "debug,b", "toggle debug mode")("report,r", "toggle report mode")(
        "reportIntervalDuringBE",
        po::value<int>(&reportIntervalDuringBE)->default_value(100000),
        "how often to report the break point to be removed during backward "
        "elimination")("inputFname,i", po::value<string>(&inputFname),
                       "input filename, gzipped or not. could be specified as "
                       "option or positional argument."
                       "It is a single column text file with no header.")(
        "outputFname,o", po::value<string>(&outputFname), "output filename");
}

void eGADA::parseCommandlineOptions()
{
    // all positional arguments are input files.
    positionOptionDescription.add("inputFname", -1);

    po::store(po::command_line_parser(argc, argv)
                  .options(optionDescription)
                  .positional(positionOptionDescription)
                  .run(),
              optionVariableMap);

    // po::store(po::parse_command_line(argc, argv, optionDescription),
    // optionVariableMap);
    po::notify(optionVariableMap);
    if (optionVariableMap.count("help") || inputFname.empty() ||
        outputFname.empty())
    {
        cout << "Usage:" << endl << usageDoc << endl;

        cout << optionDescription << endl << endl;
        cout << "Examples:" << endl << examplesDoc << endl;
        exit(1);
    }
    if (optionVariableMap.count("debug"))
    {
        debug = 1;
    }
    else
    {
        debug = 0;
    }
    if (optionVariableMap.count("report"))
    {
        report = 1;
    }
    else
    {
        report = 0;
    }
}

void eGADA::openOneInputFile(
    string &inputFname, boost::iostreams::filtering_streambuf<
                            boost::iostreams::input> &inputFilterStreamBuffer)
{

    int inputFnameLength = inputFname.length();
    if (inputFname.substr(inputFnameLength - 3, 3) == ".gz")
    {
        inputFilterStreamBuffer.push(boost::iostreams::gzip_decompressor());
        inputFile.open(inputFname.c_str(), std::ios::in | std::ios::binary);
    }
    else
    {
        inputFile.open(inputFname.c_str(), std::ios::in);
    }
    inputFilterStreamBuffer.push(inputFile);
}

void eGADA::openOutputFile()
{
    if (!outputFname.empty())
    {
        if (debug)
        {
            std::cerr << "Open file " << outputFname << " for writing ";
        }
        int outputFnameLength = outputFname.length();
        if (outputFname.substr(outputFnameLength - 3, 3) == ".gz")
        {
            // boost::iostreams::gzip_compressor gzipCompressor;
            outputFilterStreamBuffer.push(boost::iostreams::gzip_compressor());
            // outputFilterStreamBuffer.push(boost::iostreams::base64_encoder());
            // outputFile.open(outputFname.c_str(), std::ios::out |
            // std::ios::binary);
            outputFile.open(outputFname.c_str(),
                            std::ios::out | std::ios::binary);
        }
        else
        {
            outputFile.open(outputFname.c_str(), std::ios::out);
        }
        outputFilterStreamBuffer.push(outputFile);
        // outputStream.rdbuf(&outputFilterStreamBuffer);
        //(&outputFilterStreamBuffer);
        if (debug)
        {
            std::cerr << endl;
        }
    }
    else
    {
        if (debug)
        {
            std::cerr << "Warning: Output file, " << outputFname
                      << ", is an empty string." << endl;
        }
    }
}

void eGADA::closeFiles()
{
    //
    inputFile.close();
    if (!outputFname.empty())
    {
        if (debug)
        {
            std::cerr << "closing files "
                      << "...";
        }
        // delete outputFilterStreamBuffer;
        outputFile.flush();
        // outputFile.close();	// if output is a .gz, closing it
        // here would result a truncated .gz file. don't know why.
        // maybe because the buffer or the gzip compressor filter needs to be
        // closed beforehand.
    }
}

void eGADA::commandlineRun()
{
    if (debug)
    {
        std::cerr << "Entering eGADA.commandlineRun()..." << std::endl;
    }

    constructOptionDescriptionStructure();
    parseCommandlineOptions();

    if (debug)
    {
        std::cerr << "Reading data from " << inputFname << " ... ";
    }
    openOneInputFile(inputFname, inputFilterStreamBuffer);
    std::istream inputStream(&inputFilterStreamBuffer);

    i = 0;
    M = 1000;
    tn = (double *)calloc(M, sizeof(double));
    std::string line;
    std::getline(inputStream, line);
    while (!line.empty())
    {
        tn[i++] = atof(line.c_str());
        if (i >= M)
        {
            M = M + 1000;
            tn = (double *)realloc(tn, M * sizeof(double));
        }
        std::getline(inputStream, line);
    }
    M = i;
    tn = (double *)realloc(tn, M * sizeof(double));
    if (debug)
    {
        std::cerr << M << " data points." << endl;
    }

    if (debug)
    {
        std::cerr << "Running SBLandBE ... " << endl;
    }
    BaseGADA baseGADA =
        BaseGADA(tn, M, sigma2, BaseAmp, a, T, MinSegLen, debug,
                 convergenceDelta, maxNoOfIterations, convergenceMaxAlpha,
                 convergenceB, reportIntervalDuringBE);
    baseGADA.SBLandBE();
    // K = SBLandBE(tn, M, &sigma2, a, T, MinSegLen, &Iext, &Wext, debug ,
    // delta, numEMsteps, noOfBreakpointsAfterSBL, convergenceDelta,
    // maxNoOfIterations, convergenceMaxAlpha, convergenceB);
    if (debug)
    {
        std::cerr
            << boost::format(
                   " %1% breakpoints after SBL, %2% breakpoints after BE.\n") %
                   baseGADA.noOfBreakpointsAfterSBL % baseGADA.K;
    }
    if (debug)
    {
        // std::cerr<< boost::format("Backward elimination (T=%2%) and remove
        // segments that are shorter than %1% ... ") % MinSegLen % T;
        std::cerr << boost::format(
                         "Starting IextToSegLen() & IextWextToSegAmp() ... ");
    }

    // BEwTandMinLen(Wext, Iext, &K, sigma2, T, MinSegLen, debug);

    baseGADA.IextToSegLen();
    baseGADA.IextWextToSegAmp();

    if (debug)
    {
        std::cerr << " IextToSegLen() & IextWextToSegAmp() done." << endl;
    }

    if (debug)
    {
        std::cerr << "Outputting final result ... ";
    }
    openOutputFile();  // open this file here. Do not open it way
                       // before the main writing starts.
    // it will leave a long period of zero-writing-activity (due to
    // computation), which could hang the program sometimes on panfs system

    std::ostream outputStream(&outputFilterStreamBuffer);
    outputStream << "# eGADA: enhanced Genome Alteration Detection Algorithm\n";
    outputStream << "# Authors: Yu Huang polyactis@gmail.com,"
                    " Roger Pique-Regi piquereg@usc.edu\n";
    outputStream << boost::format(
                        "# Parameters: "
                        "a=%1%,T=%2%,MinSegLen=%3%,sigma2=%4%,BaseAmp=%5%, "
                        "convergenceDelta=%6%, maxNoOfIterations=%7%, "
                        "convergenceMaxAlpha=%8%, convergenceB=%9%.") %
                        baseGADA.a % baseGADA.T % baseGADA.MinSegLen %
                        baseGADA.sigma2 % baseGADA.BaseAmp %
                        baseGADA.convergenceDelta % baseGADA.maxNoOfIterations %
                        baseGADA.convergenceMaxAlpha %
                        baseGADA.convergenceB << std::endl;
    outputStream << boost::format("# Reading M=%1% probes in input file\n") %
                        baseGADA.M;
    outputStream << boost::format("# Overall mean %1%\n") % baseGADA.Wext[0];
    outputStream << boost::format("# Sigma^2=%1%\n") % baseGADA.sigma2;
    outputStream << boost::format(
                        "# Convergence: delta=%1% after %2% EM iterations.\n") %
                        baseGADA.delta % baseGADA.numEMsteps;
    outputStream << boost::format("# Found %1% breakpoints after SBL\n") %
                        baseGADA.noOfBreakpointsAfterSBL;
    outputStream << boost::format("# Kept %1% breakpoints after BE\n") %
                        baseGADA.K;

    // Collapse Segments
    if (SelectClassifySegments == 1)
    {
        if (debug)
        {
            std::cerr << " Select and Classify Segments ...";
        }
        if (SelectEstimateBaseAmp == 1)
        {
            outputStream << boost::format("# Estimating BaseAmp\n");
            baseGADA.CompBaseAmpMedianMethod();
        }
        outputStream << boost::format("# BaseAmp=%1% \n") % baseGADA.BaseAmp;
        // outputStream<< boost::format("# Classify Segments \n", BaseAmp);

        baseGADA.CollapseAmpTtest();
        if (debug)
        {
            
            std::cerr << " SelectClassifySegments done.\n";
        }
    }

    if (SelectClassifySegments == 0)
    {
        outputStream << boost::format("Start\tStop\tLength\tAmpl\n");
        for (i = 0; i < baseGADA.K + 1; i++)
            outputStream << baseGADA.Iext[i] + 1 << "\t" << baseGADA.Iext[i + 1]
                         << "\t" << baseGADA.SegLen[i] << "\t"
                         << baseGADA.SegAmp[i] << std::endl;
    }
    else if (SelectClassifySegments == 1)
    {
        outputStream << boost::format("Start\tStop\tLength\tAmpl\tState\n");
        for (i = 0; i < baseGADA.K + 1; i++)
        {
            outputStream << baseGADA.Iext[i] + 1 << "\t" << baseGADA.Iext[i + 1]
                         << "\t" << baseGADA.SegLen[i] << "\t"
                         << baseGADA.SegAmp[i] << "\t";
            if (baseGADA.SegState[i] > baseGADA.BaseAmp)
                outputStream << "G";
            else if (baseGADA.SegState[i] < baseGADA.BaseAmp)
                outputStream << "L";
            else
                outputStream << "N";
            outputStream << endl;
        }
    }
    if (debug)
    {
        std::cerr << " output done." << endl;
    }
    closeFiles();
    if (debug)
    {
        std::cerr << "Exit eGADA.commandlineRun()." << std::endl;
    }
}
#endif  // 2009-11-21 end of the whole boost python module code

#ifdef GADABIN
int main(int argc, char *argv[])
{
    eGADA instance(argc, argv);
    instance.commandlineRun();
}
#endif
