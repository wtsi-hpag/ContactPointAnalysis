// #define GNUPLOT_NO_TIDY
#include "libs/HypothesisTester/HypothesisTester.h"

#include "libs/JSL/JSL.h"
#include "src/model_UW.h"
#include "src/model_ESCUW.h"
#include "src/model_SCUW.h"
#include "src/model_MBUW.h"
#include "src/model_SAW.h"
#include "src/synthesiser.h"

std::vector<BreakPointData> loadBreaks(std::string file,int coverageThreshold, int scm)
{
	std::vector<BreakPointData> out;
	int lowestSep = -1;
	forLineVectorIn(file,' ',
		BreakPointData breaker(FILE_LINE_VECTOR);

		if (scm <= 0 || (breaker.JoinChromosome == scm && breaker.BreakChromosome == scm))
		{
			if (breaker.Coverage > coverageThreshold)
			{
				out.push_back(breaker);
				bool same = breaker.JoinChromosome == breaker.BreakChromosome;
				long long int sep = abs(breaker.JoinIdx - breaker.BreakIdx); 
				if (lowestSep == -1 || (same && sep < lowestSep))
				{
					lowestSep = sep;
					// std::cout << "New smallest sep " << breaker.JoinIdx << "  " << breaker.BreakIdx << std::endl;
				}
			}
			// std::cout<< breaker.JoinChromosome <<"  " << breaker.BreakChromosome << std::endl;
		}
	)
	// std::cout <<"The smallest separation I found was " << lowestSep << std::endl;
	return out;
}

void DetectionMode(std::string breakFile, std::string metaFile, int coverage, int singleChromosomeMode, int threads)
{
	int Nchroms = JSL::LineCount(metaFile)-1;
	ModelTester<BreakPointData> mt(threads);
	mt.Verbosity = 2;
	if (singleChromosomeMode > 0)
	{
		mt.AddHypothesis(UW(metaFile, singleChromosomeMode));
		for (int i = 2; i < 10; i+= 1)
		{
			mt.AddHypothesis(MBUW(metaFile, 1, i,singleChromosomeMode));
		}
		std::string data = "Data/HiC_Test/";
		// mt.AddHypothesis(SAW(1e3,metaFile,data,singleChromosomeMode));
		// mt.AddHypothesis(SAW(1e5,metaFile,data,singleChromosomeMode));
		// mt.AddHypothesis(SAW(1e8,metaFile,data,singleChromosomeMode));
		// mt.AddHypothesis(SAW(1e10,metaFile,data,singleChromosomeMode));
		// mt.AddHypothesis(SAW(1e30,metaFile,data,singleChromosomeMode));
	}
	else
	{
		mt.AddHypothesis(UW(metaFile));
		mt.AddHypothesis(ESCUW(metaFile));
		mt.AddHypothesis(SCUW(metaFile,Nchroms));
		// mt.AddHypothesis(MBUW(metaFile, Nchroms, 1,singleChromosomeMode));
		mt.AddHypothesis(MBUW(metaFile, Nchroms, 2,singleChromosomeMode));
		mt.AddHypothesis(MBUW(metaFile, Nchroms, 3,singleChromosomeMode));
		mt.AddHypothesis(MBUW(metaFile, Nchroms, 4,singleChromosomeMode));
		mt.AddHypothesis(MBUW(metaFile, Nchroms, 5,singleChromosomeMode));
	}
	auto breakData = loadBreaks(breakFile,coverage, singleChromosomeMode);
	auto results =	mt.BeginTest(breakData,100000);

	std::cout << "Best Fitting Model is: " << results.BestModel << std::endl;

	auto score = results.Scores;
	
	auto maxScore = *std::min_element(score.begin(),score.end());
	std::vector<double> scoreDiffs(score.size());
	double min = -1;
	for (int i = 0; i < score.size(); ++i)
	{
		double logDiff = score[i];
		scoreDiffs[i] = abs(logDiff/log(10)-1);
		if (scoreDiffs[i] > min)
		{
			min = abs(scoreDiffs[i]);
		}
	}

	JSL::gnuplot gp;
	gp.Chart(results.Models,scoreDiffs);
	gp.SetYLog(true);
	gp.SetYRange(0.2,min);
	gp.SetXLabel("Model");
	gp.SetYLabel("\\log_{10} Odds-Ratio Relative to Worst Model");
	gp.SetGrid(true);
	gp.Show();

}

void SynthesisMode(std::string syntheticOutput, int syntheticChromosomes, int synthType, int N,double param)
{
	Synthesiser s(syntheticOutput,syntheticChromosomes,synthType);
	s.Generate(N,param);
}

void GridMode(std::string breakFile, std::string metaFile)
{
	int binWidth = 3;
	int inf = 120;
	std::vector<std::vector<int>> reads;
	int trailingEdge = 0;
	std::vector<double> edge;
	double c = 0;
	while (trailingEdge < inf)
	{	
	
		trailingEdge += binWidth;
		edge.push_back(trailingEdge);
	}
	edge.push_back(130);

	int nCats = 4;
	reads.resize(nCats);
	for (int i = 0; i < nCats; ++i)
	{
		reads[i].resize(edge.size(),0);
	}
	std::vector<int> sum(nCats,0);
	forLineVectorIn(breakFile,' ',
		int breakChrom = std::stoi(FILE_LINE_VECTOR[1]);
		int joinChrom = std::stoi(FILE_LINE_VECTOR[3]);
		int coverage = std::stoi(FILE_LINE_VECTOR[5]);

		std::vector<int> cats = {0};
		// std::cout << joinChrom << "  " << breakChrom << std::endl;
		if (joinChrom == 6 || breakChrom == 6)
		{
			cats.push_back(1);
		}
		else
		{
			if (joinChrom == 9 || breakChrom == 9)
			{
				cats.push_back(2);
			}
			else
			{
				cats.push_back(3);
			}
		}

		if (cats.size() > 0)
		{
			for(int j = 0; j < edge.size(); ++j)
			{
				if (coverage < edge[j])
				{
					for (int cat : cats)
					{
						++reads[cat][j];
						++sum[cat];
					}
					break;
				}
			}
		}
	);

	JSL::gnuplot gp;
	namespace lp = JSL::LineProperties;
	gp.WindowSize(900,900);
	// gp.SetMultiplot(nCats,1);
	std::vector<std::string> titles = {"All","Chromosome 6","Chromosome 9","Remainder"};
	for (int i = 0; i < nCats; ++i)
	{
		// gp.SetAxis(i);
		std::vector<int> x;
		std::vector<double> y;
		for (int j = 0; j < edge.size(); ++j)
		{
			int end = edge[j];
			int start = 0;
			if (j > 0)
			{
				start = edge[j-1];
			}
			if (j == edge.size() -1)
			{
				end = start + binWidth;
			}
			// x.push_back(start);
			x.push_back(start);
			x.push_back(end);
			// x.push_back(end);
			// y.push_back(0);
			double v = (double)reads[i][j]/sum[i];
			y.push_back(v);
			y.push_back(v);
			// y.push_back(0);
		}
		gp.Plot(x,y,lp::Legend(titles[i]),lp::PenSize(3));
		// gp.SetTitle(titles[i]);
		gp.SetXLabel("Reads Confirming Breaks");
		gp.SetYLabel("Occurence");
		// gp.SetXLog(true);
		gp.SetXRange(2,40);
	}
	gp.SetLegend(true);
	gp.Show();
	
}

inline double erfDist(double x,double x0, double ell)
{
	return (x - x0)/(sqrt(2)*ell);
}

void PrepareHiC()
{
	std::string hiCFile = "Data/chr6.sort-117";
	std::vector<long long int> ChromLengths;
	std::string chFile = "Data/chromLengths.dat";
	forLineVectorIn(chFile,' ',
		int length = std::stoi(FILE_LINE_VECTOR[1]);
		ChromLengths.push_back(length);
		// GenomeLength += length;
	);
	int nbins = 10000;

	std::vector<int> x = JSL::Vector::linspace(0,ChromLengths[5],nbins);
	auto y = x;
	double DX = x[1] - x[0];
	std::cout<< "Bin size = " << (x[1] - x[0])/1e6 << "Mb"<< std::endl;
	std::vector<std::vector<double>> grid(nbins,std::vector<double>(nbins,0.0));

	double ell = 2e4;//DX*2;
	// int lM = 1000000;
	int l = 0;
	double pre = log(1.0/sqrt(2*M_PI*ell*ell));
	int expected = JSL::LineCount(hiCFile);
	JSL::ProgressBar pb(expected);

	int influencedSquares = std::max(ell/(DX)*5,1.0);
	std::cout << "Correlation Distance: " << ell << std::endl;
	std::cout << "Influencing " << influencedSquares << " squares " << std::endl;

	forLineVectorIn(hiCFile,' ',
	// std::cout << l << std::endl;
		int i = std::stoi(FILE_LINE_VECTOR[1]);
		int j = std::stoi(FILE_LINE_VECTOR[2]);

		if (abs(i - j) > 3.5e5)
		{
			int iBin = (double)i/ChromLengths[5] * nbins;
			int jBin = (double)j/ChromLengths[5] * nbins;
			// ++grid[iBin][jBin];
			// std::cout << "Join at " << i << "  " << j << " assigned bins " << iBin << "  " << jBin << std::endl;
			for (int q = std::max(0,iBin-influencedSquares); q < std::min(nbins,iBin + influencedSquares); ++q)
			{
				for (int p = std::max(0,jBin-influencedSquares); p< std::min(nbins,jBin + influencedSquares); ++p)
				{
					double dx = ((double)(x[p] - j))/ell;
					double dy = ((double)(y[q]- i))/ell;
					double rSq= dx*dx + dy*dy;
					// double accept = std::max(DX*2/ell,5.0);
					// std::cout << dx << "  " << dy << std::endl;
					// 
					// if (rSq < accept*accep÷
						// double v2 =  0.5*exp(2*pre -0.5*(dx*dx + dy*dy));
						
						double v = 1.0/(8*DX*DX)* (erf( erfDist(x[p]+DX,j,ell)) - erf(erfDist(x[p],j,ell))) * (erf( erfDist(x[q]+DX,i,ell)) - erf(erfDist(x[q],i,ell)));
						grid[q][p] += v;
						grid[p][q] += v;
						// std::cout << "\tAt pos " << x[p] << "  " << x[q] << " I estimate the value as " << DX*DX*v << std::endl;
						// std::cout << "Adding to " << p << "  " << q << "  " << i << "  " << j << std::endl;
					// }
				}
			}
		}
		pb.Update(l);
		++l;
		// if (l>=lM)
		// {
		// 	break;
		// }
	);

	double S = 0;
	for (int i = 0; i < nbins; ++i)
	{
		for (int j = 0; j < nbins; ++j)
		{
			S += DX*DX*grid[i][j];
		}
	}
	for (int i = 0; i < nbins; ++i)
	{
		for (int j = 0; j < nbins; ++j)
		{
			grid[i][j]/=S;
		}
	}
	std::cout << "Grid Norm: " << S << "  expected " << l << "  err = " << (S-l)/l << std::endl;
	JSL::initialiseFile("hic_map.dat");
	JSL::writeMatrixToFile("hic_map.dat",grid," ","\n");
}


int main(int argc, char ** argv)
{
	JSL::Argument<int> seed(time(NULL),"s",argc,argv);
	srand(seed);
	rand();
	JSL::Argument<std::string> breakFile("Data/breakpoint.dat","break",argc,argv);
	JSL::Argument<std::string> metaFile("Data/chromLengths.dat","meta",argc,argv);
	JSL::Toggle synthesisMode(false,"synth",argc,argv);
	JSL::Argument<std::string> syntheticOutput("Data/Synthetic/","o",argc,argv);
	JSL::Argument<int> nSynthChromosomes(5,"n",argc,argv);
	JSL::Argument<int> synthType(0,"t",argc,argv);
	JSL::Argument<int> nSynth(1000,"ng",argc,argv);
	JSL::Argument<int> SingleChromosomeMode(-1,"scm",argc,argv);
	JSL::Argument<double> sigmaMove(-999,"sigma",argc,argv);
	JSL::Argument<int> CoverageThreshold(2,"coverage",argc,argv);
	JSL::Argument<int> threads(1,"thread",argc,argv);

	loadBreaks(breakFile,10,SingleChromosomeMode);
	if (synthesisMode)
	{
		std::cout << "SYNTHESISING!" << std::endl;
		SynthesisMode(syntheticOutput.Value,nSynthChromosomes,synthType,nSynth,sigmaMove);
	}
	else
	{
		DetectionMode(breakFile,metaFile,CoverageThreshold,SingleChromosomeMode,threads);
		// GridMode(breakFile,metaFile);
		// PrepareHiC();
	}
	return 0;
}