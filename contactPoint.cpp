// #define GNUPLOT_NO_TIDY
#include "libs/HypothesisTester/HypothesisTester.h"

#include "JSL.h"
#include "src/model_UW.h"
#include "src/model_ESCUW.h"
#include "src/model_SCUW.h"
#include "src/synthesiser.h"
std::vector<BreakPointData> loadBreaks(std::string file,int coverageThreshold)
{
	std::vector<BreakPointData> out;

	forLineVectorIn(file,' ',
		BreakPointData breaker(FILE_LINE_VECTOR);

		if (breaker.Coverage > coverageThreshold)
		{
			out.push_back(breaker);
		}
	)
	return out;
}

void DetectionMode(std::string breakFile, std::string metaFile, int coverage)
{
	int Nchroms = JSL::LineCount(metaFile)-1;

	ModelTester<BreakPointData> mt;
	mt.AddHypothesis(UW(metaFile));
	mt.AddHypothesis(ESCUW(metaFile));
	mt.AddHypothesis(SCUW(metaFile,Nchroms));

	auto breakData = loadBreaks(breakFile,coverage);
	auto results =	mt.BeginTest(breakData,1000);

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
	JSL::Argument<double> sigmaMove(-999,"sigma",argc,argv);
	JSL::Argument<int> CoverageThreshold(2,"coverage",argc,argv);
	if (synthesisMode)
	{
		std::cout << "SYNTHESISING!" << std::endl;
		SynthesisMode(syntheticOutput.Value,nSynthChromosomes,synthType,nSynth,sigmaMove);
	}
	else
	{
		DetectionMode(breakFile,metaFile,CoverageThreshold);
		// GridMode(breakFile,metaFile);
	}
	return 0;
}