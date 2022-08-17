#pragma once
#include "BreakPointData.h"
#include "../libs/HypothesisTester/HypothesisTester.h"

double log_factorial(unsigned int n)
{
	if (n > 10)
	{//sterling approx
		return 0.5*log(2*M_PI) + (n+0.5) * log(n) - n;
	}
	else
	{
		double v = 0;
		for (int i = 1; i < n; ++i)
		{
			v += log(n);
		}
		return v;
	}
}

class ESCUW : public Hypothesis<BreakPointData>
{
	public:
		ESCUW(std::string metaData) : Hypothesis<BreakPointData>(1)
		{
			//compute genome length
			SetLowerBound({-10});
			SetUpperBound({-0.02});
			ExtractChromosomeLengths(metaData);


			this->Identifier = "ESCUW";
		}
		double LogProbability(const std::vector<BreakPointData> & Data, const std::vector<double> & params)
		{
			double logpMove = params[0];
			double logpStick = log(1.0 - exp(logpMove));
			// return nMove * logpMove + nStick * log(1.0 - exp(logpMove)) + LogChromSum;
			
			double lp = nStick*logpStick + nMove * logpMove + LogChromSum;
			return lp;
		}

		void LogProbabilityGradient(std::vector<double> & g, const std::vector<BreakPointData> & Data, std::vector<double> & params)
		{
			// std::cout<< "log prob" << std::endl;
			double logpMove = params[0];
			double pMove = exp(logpMove);
			double pStick = 1.0 - pMove;
			g[0] = nMove - nStick * pMove/(pStick);
		}

		squareMatrix ProbabilityHessian( const std::vector<BreakPointData> & Data, std::vector<double> & params)
		{
			squareMatrix H(params.size());
			double logpMove = params[0];
			double pMove = exp(logpMove);
			double pStick = 1.0 - pMove;
			double term = pMove/pow(pStick,2);
			H(0,0) = -nStick * term;
			// std::cout << "\tpMove " << pMove << "   " << H.Display() << "  DET =  " << H.log_LU_Determinant() << std::endl;
			return H;
		}

		double Score(const std::vector<BreakPointData> Data,int resolution)
		{
			return Score(Data);
		}
		double Score(const std::vector<BreakPointData> & Data)
		{
			PrecomputeData(Data);
	
				double t1 = log_factorial(nMove); 
				double t2 = log_factorial(nStick -1);
				double t3 = log_factorial(nMove + nStick);
				return t1 + t2 -t3  + LogChromSum;
	
		}
		// std::vector
		std::vector<unsigned int> ChromLengths;
		unsigned long long int GenomeLength;
	private:
		void ExtractChromosomeLengths(std::string metaData)
		{
			GenomeLength = 0;
			ChromLengths.resize(0);
			forLineVectorIn(metaData,' ',
				int length = std::stoi(FILE_LINE_VECTOR[1]);
				ChromLengths.push_back(length);
				GenomeLength += length;
			);
		}


		int nMove;
		int nStick;
		
		double LogChromSum;
		void PrecomputeData(const std::vector<BreakPointData> & Data)
		{
			nMove = 0;
			nStick = 0;
			
			int N = Data.size();
			LogChromSum = 0;
			for (int i = 0; i < N; ++i)
			{

				int join = Data[i].JoinChromosome-1;
				int breaker = Data[i].BreakChromosome-1;

				if (join == breaker)
				{
					++nStick;
					LogChromSum -= log(GenomeLength) + log(ChromLengths[join]);
				}
				else
				{
					++nMove;
					LogChromSum += -log(GenomeLength) -log(2) + log(1.0/(GenomeLength - ChromLengths[join]) + 1.0/(GenomeLength - ChromLengths[breaker]));
				}

				
			}

			std::cout << "\tPrecomputing, found " << nMove << "  " << nStick << "  " << LogChromSum << std::endl;
		}
};
