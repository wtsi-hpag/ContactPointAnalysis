#pragma once
#include "BreakPointData.h"
#include "../libs/HypothesisTester/HypothesisTester.h"
#include "ale.h"

class SCUW : public Hypothesis<BreakPointData>
{
	public:
		SCUW(std::string metaData, int nChroms) : Hypothesis<BreakPointData>((nChroms*(nChroms+1))/2-1,GAI)
		{
			//compute genome length
			ExtractChromosomeLengths(metaData);

			int N = nChroms;
			std::vector<double> probLower(nChroms*(nChroms+1)/2-1,-10);
			std::vector<double> probUpper(nChroms*(nChroms+1)/2-1,10);
			SetLowerBound(probLower);
			SetUpperBound(probUpper);

			this->Identifier = "SCUW";
		}
		squareMatrix ConstructMatrix(const std::vector<double> & params)
		{
			int nChroms = ChromLengths.size();
			squareMatrix display(nChroms);
			double logS = ComputeLogNormalisation(params);
			int c = 0;
			double ss =0 ;
			for (int i = 0; i < nChroms;++i)
			{
				for (int j = i; j < nChroms; ++j)
				{
					double logA = log(ChromLengths[i]) + log(ChromLengths[j]);
					double v;
					if ( i == anchor_y && j == anchor_x)
					{
						v = anchorVal;
					}
					else
					{
						v = params[c];
						++c;
					}
					double V = exp(v + logA - logS);

					display(i,j) = V;
					display(j,i) = V;
					ss += exp(v + logA - logS);
					if (i!=j)
					{
						ss += exp(v + logA - logS);
					}
				}
			}
			return display;
		}
		double ComputeLogNormalisation(const std::vector<double> & params)
		{
			int nChroms = ChromLengths.size();
			double logS = -99999999999999;
			int c = 0;
			for (int i = 0; i < nChroms;++i)
			{
				for (int j = i; j < nChroms; ++j)
				{
					double logA = log(ChromLengths[i]) + log(ChromLengths[j]);
					if (i != j)
					{
						logA += log(2);
					}					

					bool isAnchor = (i==anchor_y) && (j == anchor_x);
					double v;
					if (isAnchor)
					{
						v = logA + anchorVal;
					}
					else
					{
						v =  logA + params[c];
						++c;
					}

					logS = ale(logS,v);
				}
			}
			return logS;
		}
		double LogProbability(const std::vector<BreakPointData> & Data, const std::vector<double> & params)
		{
			int N =sqrt(this->Dimension);
			int nChroms = ChromLengths.size();
			double L = 0;

			
			double logS = ComputeLogNormalisation(params);
			int c =0;
			for (int i = 0; i < nChroms; ++i)
			{
				for (int j = i; j < nChroms; ++j)
				{
					bool isAnchor = (i==anchor_y && j ==anchor_x);
					double lambda = anchorVal;
					if (!isAnchor)
					{
						lambda = params[c];
						++c;
					}
					else
					{
						lambda = anchorVal;
					}
					L += lambda * DataBin(i,j);
				}
			}
			L -= nData * logS;
			return L;
		}

		std::vector<double> FindMaximum(const std::vector<BreakPointData> & Data, std::vector<double> & initialGuess, int N)
		{
			PrecomputeData(Data);
			double logS_opt = log(ChromLengths[anchor_y]) + log(ChromLengths[anchor_x]) + anchorVal - log(DataBin(anchor_y,anchor_x)/nData) ;

			std::vector<double> vals(this->Dimension,0);
			int c = 0;
			int nChroms = ChromLengths.size();
			double log_C_anchor =  log(ChromLengths[anchor_x]) + log(ChromLengths[anchor_y]);
			for (int i = 0; i < nChroms;++i)
			{
				for (int j = i; j < nChroms; ++j)
				{
					bool isAnchor = (i==anchor_y && j == anchor_x);
					if (!isAnchor)
					{
						double log_amp = log(ChromLengths[i]) + log(ChromLengths[j]) - log_C_anchor;
						if (i!=j)
						{
							log_amp += log(2);
						}
						if (DataBin(i,j) > 0)
						{
							vals[c] = log(DataBin(i,j)/DataBin(anchor_y,anchor_x)) - log_amp + anchorVal;
						}
						else
						{
							vals[c] = emptyVal;
						}
						++c;
					}
				}
			}
			return vals;
		}

		squareMatrix ProbabilityHessian(const std::vector<BreakPointData> & Data, std::vector<double> & params)
		{
			int N = this->Dimension;
			std::vector<int> pMap(N,0);
			std::vector<int> qMap(N,0);
			std::vector<double> logrs(N,0);
			int nChrom = ChromLengths.size();
			double logAnchor = log(ChromLengths[anchor_x]) + log(ChromLengths[anchor_y]);
			int c = 0;
			for (int i = 0; i < nChrom; ++i)
			{
				for (int j = i; j < nChrom; ++j)
				{
					bool isAnchor = (i==anchor_y && j == anchor_x);
					if (!isAnchor)
					{
						pMap[c] = i;
						qMap[c] = j;
						logrs[c] = log(ChromLengths[i]) + log(ChromLengths[j]) - logAnchor;
						if (i!=j)
						{
							logrs[c] += log(2);
						}
						++c;
					}
				}
			}
			squareMatrix H(N);
			double log_nAnchor = log(DataBin(anchor_y,anchor_x));
			for (int q = 0; q < N; ++q)
			{
				for (int p = q; p < N; ++p)
				{
					double t1 = std::max(1.0/nData,DataBin(pMap[q],qMap[q]));
					double t2 = std::max(1e-5,(double)DataBin(pMap[p],qMap[p])/nData);
					if (p == q)
					{
						t2 -=1;
					}
					double dL = t1*t2;
					H(p,q) = dL;
					H(q,p) = dL;

				}
			}
			std::cout << "\t" <<N<<"-dimensional Hessian computed" << std::endl;
			return H;
		}

		std::vector<int> ChromLengths;		
		int nData;
		int anchor_x;
		int anchor_y;
		double anchorVal = 0;
		double emptyVal = -999;
		squareMatrix DataBin;
		bool DataBinned = false;
	private:
		void ExtractChromosomeLengths(std::string metaData)
		{
			ChromLengths.resize(0);
			forLineVectorIn(metaData,' ',
				ChromLengths.push_back(std::stoi(FILE_LINE_VECTOR[1]));
			);
		}

		void PrecomputeData(const std::vector<BreakPointData> & Data)
		{
			int nChroms = ChromLengths.size();
			DataBin = squareMatrix(nChroms);
			nData = Data.size();
			int largestSum = -1;
			for (int i = 0; i < nData; ++i)
			{
				int x = Data[i].JoinChromosome-1;
				int y = Data[i].BreakChromosome-1;
				if (x > y)
				{
					++DataBin(y,x);
				}
				else
				{
					++DataBin(x,y);
				}
				if ((x + y) > largestSum)
				{
					largestSum = x + y;
					anchor_x = std::max(x,y);
					anchor_y = std::min(x,y);
				}
			}
		}
};
