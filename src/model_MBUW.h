#pragma once
#include "BreakPointData.h"
#include "../libs/HypothesisTester/HypothesisTester.h"
#include "ale.h"

class MBUW : public Hypothesis<BreakPointData>
{
	public:
		MBUW(std::string metaData, int nChroms, int blocksPerChrom, int singleTargetChrom) : Hypothesis<BreakPointData>((nChroms*blocksPerChrom*(nChroms*blocksPerChrom+1))/2-1,GAI)
		{
			TargetChrom = singleTargetChrom;
			std::cout << "MBUW initialised with target" << singleTargetChrom << std::endl;
			//compute genome length
			ExtractChromosomeLengths(metaData);
			
			BlocksPerChrom = blocksPerChrom; 
			int N = nChroms;
			int Nx = nChroms * blocksPerChrom;
			GridWidth = Nx;
			ComputeBlockSizes();
			std::vector<double> probLower(Nx*(Nx+1)/2-1,-10);
			std::vector<double> probUpper(Nx*(Nx+1)/2-1,10);
			SetLowerBound(probLower);
			SetUpperBound(probUpper);
			
			this->Identifier = "MBUW_" + std::to_string(blocksPerChrom);
		}
		
		double ComputeLogNormalisation(const std::vector<double> & params)
		{
			int nBlocks = BlockSizes.size();
			double logS = -99999999999999;
			int c = 0;
			for (int i = 0; i < nBlocks;++i)
			{
				for (int j = i; j < nBlocks; ++j)
				{
					double logA = log(BlockSizes[i]) + log(BlockSizes[j]);
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
			int nBlocks = BlockSizes.size();
			double L = 0;

			
			double logS = ComputeLogNormalisation(params);
			int c =0;
			for (int i = 0; i < nBlocks; ++i)
			{
				for (int j = i; j < nBlocks; ++j)
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
			double logS_opt = log(BlockSizes[anchor_y]) + log(BlockSizes[anchor_x]) + anchorVal - log(DataBin(anchor_y,anchor_x)/nData) ;

			std::vector<double> vals(this->Dimension,0);
			int c = 0;
			// int nChroms = ChromLengths.size();
			int nBlocks = BlockSizes.size();
			double log_C_anchor =  log(BlockSizes[anchor_x]) + log(BlockSizes[anchor_y]);
			for (int i = 0; i < nBlocks;++i)
			{
				for (int j = i; j < nBlocks; ++j)
				{
					bool isAnchor = (i==anchor_y && j == anchor_x);
					if (!isAnchor)
					{
						double log_amp = log(BlockSizes[i]) + log(BlockSizes[j]) - log_C_anchor;
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
			// std::cout << "Computing " <<N <<"-d Hessian" << std::endl;
			std::vector<int> pMap(N,0);
			std::vector<int> qMap(N,0);
			std::vector<double> logrs(N,0);
			// int nChrom = ChromLengths.size();
			int nBlocks = BlockSizes.size();
			double logAnchor = log(BlockSizes[anchor_x]) + log(BlockSizes[anchor_y]);
			int c = 0;
			for (int i = 0; i < nBlocks; ++i)
			{
				for (int j = i; j < nBlocks; ++j)
				{
					bool isAnchor = (i==anchor_y && j == anchor_x);
					if (!isAnchor)
					{
						pMap[c] = i;
						qMap[c] = j;
						logrs[c] = log(BlockSizes[i]) + log(BlockSizes[j]) - logAnchor;
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
			// std::cout << "Computed " <<N <<"-d Hessian" << std::endl;
			return H;
		}

		std::vector<int> ChromLengths;		
		std::vector<int> BlockSizes;
		int nData;
		int anchor_x;
		int anchor_y;
		double anchorVal = 0;
		double emptyVal = -999;
		squareMatrix DataBin;
		bool DataBinned = false;
	private:
		int BlocksPerChrom;
		int TargetChrom;
		int GridWidth;
		void ExtractChromosomeLengths(std::string metaData)
		{
			ChromLengths.resize(0);
			forLineVectorIn(metaData,' ',
				ChromLengths.push_back(std::stoi(FILE_LINE_VECTOR[1]));
			);
		}
		void ComputeBlockSizes()
		{
			int nChroms = ChromLengths.size();
			if (TargetChrom >0)
			{
				nChroms = 1;
			}
			BlockSizes.resize(GridWidth);
			int c = 0;
			for (int n = 0; n <nChroms; ++n)
			{
				int block = ChromLengths[n]/BlocksPerChrom;
				for (int j = 0; j < BlocksPerChrom; ++j)
				{
					if (j <BlocksPerChrom-1)
					{
						BlockSizes[c] = block;
					}
					else
					{
						BlockSizes[c] = ChromLengths[n] - j*block;
					}
					++c;
					// std::cout << n << "  " << block << "  " << c << "  " << GridWidth << "  " << TargetChrom << std::endl;
				}
			}
		}

		void PrecomputeData(const std::vector<BreakPointData> & Data)
		{
			// int nChroms = ChromLengths.size();
			int matSize = GridWidth;
			DataBin = squareMatrix(matSize);
			nData = Data.size();
			int largestSum = -1;
			// std::cout << "precomputing" << std::endl;
			for (int i = 0; i < nData; ++i)
			{
				// std::cout << "\t" << i << "/" << nData << std::endl;
				int x = Data[i].JoinChromosome-1;
				int y = Data[i].BreakChromosome-1;

				int chromBin_x = (double)Data[i].JoinIdx/ChromLengths[x] * BlocksPerChrom;
				int chromBin_y = (double) Data[i].BreakIdx/ChromLengths[y]*BlocksPerChrom;
				chromBin_x = std::min(std::max(0,chromBin_x),BlocksPerChrom-1);
				chromBin_y = std::min(std::max(0,chromBin_y),BlocksPerChrom-1);
				int xBin = chromBin_x;
				int yBin = chromBin_y;
				if (TargetChrom < 1)
				{
					xBin += x * BlocksPerChrom;
					yBin += y*BlocksPerChrom;
				}
				// std::cout << xBin << " " << yBin << "  "<< matSize << std::endl;
				if (xBin > yBin)
				{
					++DataBin(yBin,xBin);
				}
				else
				{
					++DataBin(xBin,yBin);
				}
				if ((xBin + yBin) > largestSum)
				{
					largestSum = xBin + yBin;
					anchor_x = std::max(xBin,yBin);
					anchor_y = std::min(xBin,yBin);
				}
			}
			
			// std::cout << DataBin.Display() << std::endl;
		}
};
