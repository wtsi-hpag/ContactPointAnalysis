#pragma once
#include <string>
#include "JSL.h"
#include <ostream>
#include "ale.h"
class Synthesiser
{
	public:
		Synthesiser(std::string syntheticOutput, int syntheticChromosomes,int synthType)
		{
			OutDir = syntheticOutput;
			Mode = synthType;
			NChroms = syntheticChromosomes;

			JSL::mkdir(OutDir);
			chromLengths.resize(NChroms);
			int minChromSize = 1e6;
			int maxChromSize = 1e7;
			GenomeSize = 0;
			std::vector<int> ns(NChroms);
			for (int i = 0; i < NChroms; ++i)
			{
				ns[i] = i+1;
				double r = (double)rand()/RAND_MAX;
				double q = log10(minChromSize) + r * log10(maxChromSize/minChromSize);
				int v = pow(10,q);
				chromLengths[i] = v;
				GenomeSize += v;
			}

			std::sort(chromLengths.begin(),chromLengths.end(),std::greater<>());
			std::string metaFile = OutDir + "/" + MetaFileName;
			JSL::initialiseFile(metaFile);
			JSL::writeMultiVectorToFile(metaFile," ",ns,chromLengths);
		}

		void Generate(int nBreaks,double param)
		{
			switch (Mode)
			{
			case 0:
				uw_generator(nBreaks);
				break;
			case 1:
				escuw_generator(nBreaks,param);
				break;
			case 2:
				scuw_generator(nBreaks,param);
				break;
			default:
				JSL::Error("This synthesis mode has not been activated yet");
				break;
			}
			
		}

	private:

		int Mode;
		int NChroms;
		unsigned long long int GenomeSize;
		std::string OutDir;
		std::string MetaFileName = "chromLengths.dat";
		std::string BreakFileName = "breakpoint.dat";
		std::vector<int> chromLengths;

		void uw_generator(int nBreaks)
		{
			std::string breakFile = OutDir + "/" + BreakFileName;
			JSL::initialiseFile(breakFile);
			std::ostringstream os;
			for (int i = 0; i < nBreaks; ++i)
			{
				double r1 = (double)rand()/RAND_MAX;
				double r2 = (double)rand()/RAND_MAX;
				double r3 = (double)rand()/RAND_MAX;
				unsigned long long int startIdx = r1 *GenomeSize;
				unsigned long long int endIdx = r2 *GenomeSize;
				int coverage = 2 + 48*r3;
				
				int startChrom = -1;
				int endChrom = -1;
				unsigned long long int cumulativeIdx =0;
				for (int i = 0; i < chromLengths.size(); ++i)
				{
					cumulativeIdx += chromLengths[i];
					if (startChrom == -1 && startIdx< cumulativeIdx)
					{
						startChrom = i+1;
						startIdx -= cumulativeIdx - chromLengths[i];
					}
					if (endChrom == -1 && endIdx< cumulativeIdx)
					{
						endChrom = i+1;
						endIdx -= cumulativeIdx - chromLengths[i];
					}
					
					if (startChrom > -1 && endChrom >-1)
					{
						break;
					}
				}
				spoofBreakLine(os,i,startChrom,startIdx,endChrom,endIdx,coverage);
			}
			JSL::writeStringToFile(breakFile,os.str());
		}

		void escuw_generator(int nBreaks,double param)
		{
			double sigma = param;
			double pBreak = exp(sigma);
			std::string breakFile = OutDir + "/" + BreakFileName;
			JSL::initialiseFile(breakFile);
			std::ostringstream os;

			std::vector<double> chromWeigtings(NChroms,0.0);
			chromWeigtings[0] = (double)chromLengths[0]/GenomeSize;
			for (int i = 1; i < NChroms; ++i)
			{
				chromWeigtings[i] = chromWeigtings[i-1] + (double)chromLengths[i]/GenomeSize;
			}

			for (int i = 0; i < nBreaks; ++i)
			{
				double r1 = (double)rand()/RAND_MAX;
				double r2 = (double)rand()/RAND_MAX;
				double r3 = (double)rand()/RAND_MAX;
				int startChrom;
				for (int q = 0; q < NChroms; ++q)
				{
					if (r1 < chromWeigtings[q])
					{
						startChrom = q;
						break;
					}
				}

				int endChrom = startChrom;
				int coverage = 2 + 48*r3;
				if (r2 < pBreak)
				{
					while (endChrom == startChrom)
					{
						r3 = (double)rand()/RAND_MAX;
						for (int q = 0; q < NChroms; ++q)
						{
							if (r3 < chromWeigtings[q])
							{
								endChrom = q;
								break;
							}
						}
					}
				}

				int startIdx = (double)rand()/RAND_MAX * chromLengths[startChrom];
				int endIdx = (double)rand()/RAND_MAX * chromLengths[endChrom];
				spoofBreakLine(os,i,startChrom+1,startIdx,endChrom+1,endIdx,coverage);
			}
			JSL::writeStringToFile(breakFile,os.str());
		}

		void scuw_generator(int nBreaks, double diagonalAmpify)
		{
			int nParam = (NChroms * (NChroms + 1))/2;
			std::vector<double> lambdas(nParam,0);
			//Generate ramdom params
			for (auto & lambda : lambdas)
			{
				lambda = (double)rand()/RAND_MAX*-10;
			}


			//compute the normalisation
			double S = -9999;
			int c = 0;
			for (int i = 0; i < NChroms; ++i)
			{
				for (int j = i; j < NChroms; ++j)
				{
					double logC = log(chromLengths[i]) + log(chromLengths[j]) + lambdas[c];
					if (i!=j)
					{
						logC += log(2);
					}
					else
					{
						logC += diagonalAmpify; //external parameter makes diagonal elements more preferable
					}
					++c;
					S = ale(S,logC);
				}
			}

			//populate the matrix
			std::vector<std::vector<double>>  ChromDistr(NChroms,std::vector<double>(NChroms,0));
			c = 0;
			double s = 0;
			for (int i = 0; i < NChroms;++i)
			{
				for (int j = i; j < NChroms; ++j)
				{
					double logA = log(chromLengths[i]) + log(chromLengths[j]);
					if (i == j)
					{
						logA += diagonalAmpify;
					}
					
					ChromDistr[i][j] = exp(lambdas[c] - S + logA);
					ChromDistr[j][i] = ChromDistr[i][j]; //symmetry!
					s += ChromDistr[i][j];
					if (i!=j)
					{
						s+=ChromDistr[i][j];
					}
					
					++c;
				}
			}

			//generate cumulative probability arrays
			std::vector<double> pFirstChoice(NChroms,0);
			std::vector<std::vector<double>> pSecondChoice(NChroms,std::vector<double>(NChroms,0));
			double prev1 = 0;
			for (int i = 0; i < NChroms; ++i)
			{
				double s = 0;
				for (int j = 0; j < NChroms; ++j)
				{
					s += ChromDistr[i][j];
				}
				pFirstChoice[i] = prev1 + s;
				prev1 = pFirstChoice[i];
				double prev2 = 0;
				double q = 0;
				for (int j = 0; j < NChroms; ++j)
				{
					pSecondChoice[i][j] = prev2 + ChromDistr[j][i]/s;
					prev2 = pSecondChoice[i][j];
				}
			}


			std::string breakFile = OutDir + "/" + BreakFileName;
			JSL::initialiseFile(breakFile);
			std::ostringstream os;
			for (int b = 0; b < nBreaks; ++b)
			{
				double rBreak = (double)rand()/RAND_MAX;
				double rJoin = (double)rand()/RAND_MAX;
				int coverage = 2 + 48*(double)rand()/RAND_MAX;
				int bChrom;
				int jChrom;
				for (int i = 0; i < NChroms; ++i)
				{
					if (rBreak <= pFirstChoice[i])
					{
						bChrom = i;
						break;
					}
				}
				for (int i = 0; i < NChroms; ++i)
				{
					if (rJoin <= pSecondChoice[bChrom][i])
					{
						jChrom = i;
						break;
					}
				}
				int startIdx = (double)rand()/RAND_MAX * chromLengths[bChrom];
				int endIdx = (double)rand()/RAND_MAX * chromLengths[jChrom];
				spoofBreakLine(os,b,bChrom+1,startIdx,jChrom+1,endIdx,coverage);
			}
			JSL::writeStringToFile(breakFile,os.str());
		}
		void spoofBreakLine(std::ostringstream & os,int id, int chrom1, int idx1, int chrom2, int idx2, int coverage)
		{
			os << id << " " << chrom1 << " " << idx1 << " " << chrom2 << " " << idx2 << " " << coverage << " || 0 0 synthetic_read\n";
		}
};