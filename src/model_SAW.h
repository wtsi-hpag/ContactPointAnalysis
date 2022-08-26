#pragma once
#include "../libs/HypothesisTester/HypothesisTester.h"
#include "BreakPointData.h"
#include "../libs/JSL/JSL.h"
#include <filesystem>
namespace fs = std::filesystem;
class SAW : public Hypothesis<BreakPointData>
{
	public:
		SAW(double correlation, std::string metaData,std::string hicDirectory) : Hypothesis<BreakPointData>(0)
		{
			//compute genome length
			TargetChrom = -1;
			JSL::Error("Multi-chromosomal SAW is not yet supported");
			ExtractChromosomeLengths(metaData);
			this->Identifier = "SAW Model";
			GaussianLength = correlation;
			Directory = hicDirectory;
			GetFiles();
		}
		SAW(double correlation, std::string metaData, std::string hicDirectory, int targetChrom) : Hypothesis<BreakPointData>(0)
		{
			//compute genome length
			TargetChrom = targetChrom;
			JSL::Assert("Multi-Chromosomal SAW is not yet supported",TargetChrom > 0);
			ExtractChromosomeLengths(metaData);
			double rounded = round(log10(correlation)*10)/10;
			std::string id = std::to_string((int)rounded) + "." + std::to_string((int)(10*(rounded - (int)rounded)));
			this->Identifier = "SAW_{" + id +"}";
			GaussianLength = correlation;
			Directory = hicDirectory;
			GetFiles();
		}
		double LogProbability(const std::vector<BreakPointData> & Data, const std::vector<double> & params)
		{
			double lim = 5;
			
			int N = JSL::LineCount(files[0]);



			double sumInt = 1;
			for (double i = 1; i < ChromLengths[TargetChrom-1]; ++i)
			{
				// for (double j = 1; j < ChromLengths[TargetChrom-1]; ++1)
				// {
					double r = (double)i * (double)i;
					sumInt += 2*exp( -r/(2*GaussianLength*GaussianLength));
				// }
			}


			double prefactor = -log(2*M_PI*sumInt*sumInt) - log(N);
			int l = 0;
			std::vector<double> scores(Data.size(),-9999);
			std::vector<int> associations(Data.size(),0);
			std::vector<double> closeMisses(Data.size(),-1);
			std::vector<double> bigMisses(Data.size(),-1);
			double M = 0;
			double Q = 0;
			long long int nQ = 0;
			forLineVectorIn(files[0],' ',
				int ref_x = std::stoi(FILE_LINE_VECTOR[1]);
				int ref_y = std::stoi(FILE_LINE_VECTOR[2]);

				// int ref_x = std::max(x,y);
				// int ref_y = std::min(x,y);
				double up_x = (double)((ChromLengths[TargetChrom-1] - ref_x))/(sqrt(2)*GaussianLength);
				double down_x = (double)(( ref_x))/(sqrt(2)*GaussianLength);
				double up_y = (double)((ChromLengths[TargetChrom-1] - ref_y))/(sqrt(2)*GaussianLength);
				double down_y = (double)(( ref_y))/(sqrt(2)*GaussianLength);
				double normx =  sqrt(M_PI/2) * (erf(up_x) + erf(down_x)) * GaussianLength;
				double normy = sqrt(M_PI/2) * (erf(up_y) + erf(down_y)) * GaussianLength;
				
				double N6 = ChromLengths[TargetChrom-1];
				double accuratePrefactor = -log(normx*normy) - log(N);
				M+=accuratePrefactor;
				

				for (int i = 0; i < Data.size(); ++i)
				{
					double rSq = pow(Data[i].JoinIdx - ref_x,2) + pow(Data[i].BreakIdx-ref_y,2);
					double gaussArg = rSq/(GaussianLength*GaussianLength);
					if (gaussArg < lim*lim || scores[i] < -1000)
					{


						double contribution = -0.5*gaussArg + accuratePrefactor;
						Q += contribution;
						++nQ;
						if (scores[i] == -9999)
						{
							scores[i] = contribution;
						}
						else
						{
							scores[i] = ale(scores[i],contribution);
						}
						associations[i]++;
					}

					double r = sqrt(gaussArg);
					if (closeMisses[i] == -1 || r < closeMisses[i])
					{
						closeMisses[i] = r;
					}
					if (bigMisses[i] == -1 || r > bigMisses[i])
					{
						bigMisses[i] = r;
					}
				}
				++l;
				
			);
			double v = 0;
			for (int i = 0; i < scores.size(); ++i)
			{
				v+=scores[i];
			}
			return v;
		}
		
		int TargetChrom;
	private:
		int HiCResolution;
		double GaussianLength;

		std::vector<int> ChromLengths;
		std::vector<std::string> files;
		std::string Directory;
		void ExtractChromosomeLengths(std::string metaData)
		{
			ChromLengths.resize(0);
			forLineVectorIn(metaData,' ',
				ChromLengths.push_back(std::stoi(FILE_LINE_VECTOR[1]));
			);
		}

		void GetFiles()
		{
			for (const auto & entry : fs::directory_iterator(Directory))
			{
		        auto fileArray = JSL::split(entry.path(),'/');
				std::string fileName = fileArray[fileArray.size()-1];
				fileName = JSL::split(fileName,'.')[0];
				int chrom = std::stoi(fileName.substr(3,fileName.size()-3));
				
				//do something about checking chrom vs allowed ones
				if (TargetChrom> 1 && chrom == TargetChrom)
				{
					files.push_back(entry.path());
				}
			}
		}
};
