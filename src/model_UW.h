#pragma once
#include "../libs/HypothesisTester/HypothesisTester.h"
#include "BreakPointData.h"

class UW : public Hypothesis<BreakPointData>
{
	public:
		UW(std::string metaData) : Hypothesis<BreakPointData>(0)
		{
			//compute genome length
			TargetChrom = -1;
			ExtractGenomeLength(metaData);
			this->Identifier = "UW Model";
		}
		UW(std::string metaData, int targetChrom) : Hypothesis<BreakPointData>(0)
		{
			//compute genome length
			TargetChrom = targetChrom;
			ExtractGenomeLength(metaData);
			this->Identifier = "UW Model";
		}
		double LogProbability(const std::vector<BreakPointData> & Data, const std::vector<double> & params)
		{
			return -2.0 * Data.size() * log((double)GenomeLength);
		}
		unsigned long long int GenomeLength;
		int TargetChrom;
	private:
		void ExtractGenomeLength(std::string metaData)
		{
			GenomeLength = 0;
			int i = 0;
			forLineVectorIn(metaData,' ',
				if (TargetChrom <= 0 || i == TargetChrom-1)
				{
					GenomeLength += std::stoi(FILE_LINE_VECTOR[1]);
				}
				++i;
			);
		}
};
