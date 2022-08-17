#pragma once
#include <vector>
#include <string>
struct BreakPointData
{
	int BreakChromosome;
	int BreakIdx;
	int JoinChromosome;
	int JoinIdx;
	int Coverage;


	BreakPointData(){};
	BreakPointData(std::vector<std::string> input)
	{
		BreakChromosome = std::stoi(input[1]);
		BreakIdx = std::stoi(input[2]);
		JoinChromosome = std::stoi(input[3]);
		JoinIdx = std::stoi(input[4]);
		Coverage = std::stoi(input[5]);
	}
};