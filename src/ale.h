#pragma once
#include <cmath>

double ale(double x, double y)
{
	if (x > y)
	{
		return x + log(1.0 + exp(y - x));
	}
	else
	{
		return y + log(1.0 + exp(x-y));
	}
}