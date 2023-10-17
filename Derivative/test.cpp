#include "pch.h"
#include "../Derivative/derivative.hpp"

TEST(TestCaseDerivative, Coef5) {
    std::vector<double> points5{ -2, -1, 1, 2, 3 };
    DerivativeCoef<double> dcoef5;
    double centralcoef, sumcoef = 0;

    dcoef5 = calcDerivativeCoef<double>(points5);
    for (int i = 0; i < points5.size(); i++) {
        sumcoef += dcoef5.otherCoefs[i];
    }

    sumcoef = round(sumcoef * 1000) / 1000;               // округление до 3х знаков после запятой
    centralcoef = round(dcoef5.centralCoef * 1000) / 1000;
    	
	EXPECT_EQ(centralcoef, -sumcoef);
}