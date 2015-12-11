#ifndef gmm_h_
#define gmm_h_

#include <iostream>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <fstream>

#include <cfloat>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <climits>

const static int MAX_ITERATOR = 1000;
const static double END_THR = 0.0001;
const static double SIM_THR = 0.2;
const static double W_THR = 0.1;

struct Gaussian{
    double mean, delta;
    double weight;
    Gaussian(double m=0, double v=0, double w=1.0): mean(m), delta(v), weight(w){
    }
    double getProbability(double x) const {
        return weight * std::pow(M_E, -std::pow(x-mean, 2.0) / (2*delta*delta)) / ( std::pow(2*M_E, 0.5) * delta );
    }
    private:
    friend std::ostream& operator<<(std::ostream& os, const Gaussian & x);
};


class GMM {
public:
    void gmm(const std::vector<double> & data, int mxCenter, std::vector< Gaussian > &re, double DELTA) ;
private:
    struct comp{
        std::pair<double, double> operator()(const std::pair<double, double> &a, double x) {
            return std::make_pair(a.first + x, a.second + x*x);
        }
    };

    Gaussian getGaussian(const std::vector<double> & data) ;
    double getDelta(const std::vector<double> & data) ;
    double fixCenterGmm(const std::vector<double> & data, int centers, std::vector< Gaussian > &re, double DELTA = 0.0) ;
    bool ok(const std::vector< Gaussian >& re, const std::vector< Gaussian >& tmp) ;
    double caculateBIC(const std::vector<double> &data, const std::vector< Gaussian >& gau) ;
};
#endif //end gmm_h_
