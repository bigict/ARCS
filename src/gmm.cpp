#include "gmm.h"

    std::ostream& operator<<(std::ostream& os, const Gaussian & x) {
        os << "mean: " << x.mean << " delta: " << x.delta << " weight: " << x.weight;
        return os;
    }

    void GMM::gmm(const std::vector<double> & data, int mxCenter, std::vector< Gaussian > &re, double DELTA) {
        double BIC = DBL_MAX;
        std::vector< Gaussian > tmpResult;
        for(int i = 1; i <= mxCenter; ++i) {
            std::vector< Gaussian > tmp;
            double newBIC = fixCenterGmm(data, i, tmp, DELTA);
            if( newBIC < BIC) {
                BIC = newBIC;
                tmpResult = tmp; 
            }
        }
        for(int i = 0; i < tmpResult.size(); ++i) {
            bool ok = true;
            if( tmpResult[i].weight < W_THR ) {
                continue;
            }
            for(int j = i+1; j < tmpResult.size(); ++j) {
                if( fabs(tmpResult[i].mean - tmpResult[j].mean) < tmpResult[i].mean * SIM_THR) {
                    ok = false;
                    tmpResult[j].weight += tmpResult[i].weight;
                    break;
                }
            }
            if(ok) {
                re.push_back(tmpResult[i]);
            }
        }
        return ;
    }

    Gaussian GMM::getGaussian(const std::vector<double> & data) {
        std::pair<double, double> re = accumulate(data.begin(), data.end(), std::make_pair(0.0, 0.0), comp());
        return Gaussian(re.first / data.size(), std::pow( re.second / data.size() - std::pow(re.first / data.size(), 2.0), 0.5), 1.0);
    }

    double GMM::getDelta(const std::vector<double> & data) {
        std::pair<double, double> re = accumulate(data.begin(), data.end(), std::make_pair(0.0, 0.0), comp());
        return std::pow( re.second / data.size() - std::pow(re.first / data.size(), 2.0), 0.5);
    }

    double GMM::fixCenterGmm(const std::vector<double> & data, int centers, std::vector< Gaussian > &re, double DELTA ) {
        if( centers <= 1 ) {
            re.push_back( getGaussian(data) );
            return caculateBIC(data, re);
        }
        bool fixDelta = true;
        if(DELTA <= 0.01) {
            fixDelta = false;
        }
        double mx = *max_element(data.begin(), data.end());
        double mn = *min_element(data.begin(), data.end());
        double diff = mx - mn;
        double delta = getDelta(data);
        for(int i = 0; i < centers; ++i) {
            if(fixDelta) {
                re.push_back( Gaussian(mn + i*diff/(centers-1), DELTA, 1.0 / centers) );
            } else {
                re.push_back( Gaussian(mn + i*diff/(centers-1), delta, 1.0 / centers) );
            }
        }
        std::vector< std::vector<double> > beta( data.size(), std::vector<double>(centers, 0.0) );
        std::vector< Gaussian > tmp(centers, Gaussian() );
        int itera = 0;
        while( itera++ < MAX_ITERATOR && !ok(re, tmp) ) {
            tmp = re;
            for(int i = 0; i < data.size(); ++i) {
                for(int j = 0; j < centers; ++j) {
                    beta[i][j] = re[j].getProbability(data[i]);
                }
                double sum = accumulate(beta[i].begin(), beta[i].end(), 0.0);
                for(int j = 0; j < centers; ++j) {
                    beta[i][j] /= sum;
                }
            }
            for(int j = 0; j < centers; ++j) {
                double sumBeta = 0.0, sumweightBeta = 0.0, sumVar = 0.0;
                for(int i = 0; i < data.size(); ++i) {
                    sumBeta += beta[i][j];
                    sumweightBeta += data[i] * beta[i][j];
                }
                re[j].weight = sumBeta / data.size();
                re[j].mean = sumweightBeta / sumBeta;
                for(int i = 0; i < data.size(); ++i) {
                    sumVar += beta[i][j] * (data[i] - re[j].mean) * (data[i] - re[j].mean);
                }
                if(!fixDelta) {
                    re[j].delta = std::pow( sumVar / sumBeta, 0.5);
                }
            }
        }
        return caculateBIC(data, re);
    }

    bool GMM::ok(const std::vector< Gaussian >& re, const std::vector< Gaussian >& tmp) {
        double diff = 0.0;
        double sum = 0.0;
        for(int i = 0; i < re.size(); ++i) {
            diff += fabs( re[i].mean - tmp[i].mean );
            sum += re[i].mean;
        }
        return diff / sum < END_THR;
    }

    double GMM::caculateBIC(const std::vector<double> &data, const std::vector< Gaussian >& gau) {
        double BIC = (2 * gau.size() ) * log( data.size() );
        for(int i = 0; i < data.size(); ++i) {
            double pro = 0.0;
            for(int j = 0; j < gau.size(); ++j) {
                pro += gau[j].getProbability(data[i]);
            }
            BIC -= 2*log(pro);
        }
        return BIC;
    }
