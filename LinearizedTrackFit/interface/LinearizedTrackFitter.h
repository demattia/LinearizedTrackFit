//
// Created by Marco De Mattia on 4/14/15.
//

#ifndef REMOTEPROJECTS_LINEARIZEDTRACKFITTER_H
#define REMOTEPROJECTS_LINEARIZEDTRACKFITTER_H

#include <vector>
#include <memory>
#include <bitset>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearizedTrackFitterBase.h"
#include <TH2D.h>
#include <TFile.h>

class LinearizedTrackFitter : public LinearizedTrackFitterBase
{
 public:
  LinearizedTrackFitter(const std::string & baseDir, const bool inputExtrapolateR,
                        const int inputExtrapolatedRPrecision,
                        const bool inputCorrectNonRadialStrips, const int regionsNumber,
                        const std::string & preEstimatePtDirName = "",
                        const std::string & preEstimateCotThetaDirName = "",
                        const std::string & linearFitLowPtDirName = "",
                        const std::string & linearFitHighPtDirName = "",
                        const std::string & linearFitLongitudinalDirName = "",
                        const bool alignPrincipals = true);
  virtual ~LinearizedTrackFitter() = default;
  virtual double fit(const std::vector<double> & vars, const std::vector<int> & layers);
  virtual double fit(const std::vector<double> & vars, const std::vector<int> & layers,
                     const std::vector<int> & stripIndexes) { return fit(vars, layers); }
  double fit(const std::vector<double> & vars, const int bits);
  virtual std::vector<double> principalComponents();
  virtual std::vector<double> normalizedPrincipalComponents();

  std::vector<double> getFirstOrderTerm() { return firstOrderTerm_; }
  std::vector<double> getSecondOrderTerm1() { return secondOrderTerm1_; }
  std::vector<double> getSecondOrderTerm2() { return secondOrderTerm2_; }
  std::vector<double> getSecondOrderTerm3() { return secondOrderTerm3_; }

 protected:
  virtual double fit(const double & chargeOverTwoRho, const double & cotTheta, const double & tgTheta);

  double preEstimatedPt_;
  std::vector<double> varsR_;
  std::vector<double> extrapolatedR_;
  Matrix<long double, Dynamic, 1> correctedVarsPhi_;
  Matrix<long double, Dynamic, 1> correctedVarsZ_;
  std::unordered_map<unsigned long, EstimatorSimple> chargeOverPtEstimator_;
  std::unordered_map<unsigned long, EstimatorSimple> cotThetaEstimator_;
  std::unordered_map<unsigned long, EstimatorSimple> tgThetaEstimator_;
  std::unordered_map<unsigned long, MatrixReader> linearFitLowPt_;
  std::unordered_map<unsigned long, MatrixReader> linearFitHighPt_;
  std::unordered_map<unsigned long, MatrixReader> linearFitLongitudinal_;
  CorrectPhiForNonRadialStripsLookup correctPhiForNonRadialStripsLookup_;
  std::vector<double> firstOrderTerm_;
  std::vector<double> secondOrderTerm1_;
  std::vector<double> secondOrderTerm2_;
  std::vector<double> secondOrderTerm3_;
};

#endif //REMOTEPROJECTS_LINEARIZEDTRACKFITTER_H
