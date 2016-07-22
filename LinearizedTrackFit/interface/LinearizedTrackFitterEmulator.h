//
// Created by Marco De Mattia on 9/2/15.
//

#ifndef REMOTEPROJECTS_LINEARIZEDTRACKFITTEREMULATOR_H
#define REMOTEPROJECTS_LINEARIZEDTRACKFITTEREMULATOR_H

#include <vector>
#include <memory>
#include <bitset>
// #include <boost/multiprecision/cpp_int.hpp>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GetVariables.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixReader.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/CombinationIndexListBuilder.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/BuildTestFunctions.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearizedTrackFitterBase.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixReaderEmulator.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/EmulatorTools.h"

class LinearizedTrackFitterEmulator: public LinearizedTrackFitterBase
{
 public:
  LinearizedTrackFitterEmulator(const std::string & baseDir, const bool inputExtrapolateR,
                                const int inputExtrapolatedRPrecision,
                                const bool inputCorrectNonRadialStrips, const int regionsNumber,
                                const int registerBits = 48, const int bitsX = 27, const int bitsA = 18,
                                const double & deltaPhi = 2., const double & deltaR = 1024., const double & deltaZ = 1024.,
                                const int maxBitsMultiplyUnitX = 27, const int maxBitsMultiplyUnitY = 18,
                                const std::string & preEstimatePtDirName = "",
                                const std::string & preEstimateCotThetaDirName = "",
                                const std::string & linearFitLowPtDirName = "",
                                const std::string & linearFitHighPtDirName = "",
                                const std::string & linearFitLongitudinalDirName = "",
                                const int reducedAccuracyBitsPhi = 18,
                                const int reducedAccuracyBitsR = 18,
                                const int reducedAccuracyBitsZ = 18,
                                const int reducedAccuracyBitsStripIndex = 3,
                                const bool normalizeMatrices = true,
                                const bool commonLowHighPtFix = true,
                                const bool powerTwoRangeChargeOverPt = false,
                                const bool saveOutput = false,
                                const bool alignPrincipals = true);

  virtual ~LinearizedTrackFitterEmulator() = default;

  double fit(const std::vector<double> & vars, const std::vector<int> & layers,
             const std::vector<int> & stripIndexes);
  double fit(const std::vector<bigInt> & vars, const std::vector<int> & layers,
             const std::vector<int> & stripIndexes);
  virtual std::vector<double> normalizedPrincipalComponents();

  template <class T>
  double fit(const std::vector<T> & vars, const std::vector<int> & stripIndexes, const int bits)
  {
    std::vector<T> cleanedVars;
    std::vector<int> layers;
    fillVariablesFromBits(vars, bits, cleanedVars, layers);
    return fit(cleanedVars, layers, stripIndexes);
  }

  std::vector<bigInt> estimatedParsInt() { return estimatedParsInt_; }
  std::vector<bigInt> chi2TermsInt() { return chi2TermsInt_; }

 private:
  virtual double fit(const bigInt & chargeOverTwoRho, const bigInt & cotTheta, const bigInt & tgTheta);
  double fit();

  // typedef boost::multiprecision::int128_t bigInt;
  // typedef long long int bigInt;
  // typedef int64_t bigInt;
  std::vector<bigInt> varsPhiInt_;
  std::vector<bigInt> varsRInt_;
  std::vector<bigInt> varsZInt_;
  std::vector<bigInt> extrapolatedRInt_;
  std::vector<bigInt> correctedVarsPhiInt_;
  std::vector<bigInt> correctedVarsZInt_;
  std::vector<bigInt> stripCorrectionTermsInt_;
  std::unordered_map<unsigned long, EstimatorSimpleEmulator > chargeOverPtEstimator_;
  std::unordered_map<unsigned long, EstimatorSimpleEmulator > cotThetaEstimator_;
  std::unordered_map<unsigned long, EstimatorSimpleEmulator > tgThetaEstimator_;
  std::unordered_map<unsigned long, MatrixReaderEmulator> linearFitLowPt_;
  std::unordered_map<unsigned long, MatrixReaderEmulator> linearFitHighPt_;
  std::unordered_map<unsigned long, MatrixReaderEmulator> linearFitLongitudinal_;
  std::unordered_map<unsigned long, std::vector<bigInt> > meanRadiusInt_;
  int registerBits_;
  int bitsX_;
  int bitsA_;
  int shiftBitsPhi_;
  int shiftBitsR_;
  int shiftBitsZ_;
  double deltaPhi_;
  int reductionPhi_;
  double deltaR_;
  double deltaZ_;
  int reductionZ_;
  int reductionPt_;
  int reductionCotTheta_;
  int reductionTgTheta_;
  int maxBitsMultiplyUnitX_;
  int maxBitsMultiplyUnitY_;
  int reducedAccuracyBitsStripIndex_;
  bigInt oneOverTwoRhoSplitValueInt_;
  bigInt Rcut_;
  bigInt RadiusCut1_;
  bigInt RadiusCut2_;
  CorrectPhiForNonRadialStripsEmulator correctPhiForNonRadialStripsEmulator_;
  std::vector<bigInt> oneOverRInt_;
  std::vector<bigInt> oneOverRSquaredInt_;
  std::vector<bigInt> stripIndex_;
  bigInt preChargeOverTwoRho_;
  std::vector<bigInt> estimatedParsInt_;
  std::vector<bigInt> chi2TermsInt_;

  // Output files to store variables for input to the firmware and output for comparison with the firmware
  bool saveOutput_;
  std::ofstream varOut_;
  std::ofstream emulatorOutput_;


  template <class T>
  void initialize(const std::vector<T> & vars, const std::vector<int> & layers, const std::vector<int> & stripIndexes)
  {
    consistencyCheck(vars, layers);
    varsNum_ = vars.size()/3;
    varsPhiInt_.clear();
    varsPhiInt_.reserve(varsNum_);
    varsRInt_.clear();
    varsRInt_.reserve(varsNum_);
    varsZInt_.clear();
    varsZInt_.reserve(varsNum_);

    stripCorrectionTermsInt_.clear();
    stripCorrectionTermsInt_.reserve(varsNum_);
    // Make sure that the default strip corrections are zero.
    std::fill(stripCorrectionTermsInt_.begin(), stripCorrectionTermsInt_.end(), 0);

    oneOverRInt_.clear();
    oneOverRInt_.reserve(varsNum_);
    oneOverRSquaredInt_.clear();
    oneOverRSquaredInt_.reserve(varsNum_);

    uniqueLayers_ = layers;
    std::sort(uniqueLayers_.begin(), uniqueLayers_.end());
    uniqueLayers_.erase(std::unique(uniqueLayers_.begin(), uniqueLayers_.end()), uniqueLayers_.end());

    if (stripIndexes.size() > 6) {
      std::cout << "Error: strip indexes vector size is greater than 6. Does it contain duplicates?" << std::endl;
      throw;
    }
  }
};

#endif //REMOTEPROJECTS_LINEARIZEDTRACKFITTEREMULATOR_H
