//
// Created by Marco De Mattia on 9/7/15.
//

#ifndef REMOTEPROJECTS_LINEARIZEDTRACKFITTERBASE_H
#define REMOTEPROJECTS_LINEARIZEDTRACKFITTERBASE_H

#include <string>
#include <vector>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GetVariables.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixReader.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/CombinationIndexListBuilder.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/BuildTestFunctions.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/EmulatorTools.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixReaderEmulator.h"

class LinearizedTrackFitterBase
{
 public:
  LinearizedTrackFitterBase(const std::string & baseDir, const bool inputExtrapolateR,
                            const int inputExtrapolatedRPrecision,
                            const bool inputCorrectNonRadialStrips, const int regionsNumber,
                            const std::string & preEstimatePtDirName = "",
                            const std::string & preEstimateCotThetaDirName = "",
                            const std::string & linearFitLowPtDirName = "",
                            const std::string & linearFitHighPtDirName = "",
                            const std::string & linearFitLongitudinalDirName = "",
                            const bool alignPrincipals = true);

  std::vector<double> estimatedPars() { return estimatedPars_; }
  virtual std::vector<double> principalComponents() { return std::vector<double>(varsNum_, 0.); }
  virtual std::vector<double> normalizedPrincipalComponents() { return std::vector<double>(varsNum_, 0.); }
  double chi2Transverse() const { return chi2Transverse_; }
  int ndofTransverse() const { return ndofTransverse_; }
  double chi2Longitudinal() const { return chi2Longitudinal_; }
  int ndofLongitudinal() const { return ndofLongitudinal_; }
  int ndof() const { return (ndofTransverse_+ndofLongitudinal_); }

 protected:

  void fillMatrix(std::unordered_map<unsigned long, EstimatorSimple> * matrices, const unsigned long index,
                  const std::string & fullFileName, const double & deltaPhi, const double & deltaA,
                  const int registerBits, const int bitsPhi, const int bitsA,
                  const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                  const MaxDeltaAndFactors & maxDeltaAndFactors,
                  const double & scaleFactor, const double & ptSplitValue = 10.,
                  const std::vector<bool> & powerTwoRanges = {});
  void fillMatrix(std::unordered_map<unsigned long, EstimatorSimpleEmulator> * matrices, const unsigned long index,
                  const std::string & fullFileName, const double & deltaPhi, const double & deltaA,
                  const int registerBits, const int bitsPhi, const int bitsA,
                  const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                  const MaxDeltaAndFactors & maxDeltaAndFactors,
                  const double & scaleFactor, const double & ptSplitValue = 10.,
                  const std::vector<bool> & powerTwoRanges = {});
  void fillMatrix(std::unordered_map<unsigned long, MatrixReader> * matrices, const unsigned long index,
                  const std::string & fullFileName, const double & deltaPhi, const double & deltaA,
                  const int registerBits, const int bitsPhi, const int bitsA,
                  const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                  const MaxDeltaAndFactors & maxDeltaAndFactors,
                  const double & scaleFactor, const double & ptSplitValue = 10.,
                  const std::vector<bool> & powerTwoRanges = {});
  void fillMatrix(std::unordered_map<unsigned long, MatrixReaderEmulator> * matrices, const unsigned long index,
                  const std::string & fullFileName, const double & deltaPhi, const double & deltaA,
                  const int registerBits, const int bitsPhi, const int bitsA,
                  const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                  const MaxDeltaAndFactors & maxDeltaAndFactors,
                  const double & scaleFactor, const double & ptSplitValue = 10.,
                  const std::vector<bool> & powerTwoRanges = {});
  MaxDeltaAndFactors computeMaxDeltaAndFactorLowHighPt(std::unordered_map<unsigned long, MatrixReaderEmulator> * matricesLowPt,
                                                       std::unordered_map<unsigned long, MatrixReaderEmulator> * matricesHighPt);
  std::string buildFullFileName(const std::string & fileName, const std::string & baseDir, const unsigned long & index);
  void fillMatrices(const std::string & baseDirLowPt, const std::string & fileNameLowPt,
                    std::unordered_map<unsigned long, MatrixReaderEmulator> * matricesLowPt,
                    const std::string & baseDirHighPt, const std::string & fileNameHighPt,
                    std::unordered_map<unsigned long, MatrixReaderEmulator> * matricesHighPt,
                    const double & deltaX = 0., const double & deltaA = 0.,
                    const int registerBits = 0, const int bitsX = 0, const int bitsA = 0,
                    const int maxBitsMultiplyUnitX = 18, const int maxBitsMultiplyUnitY = 27,
                    const double & scaleFactor = 1., const double & ptSplitValue = 10.,
                    const std::vector<bool> & powerTwoRanges = {});

  std::string preEstimatePtDirName_;
  std::string preEstimateCotThetaDirName_;
  std::string linearFitLowPtDirName_;
  std::string linearFitHighPtDirName_;
  std::string linearFitLongitudinalDirName_;
  double ptSplitValue_;
  std::vector<double> principalComponents_;
  std::vector<double> normalizedPrincipalComponents_;
  std::unordered_map<unsigned long, std::vector<double> > meanRadius_;
  unsigned long combinationIndex_;
  std::vector<int> uniqueLayers_;
  unsigned int varsNum_;
  std::string baseDir_;
  CombinationIndexListBuilder combinationIndexListBuilder_;
  bool extrapolateR_;
  int extrapolatedRPrecision_;
  bool correctNonRadialStrips_;
  int regionsNumber_;
  double chi2Transverse_;
  int ndofTransverse_;
  double chi2Longitudinal_;
  int ndofLongitudinal_;
  std::vector<double> estimatedPars_;
  double rotationFactor_;
  bool alignPrincipals_;


  MaxDeltaAndFactors computeMaxDeltaAndFactor(std::unordered_map<unsigned long, MatrixReaderEmulator> * matrices,
                                              const bool excludeBarrelFromNormalizedMatrices) const;

  MaxDeltaAndFactors computeMaxDeltaAndFactor(std::unordered_map<unsigned long, EstimatorSimpleEmulator> * matrices,
                                              const bool excludeBarrelFromNormalizedMatrices) const;

  template <class T>
  MaxDeltaAndFactors computeMaxDeltaAndFactor(std::unordered_map<unsigned long, T> * matrices,
                                              const bool excludeBarrelFromNormalizedMatrices) const
  {
    return MaxDeltaAndFactors();
  }


  template <class T>
  void fillMatrices(const std::string & baseDir, const std::string & fileName,
                    std::unordered_map<unsigned long, T> * matrices,
                    const double & deltaX = 0., const double & deltaA = 0.,
                    const int registerBits = 0, const int bitsX = 0, const int bitsA = 0,
                    const int maxBitsMultiplyUnitX = 18, const int maxBitsMultiplyUnitY = 27,
                    const double & scaleFactor = 1., const double & ptSplitValue = 10.,
                    const bool normalizeMatrices = false, const bool excludeBarrelFromNormalizedMatrices = false,
                    const std::vector<bool> & powerTwoRanges = {})
  {
    bool fiveOutOfSix = true;

    std::vector<unsigned long> combinationIndexList;
    combinationIndexListBuilder_.fillDefaultIndexList(combinationIndexList, fiveOutOfSix, regionsNumber_);

    for (auto index : combinationIndexList) {
      try {
        std::string fullFileName(fileName);
        fullFileName.replace(fullFileName.find("0"), 1, std::to_string(index));
        fullFileName = baseDir + "/" + fullFileName;
        fillMatrix(matrices, index, fullFileName, deltaX, deltaA, registerBits, bitsX, bitsA,
                   maxBitsMultiplyUnitX, maxBitsMultiplyUnitY, MaxDeltaAndFactors(), scaleFactor, ptSplitValue,
                   powerTwoRanges);
      }
      catch (int exception) {
        std::cout << "Error: Matrix for combination = " << index << " not found" << std::endl;
        throw;
      }
    }

    // Normalize matrices so that each row is encoded with the same range (it can be different for different rows)
    if (normalizeMatrices) {
      MaxDeltaAndFactors maxDeltaAndFactors(computeMaxDeltaAndFactor(matrices, excludeBarrelFromNormalizedMatrices));
      matrices->clear();
      for (auto index : combinationIndexList) {
        try {
          std::string fullFileName(fileName);
          fullFileName.replace(fullFileName.find("0"), 1, std::to_string(index));
          fullFileName = baseDir + "/" + fullFileName;
          fillMatrix(matrices, index, fullFileName, deltaX, deltaA, registerBits, bitsX, bitsA,
                     maxBitsMultiplyUnitX, maxBitsMultiplyUnitY, maxDeltaAndFactors, scaleFactor, ptSplitValue);
        }
        catch (int exception) {
          std::cout << "Error: Matrix for combination = " << index << " unable to normalize" << std::endl;
          throw;
        }
      }
    }

    bool writeEncodedCoefficients = false;
    if (writeEncodedCoefficients) {
      for (auto index : combinationIndexList) {
        try {
          matrices->find(index)->second.write();
        }
        catch (int exception) {
          std::cout << "Error: Matrix for combination = " << index << " unable to write" << std::endl;
          throw;
        }
      }
    }
  }


  template <class T>
  void consistencyCheck(const std::vector<T> & vars, const std::vector<int> & layers)
  {
    if (vars.size() < 15) {
      std::cout << "Error: number of input variables is less than 15. Please provide 5 or 6 sets of (phi, R, z) ordered from the innermost to the outermost layer." << std::endl;
      std::cout << "Number of input variables = " << vars.size() << std::endl;
      throw;
    }
    if (layers.size()*3 != vars.size()) {
      std::cout << "Error: inconsistent number of layers and number of variables. They should be in a ratio of 1/3." << std::endl;
      std::cout << "Number of layers = " << layers.size() << std::endl;
      for (auto l : layers) std::cout << l << ", ";
      std::cout << std::endl;
      std::cout << "Number of variables = " << vars.size() << std::endl;
      throw;
    }
  }


  /// This is used by the full simulation
  template <class T>
  void fillVariablesFromBits(const std::vector<T> & vars, const int bits,
                             std::vector<T> & cleanedVars, std::vector<int> & layers)
  {
    if (bits == 0) layers = {5, 6, 7, 8, 9, 10};
    else if (bits == 1) layers = {6, 7, 8, 9, 10};
    else if (bits == 2) layers = {5, 7, 8, 9, 10};
    else if (bits == 3) layers = {5, 6, 8, 9, 10};
    else if (bits == 4) layers = {5, 6, 7, 9, 10};
    else if (bits == 5) layers = {5, 6, 7, 8, 10};
    else if (bits == 6) layers = {5, 6, 7, 8, 9};
    else {
      std::cout << "Error: unknown bits = " << bits << std::endl;
      throw;
    }

    // Clean the variables removing the 0 values corresponding to the missing layer
    if (bits > 0) {
      cleanedVars.clear();
      // std::vector<T> cleanedVars;
      for (size_t i = 0; i < vars.size() / 3; ++i) {
        if (i != size_t(bits - 1)) {
          cleanedVars.push_back(vars.at(i * 3));
          cleanedVars.push_back(vars.at(i * 3 + 1));
          cleanedVars.push_back(vars.at(i * 3 + 2));
        }
      }
    }
    else {
      // We should avoid this copy
      cleanedVars = vars;
    }
  }


  /// Align the principal components when there are missing variables
  template <class T>
  void alignPrincipals(std::vector<T> & principals, const int nDof)
  {
    for (int i=0; i<4-nDof; ++i) {
      principals.insert(principals.begin(), 0);
    }
  }
};



#endif //REMOTEPROJECTS_LINEARIZEDTRACKFITTERBASE_H
