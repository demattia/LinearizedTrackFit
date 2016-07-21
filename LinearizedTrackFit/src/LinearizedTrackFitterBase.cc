//
// Created by Marco De Mattia on 9/7/15.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearizedTrackFitterBase.h"

LinearizedTrackFitterBase::LinearizedTrackFitterBase(const std::string & baseDir, const bool inputExtrapolateR,
                                                     const int inputExtrapolatedRPrecision,
                                                     const bool inputCorrectNonRadialStrips, const int regionsNumber,
                                                     const std::string & preEstimatePtDirName,
                                                     const std::string & preEstimateCotThetaDirName,
                                                     const std::string & linearFitLowPtDirName,
                                                     const std::string & linearFitHighPtDirName,
                                                     const std::string & linearFitLongitudinalDirName,
                                                     const bool alignPrincipals) :
    preEstimatePtDirName_(preEstimatePtDirName),
    preEstimateCotThetaDirName_(preEstimateCotThetaDirName),
    linearFitLowPtDirName_(linearFitLowPtDirName),
    linearFitHighPtDirName_(linearFitHighPtDirName),
    linearFitLongitudinalDirName_(linearFitLongitudinalDirName),
    ptSplitValue_(10.),
    combinationIndex_(0),
    baseDir_(baseDir),
    extrapolateR_(inputExtrapolateR),
    extrapolatedRPrecision_(inputExtrapolatedRPrecision),
    correctNonRadialStrips_(inputCorrectNonRadialStrips),
    regionsNumber_(regionsNumber),
    chi2Transverse_(-1.),
    ndofTransverse_(-1),
    chi2Longitudinal_(-1.),
    ndofLongitudinal_(-1),
    rotationFactor_(0.),
    alignPrincipals_(alignPrincipals)
{
  if (extrapolatedRPrecision_ < 0 || extrapolatedRPrecision_ > 4) {
    std::cout << "Error: extrapolated R precision can only be a number between 0 and 4 (included). ";
    std::cout << "Extrapolated R precision requested = " << extrapolatedRPrecision_ << std::endl;
  }

  if (correctNonRadialStrips_) extrapolateR_ = true;
  if (linearFitLowPtDirName_ == "") {
    preEstimatePtDirName_ = baseDir_+"/FlatPt/PreEstimate_Transverse";
    preEstimateCotThetaDirName_ = baseDir_+"/FlatPt/PreEstimate_Longitudinal_Rz";
    if (correctNonRadialStrips_) {
      linearFitLowPtDirName_ = baseDir_+"/FlatOneOverPt/Combinations_FullCorrections_2_10_SecondFirst_flatPtTgThetaAndPtPreEstimate";
      linearFitHighPtDirName_ = baseDir_+"/FlatPt/Combinations_FullCorrections_10_more_SecondFirst_flatPtTgThetaAndPtPreEstimate";
      linearFitLongitudinalDirName_ = baseDir_+"/FlatPt/Combinations_Longitudinal_Rz_flatPtTgThetaAndPtPreEstimate";
    }
    else if (extrapolateR_) {
      linearFitLowPtDirName_ = baseDir_+"/Combinations_Transverse_SecondOrder_ExtrapolatedRSecondOrder_2_10";
      linearFitHighPtDirName_ = baseDir_+"/Combinations_Transverse_SecondOrder_ExtrapolatedRSecondOrder_10_more";
    }
    else {
      linearFitLowPtDirName_ = baseDir_+"/NineRegions/Combinations_OldBaseline_2_10";
      linearFitHighPtDirName_ = baseDir_+"/NineRegions/Combinations_OldBaseline_10_more";
      linearFitLongitudinalDirName_ = baseDir_+"/NineRegions/Combinations_Longitudinal_Rz_SecondOrder";
    }
  }
}


void LinearizedTrackFitterBase::fillMatrix(std::unordered_map<unsigned long, EstimatorSimple> * matrices,
                                           const unsigned long index, const std::string & fullFileName,
                                           const double & deltaPhi, const double & deltaR, const int registerBits,
                                           const int bitsPhi, const int bitsA,
                                           const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                                           const MaxDeltaAndFactors & maxDeltaAndFactors,
                                           const double & scaleFactor, const double & ptSplitValue,
                                           const std::vector<bool> & powerTwoRanges)
{
  matrices->insert(std::make_pair(index, EstimatorSimple(fullFileName, scaleFactor)));
}


void LinearizedTrackFitterBase::fillMatrix(std::unordered_map<unsigned long, EstimatorSimpleEmulator> * matrices,
                                           const unsigned long index, const std::string & fullFileName,
                                           const double & deltaPhi, const double & deltaR, const int registerBits,
                                           const int bitsPhi, const int bitsA,
                                           const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                                           const MaxDeltaAndFactors & maxDeltaAndFactors,
                                           const double & scaleFactor, const double & ptSplitValue,
                                           const std::vector<bool> & powerTwoRanges)
{
  matrices->insert(std::make_pair(index, EstimatorSimpleEmulator(fullFileName, deltaPhi, deltaR, registerBits, bitsPhi,
                                                                 bitsA, maxBitsMultiplyUnitX, maxBitsMultiplyUnitY,
                                                                 maxDeltaAndFactors.maxDeltaD, maxDeltaAndFactors.maxFactorD,
                                                                 scaleFactor, ptSplitValue)));
}


void LinearizedTrackFitterBase::fillMatrix(std::unordered_map<unsigned long, MatrixReader> * matrices,
                                           const unsigned long index, const std::string & fullFileName,
                                           const double & deltaPhi, const double & deltaR, const int registerBits,
                                           const int bitsPhi, const int bitsA,
                                           const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                                           const MaxDeltaAndFactors & maxDeltaAndFactors,
                                           const double & scaleFactor, const double & ptSplitValue,
                                           const std::vector<bool> & powerTwoRanges)
{
  matrices->insert(std::make_pair(index, MatrixReader(fullFileName)));
}


void LinearizedTrackFitterBase::fillMatrix(std::unordered_map<unsigned long, MatrixReaderEmulator> * matrices,
                                           const unsigned long index, const std::string & fullFileName,
                                           const double & deltaPhi, const double & deltaR, const int registerBits,
                                           const int bitsPhi, const int bitsA,
                                           const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                                           const MaxDeltaAndFactors & maxDeltaAndFactors,
                                           const double & scaleFactor, const double & ptSplitValue,
                                           const std::vector<bool> & powerTwoRanges)
{
  matrices->insert(std::make_pair(index, MatrixReaderEmulator(fullFileName, deltaPhi, registerBits, bitsPhi, bitsA,
                                                              maxBitsMultiplyUnitX, maxBitsMultiplyUnitY,
                                                              maxDeltaAndFactors, powerTwoRanges)));
}


template <>
MaxDeltaAndFactors LinearizedTrackFitterBase::computeMaxDeltaAndFactor(std::unordered_map<unsigned long, MatrixReaderEmulator> * matrices,
                                                                       const bool excludeBarrelFromNormalizedMatrices) const
{
  std::vector<double> minDeltaD(2, 9999999.);
  std::vector<double> maxDeltaD(2, 0.);
  std::vector<int> minFactorD(2, 9999999);
  std::vector<int> maxFactorD(2, 0);

  std::vector<double> minDeltaV(6, 9999999.);
  std::vector<double> maxDeltaV(6, 0.);
  std::vector<int> minFactorV(6, 9999999);
  std::vector<int> maxFactorV(6, 0);

  for (auto it = matrices->begin(); it != matrices->end(); ++it) {
    size_t vSize = it->second.getDeltaVSize();
    int shift = 0;
    if (vSize == 5) shift = 1;
    for (size_t i = 0; i < vSize; ++i) {
      if (minDeltaV[i + shift] > it->second.getDeltaV(i)) {
        minDeltaV[i + shift] = it->second.getDeltaV(i);
        minFactorV[i + shift] = it->second.getFactorV(i);
      }
      if (maxDeltaV[i + shift] < it->second.getDeltaV(i)) {
        maxDeltaV[i + shift] = it->second.getDeltaV(i);
        maxFactorV[i + shift] = it->second.getFactorV(i);
      }
    }
//    for (size_t i=0; i<it->second.getDeltaDSize(); ++i) {
    for (size_t i = 0; i < 2; ++i) {
      if (minDeltaD[i] > it->second.getDeltaD(i)) {
        minDeltaD[i] = it->second.getDeltaD(i);
        minFactorD[i] = it->second.getFactorD(i);
      }
      if (maxDeltaD[i] < it->second.getDeltaD(i)) {
        maxDeltaD[i] = it->second.getDeltaD(i);
        maxFactorD[i] = it->second.getFactorD(i);
      }
    }
  }
//    for (int i = 0; i < 2; ++i) {
//      std::cout << "Maximum variation of bit-shift factor for parameter " << i << " = " <<
//      maxFactorD[i] - minFactorD[i] << std::endl;
//    }
//    for (int i = 0; i < 6; ++i) {
//      std::cout << "Maximum variation of bit-shift factor for chi2 term " << i << " = " <<
//      maxFactorV[i] - minFactorV[i] << std::endl;
//    }

  // Using move semantics in the constructor since the max vectors go out of scope
  return MaxDeltaAndFactors(maxDeltaD, maxFactorD, maxDeltaV, maxFactorV);
}


template <>
MaxDeltaAndFactors LinearizedTrackFitterBase::computeMaxDeltaAndFactor(std::unordered_map<unsigned long, EstimatorSimpleEmulator> * matrices,
                                                                       const bool excludeBarrelFromNormalizedMatrices) const
{
  // double minDeltaA = 9999999.;
  double maxDeltaA = 0.;
  // int minFactorA = 9999999;
  int maxFactorA = 0;
  for (auto it = matrices->begin(); it != matrices->end(); ++it) {
    unsigned long combIndex = it->first;
    // Exclude the barrel if requested. This is useful to avoid the barrel coefficients for the pre-estimate of
    // tgTheta to affect the range of the coefficients for the forward. The pre-tgTheta is not used in the barrel.
    if (excludeBarrelFromNormalizedMatrices && (combIndex == 2016 || combIndex == 1504 || combIndex == 1760 ||
                                                combIndex == 1888 || combIndex == 1952 || combIndex == 1984 || combIndex == 992)) {
      continue;
    }

//    if (minDeltaA > it->second.deltaA()) {
//      minDeltaA = it->second.deltaA();
//      minFactorA = it->second.reduction();
//    }
    if (maxDeltaA < it->second.deltaA()) {
      maxDeltaA = it->second.deltaA();
      maxFactorA = it->second.reduction();
    }
  }
//    std::cout << "Maximum variation of bit-shift factor for pre-estimate = " << maxFactorA - minFactorA << std::endl;
  // Using move semantics in the constructor since the max vectors go out of scope
  return MaxDeltaAndFactors(maxDeltaA, maxFactorA);
}


MaxDeltaAndFactors LinearizedTrackFitterBase::computeMaxDeltaAndFactor(std::unordered_map<unsigned long, MatrixReaderEmulator> * matricesLowPt,
                                                                       std::unordered_map<unsigned long, MatrixReaderEmulator> * matricesHighPt)
{
  MaxDeltaAndFactors maxDeltaAndFactorsLowPt(computeMaxDeltaAndFactor(matricesLowPt, false));
  MaxDeltaAndFactors maxDeltaAndFactorsHighPt(computeMaxDeltaAndFactor(matricesHighPt, false));
  std::vector<double> maxDeltaD(2, 0.);
  std::vector<int> maxFactorD(2, 0);
  std::vector<double> maxDeltaV(6, 0.);
  std::vector<int> maxFactorV(6, 0);

  for (size_t i=0; i<maxDeltaD.size(); ++i) {
    // std::cout << "maxDeltaDLowPt = " << maxDeltaAndFactorsLowPt.maxDeltaD.at(i) << std::endl;
    // std::cout << "maxFactorDLowPt = " << maxDeltaAndFactorsLowPt.maxFactorD.at(i) << std::endl;
    // std::cout << "maxDeltaDHighPt = " << maxDeltaAndFactorsHighPt.maxDeltaD.at(i) << std::endl;
    // std::cout << "maxFactorDHighPt = " << maxDeltaAndFactorsHighPt.maxFactorD.at(i) << std::endl;
    if (maxDeltaAndFactorsLowPt.maxDeltaD.at(i) > maxDeltaAndFactorsHighPt.maxDeltaD.at(i)) {
      maxDeltaD.at(i) = maxDeltaAndFactorsLowPt.maxDeltaD.at(i);
      maxFactorD.at(i) = maxDeltaAndFactorsLowPt.maxFactorD.at(i);
    }
    else {
      maxDeltaD.at(i) = maxDeltaAndFactorsHighPt.maxDeltaD.at(i);
      maxFactorD.at(i) = maxDeltaAndFactorsHighPt.maxFactorD.at(i);
    }
  }
  for (size_t i=0; i<maxDeltaV.size(); ++i) {
    // std::cout << "maxDeltaVLowPt = " << maxDeltaAndFactorsLowPt.maxDeltaV.at(i) << std::endl;
    // std::cout << "maxFactorVLowPt = " << maxDeltaAndFactorsLowPt.maxFactorV.at(i) << std::endl;
    // std::cout << "maxDeltaVHighPt = " << maxDeltaAndFactorsHighPt.maxDeltaV.at(i) << std::endl;
    // std::cout << "maxFactorVHighPt = " << maxDeltaAndFactorsHighPt.maxFactorV.at(i) << std::endl;
    if (maxDeltaAndFactorsLowPt.maxDeltaV.at(i) > maxDeltaAndFactorsHighPt.maxDeltaV.at(i)) {
      maxDeltaV.at(i) = maxDeltaAndFactorsLowPt.maxDeltaV.at(i);
      maxFactorV.at(i) = maxDeltaAndFactorsLowPt.maxFactorV.at(i);
    }
    else {
      maxDeltaV.at(i) = maxDeltaAndFactorsHighPt.maxDeltaV.at(i);
      maxFactorV.at(i) = maxDeltaAndFactorsHighPt.maxFactorV.at(i);
    }
  }

  return MaxDeltaAndFactors(maxDeltaD, maxFactorD, maxDeltaV, maxFactorV);
}


std::string LinearizedTrackFitterBase::buildFullFileName(const std::string & fileName, const std::string & baseDir,
                                                         const unsigned long & index)
{
  std::string fullFileName(fileName);
  fullFileName.replace(fullFileName.find("0"), 1, std::to_string(index));
  fullFileName = baseDir + "/" + fullFileName;
  return fullFileName;
}


void LinearizedTrackFitterBase::fillMatrices(const std::string & baseDirLowPt, const std::string & fileNameLowPt,
                                             std::unordered_map<unsigned long, MatrixReaderEmulator> * matricesLowPt,
                                             const std::string & baseDirHighPt, const std::string & fileNameHighPt,
                                             std::unordered_map<unsigned long, MatrixReaderEmulator> * matricesHighPt,
                                             const double & deltaX, const double & deltaA,
                                             const int registerBits, const int bitsX, const int bitsA,
                                             const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                                             const double & scaleFactor, const double & ptSplitValue,
                                             const std::vector<bool> & powerTwoRanges)
{
  bool fiveOutOfSix = true;

  std::vector<unsigned long> combinationIndexList;
  combinationIndexListBuilder_.fillDefaultIndexList(combinationIndexList, fiveOutOfSix, regionsNumber_);

  for (auto index : combinationIndexList) {
    try {
      fillMatrix(matricesLowPt, index, buildFullFileName(fileNameLowPt, baseDirLowPt, index), deltaX, deltaA,
                 registerBits, bitsX, bitsA, maxBitsMultiplyUnitX, maxBitsMultiplyUnitY, MaxDeltaAndFactors(),
                 scaleFactor, ptSplitValue, powerTwoRanges);
      fillMatrix(matricesHighPt, index, buildFullFileName(fileNameHighPt, baseDirHighPt, index), deltaX, deltaA,
                 registerBits, bitsX, bitsA, maxBitsMultiplyUnitX, maxBitsMultiplyUnitY, MaxDeltaAndFactors(),
                 scaleFactor, ptSplitValue, powerTwoRanges);
    }
    catch (int exception) {
      std::cout << "Error: Matrix for combination = " << index << " not found" << std::endl;
      throw;
    }
  }

  // Normalize matrices so that each row is encoded with the same range (it can be different for different rows)
  // Make it so that low and high pT coefficients are normalized consistently since they will use the same slices
  // in the firmware.
  MaxDeltaAndFactors maxDeltaAndFactors(computeMaxDeltaAndFactor(matricesLowPt, matricesHighPt));
  matricesLowPt->clear();
  matricesHighPt->clear();
  for (auto index : combinationIndexList) {
    try {
      fillMatrix(matricesLowPt, index, buildFullFileName(fileNameLowPt, baseDirLowPt, index), deltaX, deltaA,
                 registerBits, bitsX, bitsA, maxBitsMultiplyUnitX, maxBitsMultiplyUnitY, maxDeltaAndFactors,
                 scaleFactor, ptSplitValue);
      fillMatrix(matricesHighPt, index, buildFullFileName(fileNameHighPt, baseDirHighPt, index), deltaX, deltaA,
                 registerBits, bitsX, bitsA, maxBitsMultiplyUnitX, maxBitsMultiplyUnitY, maxDeltaAndFactors,
                 scaleFactor, ptSplitValue);
    }
    catch (int exception) {
      std::cout << "Error: Matrix for combination = " << index << " unable to normalize" << std::endl;
      throw;
    }
  }

  bool writeMatrices = false;
  if (writeMatrices) {
    for (auto index : combinationIndexList) {
      try {
        matricesLowPt->find(index)->second.write();
        matricesHighPt->find(index)->second.write();
      }
      catch (int exception) {
        std::cout << "Error: Matrix for combination = " << index << " unable to write" << std::endl;
        throw;
      }
    }
  }
}
