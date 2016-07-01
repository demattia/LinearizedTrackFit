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
                                                     const std::string & linearFitLongitudinalDirName) :
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
    rotationFactor_(0.)
{
  if (extrapolatedRPrecision_ < 0 || extrapolatedRPrecision_ > 4) {
    std::cout << "Error: extrapolated R precision can only be a number between 0 and 4 (included). ";
    std::cout << "Extrapolated R precision requested = " << extrapolatedRPrecision_ << std::endl;
  }

  if (correctNonRadialStrips_) extrapolateR_ = true;
  if (linearFitLowPtDirName_ == "") {
    preEstimatePtDirName_ = baseDir_+"/FlatPt/PreEstimate_Transverse";
    preEstimateCotThetaDirName_ = baseDir_+"/FlatPt/PreEstimate_Longitudinal_Rz";
    linearFitLowPtDirName_ = baseDir_+"/FlatOneOverPt/Combinations_FullCorrections_2_10_SecondFirst_flatPtTgThetaAndPtPreEstimate";
    linearFitHighPtDirName_ = baseDir_+"/FlatPt/Combinations_FullCorrections_10_more_SecondFirst_flatPtTgThetaAndPtPreEstimate";
    linearFitLongitudinalDirName_ = baseDir_+"/FlatPt/Combinations_Longitudinal_Rz_flatPtTgThetaAndPtPreEstimate";
  }
}


void LinearizedTrackFitterBase::consistencyCheck(const std::vector<double> & vars, const std::vector<int> & layers)
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
  }
}


void LinearizedTrackFitterBase::fillMatrix(std::unordered_map<unsigned long, EstimatorSimple> * matrices,
                                           const unsigned long index, const std::string & fullFileName,
                                           const double & deltaPhi, const double & deltaR, const int registerBits,
                                           const int bitsPhi, const int bitsA,
                                           const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                                           const double & scaleFactor, const double & ptSplitValue,
                                           const std::vector<bool> & powerTwoRanges)
{
  matrices->insert(std::make_pair(index, EstimatorSimple(fullFileName, scaleFactor)));
}


void LinearizedTrackFitterBase::fillMatrix(std::unordered_map<unsigned long, MatrixReader> * matrices,
                                           const unsigned long index, const std::string & fullFileName,
                                           const double & deltaPhi, const double & deltaR, const int registerBits,
                                           const int bitsPhi, const int bitsA,
                                           const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                                           const double & scaleFactor, const double & ptSplitValue,
                                           const std::vector<bool> & powerTwoRanges)
{
  matrices->insert(std::make_pair(index, MatrixReader(fullFileName)));
}


std::string LinearizedTrackFitterBase::buildFullFileName(const std::string & fileName, const std::string & baseDir,
                                                         const unsigned long & index)
{
  std::string fullFileName(fileName);
  fullFileName.replace(fullFileName.find("0"), 1, std::to_string(index));
  fullFileName = baseDir + "/" + fullFileName;
  return fullFileName;
}
