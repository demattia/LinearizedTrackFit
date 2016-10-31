//
// Created by Marco De Mattia on 9/2/15.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearizedTrackFitterEmulator.h"
// typedef boost::multiprecision::int128_t bigInt;
// typedef long long int bigInt;
// typedef int64_t bigInt;

LinearizedTrackFitterEmulator::LinearizedTrackFitterEmulator(const std::string & baseDir, const bool inputExtrapolateR,
                                                             const int inputExtrapolatedRPrecision,
                                                             const bool inputCorrectNonRadialStrips, const int regionsNumber,
                                                             const int registerBits, const int bitsX, const int bitsA,
                                                             const double & deltaPhi, const double & deltaR,
                                                             const double & deltaZ, const int maxBitsMultiplyUnitX,
                                                             const int maxBitsMultiplyUnitY,
                                                             const std::string & preEstimatePtDirName,
                                                             const std::string & preEstimateCotThetaDirName,
                                                             const std::string & linearFitLowPtDirName,
                                                             const std::string & linearFitHighPtDirName,
                                                             const std::string & linearFitLongitudinalDirName,
                                                             const int reducedAccuracyBitsPhi,
                                                             const int reducedAccuracyBitsR,
                                                             const int reducedAccuracyBitsZ,
                                                             const int reducedAccuracyBitsStripIndex,
                                                             const bool normalizeMatrices,
                                                             const bool commonLowHighPtFix,
                                                             const bool powerTwoRangeChargeOverPt,
                                                             const bool saveOutput,
                                                             const bool alignPrincipals) :
    LinearizedTrackFitterBase(baseDir, inputExtrapolateR, inputExtrapolatedRPrecision,
                              inputCorrectNonRadialStrips, regionsNumber,
                              preEstimatePtDirName, preEstimateCotThetaDirName, linearFitLowPtDirName,
                              linearFitHighPtDirName, linearFitLongitudinalDirName, alignPrincipals),
    registerBits_(registerBits), bitsX_(bitsX), bitsA_(bitsA),
    shiftBitsPhi_(0), shiftBitsR_(0), shiftBitsZ_(0),
    deltaPhi_(deltaPhi), deltaR_(deltaR), deltaZ_(deltaZ),
    reductionPt_(0), reductionCotTheta_(0), reductionTgTheta_(0), maxBitsMultiplyUnitX_(maxBitsMultiplyUnitX),
    maxBitsMultiplyUnitY_(maxBitsMultiplyUnitY),
    reducedAccuracyBitsStripIndex_(reducedAccuracyBitsStripIndex),
    // varOut_("var_out.txt", std::ios_base::app | std::ios_base::out),
    preChargeOverTwoRho_(0),
    saveOutput_(saveOutput)
{
  if (saveOutput_) {
    varOut_ .open("var_out.txt", std::ios_base::out);
    emulatorOutput_.open("emulator_output.txt", std::ios_base::out);
  }

  shiftBitsPhi_ = computeShiftBits(bitsX, reducedAccuracyBitsPhi);
  shiftBitsR_ = computeShiftBits(bitsX, reducedAccuracyBitsR);
  shiftBitsZ_ = computeShiftBits(bitsX, reducedAccuracyBitsZ);
  // Evaluate the reduction factors.
  reductionPhi_ = 1;
  while (std::pow(2., reductionPhi_) < deltaPhi) {
    ++reductionPhi_;
  }
  if (deltaPhi - int(std::pow(2, reductionPhi_)) != 0) {
    std::cout << "Error: deltaPhi must be a power of 2. The value is instead = " << deltaPhi << std::endl;
    throw;
  }
  reductionZ_ = 1;
  while (std::pow(2., reductionZ_) < deltaZ) {
    ++reductionZ_;
  }
  if (deltaZ - int(std::pow(2, reductionZ_)) != 0) {
    std::cout << "Error: deltaZ must be a power of 2. The value is instead = " << deltaPhi << std::endl;
    throw;
  }

  // Store the encoded value of the PS-2S split in the disks
  Rcut_ = encode(61., deltaR_, bitsX_);
  // Store the encoded values of the 2S split in the disks
  RadiusCut1_ = encode(82.5, deltaR, bitsX);
  RadiusCut2_ = encode(77., deltaR, bitsX);

  // Fill all pre-estimates
  fillMatrices(preEstimatePtDirName_, "matrixVD_0_pre_chargeOverPt.txt", &chargeOverPtEstimator_, deltaPhi_, deltaR_,
               registerBits_, bitsX_, bitsA_, maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_,
               (3.8114*0.003)/2., ptSplitValue_, normalizeMatrices);

  // R and z are assumed to have the same number of layers. If not the estimator needs to be modified.
  // We change the sign of the cot(theta) pre-estimate coefficients to match the firmware where we estimate
  // -cot(theta) to obtain the correction terms with the right sign and perform only additions.
  fillMatrices(preEstimateCotThetaDirName_, "matrixVD_0_pre_cotTheta.txt", &cotThetaEstimator_, deltaZ_, deltaR_,
               registerBits_, bitsX_, bitsA_, maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_, -1., 10., normalizeMatrices);
  // Exclude the barrel for the normalization of the ranges for tgTheta only since tgTheta is not used in the barrel.
  bool excludeBarrelFromNormalizedMatrices = true;
  fillMatrices(preEstimateCotThetaDirName_, "matrixVD_0_pre_tgTheta.txt", &tgThetaEstimator_, deltaZ_, deltaR_,
               registerBits_, bitsX_, bitsA_, maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_, 1., 10., normalizeMatrices,
               excludeBarrelFromNormalizedMatrices);

  std::vector<bool> powerTwoRanges(6, true);
  // This requires the c/pT to be the first estimated parameter
  powerTwoRanges[0] = powerTwoRangeChargeOverPt;

  // Fill all PCA coefficients for parameters and chi2 estimates
  if (commonLowHighPtFix) {
    // Normalize low pT and high pT coefficients consistently since they will use the same slices in the firmware.
    fillMatrices(linearFitLowPtDirName_, "matrixVD_0.txt", &linearFitLowPt_,
                 linearFitHighPtDirName_, "matrixVD_0.txt", &linearFitHighPt_,
                 deltaPhi_, deltaR_, registerBits_, bitsX_, bitsA_, maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_,
                 1., 10., powerTwoRanges);
  }
  else {
    fillMatrices(linearFitLowPtDirName_, "matrixVD_0.txt", &linearFitLowPt_, deltaPhi_, deltaR_,
                 registerBits_, bitsX_, bitsA_, maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_, 1., 10.,
                 normalizeMatrices, false, powerTwoRanges);
    fillMatrices(linearFitHighPtDirName_, "matrixVD_0.txt", &linearFitHighPt_, deltaPhi_, deltaR_,
                 registerBits_, bitsX_, bitsA_, maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_, 1., 10.,
                 normalizeMatrices, false, powerTwoRanges);
  }

  fillMatrices(linearFitLongitudinalDirName_, "matrixVD_0.txt", &linearFitLongitudinal_, deltaZ_, deltaR_,
               registerBits_, bitsX_, bitsA_, maxBitsMultiplyUnitX, maxBitsMultiplyUnitY_, 1., 10., normalizeMatrices);
}


double LinearizedTrackFitterEmulator::fit(const std::vector<double> & vars, const std::vector<int> & layers,
                                          const std::vector<int> & stripIndexes)
{
  initialize(vars, layers, stripIndexes);

  // Perform floating point quantization (encoding)

  rotationFactor_ = computeRotationFactor(vars);

  for (unsigned int i=0; i<varsNum_; ++i) {
    bigInt tempVarPhiInt = encode(vars[i*3]-rotationFactor_, deltaPhi_, bitsX_);

    // Rotate phi accounting for the discontinuity at +/-pi. See LinearizedTrackFitter.cc for more details.
    if (rotationFactor_ < -1.2 && vars[i*3] > 0.) tempVarPhiInt = encode(vars[i*3]-2*M_PI - rotationFactor_, deltaPhi_, bitsX_);
    else if (rotationFactor_ > 1.2 && vars[i*3] < 0.) tempVarPhiInt = encode(vars[i*3]+2*M_PI - rotationFactor_, deltaPhi_, bitsX_);

    if (shiftBitsPhi_ != 0) {
      truncate(tempVarPhiInt, shiftBitsPhi_);
      alignBits(tempVarPhiInt, -shiftBitsPhi_);
    }
    varsPhiInt_.push_back(tempVarPhiInt);
  }
  for (unsigned int i=0; i<varsNum_; ++i) {
    bigInt tempVarRInt = encode(vars[i*3+1], deltaR_, bitsX_);
    if (shiftBitsR_ != 0) {
      truncate(tempVarRInt, shiftBitsR_);
      alignBits(tempVarRInt, -shiftBitsR_);
    }
    varsRInt_.push_back(tempVarRInt);
    // oneOverRInt_.push_back(encode(1./vars[i*3+1], deltaR_, bitsX_));
    // double value = 0.009/(vars[i*3+1]*vars[i*3+1]);
    // Encode with a range of -8 so that there it fits in 18 bits after the product with DeltaStripIndex.
    // oneOverRSquaredInt_.push_back(encode(0.009/(vars[i*3+1]*vars[i*3+1]), std::pow(2, -8), 18));
  }
  for (unsigned int i=0; i<varsNum_; ++i) {
    bigInt tempVarZInt = encode(vars[i*3+2], deltaZ_, bitsX_);
    if (shiftBitsZ_ != 0) {
      truncate(tempVarZInt, shiftBitsZ_);
      alignBits(tempVarZInt, -shiftBitsZ_);
    }
    varsZInt_.push_back(tempVarZInt);
  }
  extrapolatedRInt_ = varsRInt_;

  if (stripIndexes.size() > 6) {
    std::cout << "Error: strip indexes vector size is greater than 6. Does it contain duplicates?" << std::endl;
    throw;
  }
  stripIndex_.clear();
  for (auto index : stripIndexes) {
    // Reduce accuracy to the required number of bits and add zeros to ensure that the result is still a number
    // in the same range as the original (0 to 1015).
    bigInt bigIntStripIndex(index);
    if (reducedAccuracyBitsStripIndex_ != 0) {
      truncate(bigIntStripIndex, reducedAccuracyBitsStripIndex_);
      alignBits(bigIntStripIndex, -reducedAccuracyBitsStripIndex_);
    }
    // std::cout << bigIntStripIndex << std::endl;
    // stripIndex_.push_back(index);
    stripIndex_.push_back(bigIntStripIndex);
  }

  // Keep floating point R for computing the combination index. This is only for convenience since utilizing
  // the encoded value would require to change (encode) also the cut values. The distance among strip centers
  // in the disks is about 5 cm. This is much larger than the accuracy of the encoded R except for very low
  // number of bits (that are not realistic).
  std::vector<double> varsR;
  for (unsigned int i=0; i<varsNum_; ++i) { varsR.push_back(vars[i*3+1]); }
  combinationIndex_ = combinationIndex(uniqueLayers_, varsR, regionsNumber_);
  return fit();
}


double LinearizedTrackFitterEmulator::fit(const std::vector<bigInt> & varsInt, const std::vector<int> & layers,
                                          const std::vector<int> & stripIndexes)
{
  initialize(varsInt, layers, stripIndexes);
  // Extract the variables from the vector
  for (unsigned int i=0; i<varsNum_; ++i) {
    bigInt tempPhiInt = varsInt.at(i*3);
    if (shiftBitsPhi_ != 0) alignBits(tempPhiInt, -shiftBitsPhi_);
    varsPhiInt_.push_back(tempPhiInt);

    bigInt tempRInt = varsInt.at(i*3+1);
    if (shiftBitsR_ != 0) alignBits(tempRInt, -shiftBitsR_);
    varsRInt_.push_back(tempRInt);

    bigInt tempZInt = varsInt.at(i*3+2);
    if (shiftBitsZ_ != 0) alignBits(tempZInt, -shiftBitsZ_);
    varsZInt_.push_back(tempZInt);
  }
  extrapolatedRInt_ = varsRInt_;

  stripIndex_.clear();
  for (auto index : stripIndexes) stripIndex_.push_back(bigInt(index));

  combinationIndex_ = combinationIndex(uniqueLayers_, varsRInt_, regionsNumber_, Rcut_, RadiusCut1_, RadiusCut2_);
  return fit();
}


double LinearizedTrackFitterEmulator::fit()
{
  // std::cout << "combination index = " << combinationIndex_ << std::endl;
  if (saveOutput_) {
    varOut_ << "combination index = " << combinationIndex_ << std::endl;
    emulatorOutput_ << "combination index = " << combinationIndex_ << std::endl;

    // std::cout << "varsPhiInt_ = " << std::endl;
    varOut_ << "varsPhiInt_ = " << std::endl;
    emulatorOutput_ << "varsPhiInt_ = " << std::endl;
    for (const auto &phiInt : varsPhiInt_) {
      // std::cout << phiInt << std::endl;
      varOut_ << phiInt << std::endl;
      emulatorOutput_ << phiInt << std::endl;
    }
    // std::cout << "varsRInt_ = " << std::endl;
    varOut_ << "varsRInt_ = " << std::endl;
    emulatorOutput_ << "varsRInt_ = " << std::endl;
    for (const auto &RInt : varsRInt_) {
      // std::cout << RInt << std::endl;
      varOut_ << RInt << std::endl;
      emulatorOutput_ << RInt << std::endl;
    }
    // std::cout << "varsZInt_ = " << std::endl;
    varOut_ << "varsZInt_ = " << std::endl;
    emulatorOutput_ << "varsZInt_ = " << std::endl;
    for (const auto &zInt : varsZInt_) {
      // std::cout << zInt << std::endl;
      varOut_ << zInt << std::endl;
      emulatorOutput_ << zInt << std::endl;
    }
  }

  auto iterPt = chargeOverPtEstimator_.find(combinationIndex_);
  if (iterPt == chargeOverPtEstimator_.end()) {
    return -1.;
  }

  if (!readMean(preEstimateCotThetaDirName_, "MeanRadius_", combinationIndex_, meanRadius_)) {
    std::cout << "Error: mean radii not found for combination = " << combinationIndex_ << std::endl;
    throw;
  }
  else {
    std::vector<bigInt> meanRadiusIntVector_;
    encodeVector(meanRadius_[combinationIndex_], meanRadiusIntVector_, deltaR_, bitsX_);
    meanRadiusInt_.insert(std::make_pair(combinationIndex_, meanRadiusIntVector_));
  }

  // Correct the input variables and split them between phi and z vectors
  EstimatorSimpleEmulator & chargeOverTwoRhoEstimator = iterPt->second;
  bigInt preChargeOverTwoRho = chargeOverTwoRhoEstimator.estimate(varsPhiInt_);

  reductionPt_ = chargeOverTwoRhoEstimator.reduction();
  oneOverTwoRhoSplitValueInt_ = chargeOverTwoRhoEstimator.chargeOverTwoRhoSplitValue();

  auto iterCotTheta = cotThetaEstimator_.find(combinationIndex_);
  EstimatorSimpleEmulator & cotThetaEstimator = iterCotTheta->second;

  // Retake it here because we need it with the charge
  bigInt cotTheta = cotThetaEstimator.estimate(varsRInt_, varsZInt_);
  reductionCotTheta_ = cotThetaEstimator.reduction();
  // std::cout << "decoded cotTheta = " << decode(cotTheta, deltaZ_*cotThetaEstimator.deltaA(), bitsX_) << std::setprecision(20) << std::endl;

  // Extrapolate R if required
  // Warning: do not put in the following loop or the correctedVarsZ will be modified for all the elements after the
  // first one and the results will be incorrect.
  bigInt tgTheta = 0;
  if (extrapolateR_) {
    auto iterTgTheta = tgThetaEstimator_.find(combinationIndex_);
    EstimatorSimpleEmulator &tgThetaEstimator = iterTgTheta->second;
    reductionTgTheta_ = tgThetaEstimator.reduction();
    tgTheta = tgThetaEstimator.estimate(varsRInt_, varsZInt_);
    // std::cout << "tgThetaEstimator.deltaA() = " << tgThetaEstimator.deltaA() << std::endl;
    // std::cout << "decoded tgTheta = " << decode(tgTheta, deltaZ_*tgThetaEstimator.deltaA(), bitsX_) << std::setprecision(20) << std::endl;
  }

  preChargeOverTwoRho_ = preChargeOverTwoRho;
  return fit(preChargeOverTwoRho, cotTheta, tgTheta);
}


double LinearizedTrackFitterEmulator::fit(const bigInt & chargeOverTwoRho, const bigInt & cotTheta, const bigInt & tgTheta)
{
  // Extrapolate R if required
  // Warning: do not put in the following loop or the correctedVarsZ will be modified for all the elements after the
  // first one and the results will be incorrect.
  if (extrapolateR_) {
    for (unsigned int i=0; i<varsNum_; ++i) {
      if (extrapolatedRPrecision_ == 0) {
        extrapolatedRInt_[i] = extrapolateRFirstOrderEmulator(varsRInt_[i], varsZInt_[i], uniqueLayers_[i],
                                                              tgTheta, chargeOverTwoRho, uniqueLayers_, varsRInt_,
                                                              varsZInt_, Rcut_, registerBits_, bitsX_, deltaZ_,
                                                              reductionPhi_, reductionZ_,
                                                              reductionPt_, reductionTgTheta_,
                                                              maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_);
      }
      else if (extrapolatedRPrecision_ == 1) {
        extrapolatedRInt_[i] = extrapolateRSecondOrderFirstTermOnlyEmulator(varsRInt_[i], varsZInt_[i], uniqueLayers_[i],
                                                                            tgTheta, chargeOverTwoRho, uniqueLayers_, varsRInt_,
                                                                            varsZInt_, Rcut_, registerBits_, bitsX_, deltaZ_,
                                                                            reductionPhi_, reductionZ_,
                                                                            reductionPt_, reductionTgTheta_,
                                                                            maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_);
      }
      else if (extrapolatedRPrecision_ == 2) {
        extrapolatedRInt_[i] = extrapolateRSecondOrderFirstTwoTermsOnlyEmulator(varsRInt_[i], varsZInt_[i], uniqueLayers_[i],
                                                                                tgTheta, chargeOverTwoRho, uniqueLayers_, varsRInt_,
                                                                                varsZInt_, Rcut_, registerBits_, bitsX_, deltaZ_,
                                                                                reductionPhi_, reductionZ_,
                                                                                reductionPt_, reductionTgTheta_,
                                                                                maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_);
      }
      else if (extrapolatedRPrecision_ == 3 || extrapolatedRPrecision_ == 4) {
        extrapolatedRInt_[i] = extrapolateRSecondOrderEmulator(varsRInt_[i], varsZInt_[i], uniqueLayers_[i],
                                                               tgTheta, chargeOverTwoRho, uniqueLayers_, varsRInt_,
                                                               varsZInt_, Rcut_, registerBits_, bitsX_, deltaZ_,
                                                               reductionPhi_, reductionZ_,
                                                               reductionPt_, reductionTgTheta_,
                                                               maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_);
      }
      if (correctNonRadialStrips_) {
        stripCorrectionTermsInt_[i] = correctPhiForNonRadialStripsEmulator_.correctPhiForNonRadialStrips(varsPhiInt_[i], stripIndex_[i],
                                                                                                         extrapolatedRInt_[i], varsRInt_[i],
                                                                                                         Rcut_, // oneOverRInt_[i], oneOverRSquaredInt_[i],
                                                                                                         uniqueLayers_[i], registerBits_, bitsX_,
                                                                                                         maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_,
                                                                                                         reductionZ_, reductionPhi_);
//        varsPhiInt_[i] = correctPhiForNonRadialStripsEmulator_.correctPhiForNonRadialStrips(varsPhiInt_[i], stripIndex_[i],
//                                                                                            extrapolatedRInt_[i], varsRInt_[i],
//                                                                                            Rcut_, oneOverRInt_[i], oneOverRSquaredInt_[i],
//                                                                                            uniqueLayers_[i], registerBits_, bitsX_,
//                                                                                            maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_,
//                                                                                            reductionZ_, reductionPhi_);
      }
    }
  }

  std::vector<bigInt> deltaRInt;
  preDiffer(varsRInt_, meanRadiusInt_[combinationIndex_], deltaRInt, registerBits_);
  std::vector<bigInt> deltaExtrapolatedRInt;
  preDiffer(extrapolatedRInt_, meanRadiusInt_[combinationIndex_], deltaExtrapolatedRInt, registerBits_);

  if (saveOutput_) {
    emulatorOutput_ << "meanRadiusInt_ = " << std::endl;
    for (auto meanR : meanRadiusInt_[combinationIndex_]) emulatorOutput_ << meanR << std::endl;
    emulatorOutput_ << "pre cOverTwoRho = " << chargeOverTwoRho << std::endl;
    emulatorOutput_ << "pre cotTheta = " << cotTheta << std::endl;
    emulatorOutput_ << "pre tgTheta = " << tgTheta << std::endl;
  }

  correctedVarsPhiInt_.clear();
  correctedVarsPhiInt_.reserve(varsNum_);
  phiCorrectionUnit(varsPhiInt_, deltaPhi_, bitsX_, extrapolatedRInt_, deltaExtrapolatedRInt, chargeOverTwoRho, registerBits_,
                    reductionPt_, reductionPhi_, reductionZ_, maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_,
                    stripCorrectionTermsInt_, correctedVarsPhiInt_);
//  std::cout << "decoded correctedVarsPhiInt_ = " << std::endl;
//  for (auto v : correctedVarsPhiInt_) std::cout << decode(v, deltaPhi_, bitsX_) << std::endl;
//  std::cout << std::endl;

  correctedVarsZInt_.clear();
  correctedVarsZInt_.reserve(varsNum_);
//  zCorrectionUnit(varsZInt_, deltaZ_, bitsX_, varsRInt_, deltaRInt, chargeOverTwoRho, cotTheta,
//                  registerBits_, reductionPt_, reductionCotTheta_, reductionPhi_,
//                  maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_, correctedVarsZInt_);
  zCorrectionUnit(varsZInt_, deltaZ_, bitsX_, extrapolatedRInt_, deltaExtrapolatedRInt, chargeOverTwoRho, cotTheta,
                  registerBits_, reductionPt_, reductionCotTheta_, reductionPhi_,
                  maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_, correctedVarsZInt_);
  if (saveOutput_) {
    emulatorOutput_ << "correctedVarsPhiInt_ = " << std::endl;
    for (auto v : correctedVarsPhiInt_) emulatorOutput_ << v << std::endl;
    emulatorOutput_ << "correctedVarsZInt_ = " << std::endl;
    for (auto v : correctedVarsZInt_) emulatorOutput_ << v << std::endl;
  }
//  for (auto v : correctedVarsZInt_) std::cout << v << std::endl;
//  std::cout << std::endl;

  // Evaluate the chi2/ndof
  MatrixReaderEmulator * linearFitLongitudinal = &(linearFitLongitudinal_.find(combinationIndex_)->second);
  ndofLongitudinal_ = linearFitLongitudinal->nDof();

  chi2Longitudinal_ = linearFitLongitudinal->normChi2(correctedVarsZInt_, deltaZ_, registerBits_, bitsX_, bitsA_, emulatorOutput_, saveOutput_)*ndofLongitudinal_;

  MatrixReaderEmulator * linearFitTransverse = nullptr;
  if (saveOutput_) {
    emulatorOutput_ << "chargeOverTwoRho = " << chargeOverTwoRho << std::endl;
    emulatorOutput_ << "oneOverTwoRhoSplitValueInt_ = " << oneOverTwoRhoSplitValueInt_ << std::endl;
  }
  // The comparison with the positive edge is >= to match the behavior of the firmware

  // Matches the behavior of the old version using the XOR of the sign of the two differences
  // if ((chargeOverTwoRho < -oneOverTwoRhoSplitValueInt_) || (chargeOverTwoRho >= oneOverTwoRhoSplitValueInt_)) {
  //   linearFitTransverse = &(linearFitLowPt_.find(combinationIndex_)->second);
  // }
  // New version utilizing a single DSP
  if ((chargeOverTwoRho < -oneOverTwoRhoSplitValueInt_) || (chargeOverTwoRho > oneOverTwoRhoSplitValueInt_)) {
    linearFitTransverse = &(linearFitLowPt_.find(combinationIndex_)->second);
  }
  else linearFitTransverse = &(linearFitHighPt_.find(combinationIndex_)->second);
  ndofTransverse_ = linearFitTransverse->nDof();

  chi2Transverse_ = linearFitTransverse->normChi2(correctedVarsPhiInt_, deltaPhi_, registerBits_, bitsX_, bitsA_, emulatorOutput_, saveOutput_)*ndofTransverse_;
  // Must be called after the normChi2 call.
  // chi2TermsInt_.clear();
  chi2TermsInt_ = linearFitTransverse->getChiTermsInt();
  if (alignPrincipals_) alignPrincipals(chi2TermsInt_, linearFitTransverse->nDof());
  auto tempChi2TermsInt = linearFitLongitudinal->getChiTermsInt();
  if (alignPrincipals_) alignPrincipals(tempChi2TermsInt, linearFitLongitudinal->nDof());
  chi2TermsInt_.insert(chi2TermsInt_.end(), tempChi2TermsInt.begin(), tempChi2TermsInt.end());

  // Estimate the track parameters
  // estimatedPars_.clear();
  // estimatedParsInt_.clear();

  if (saveOutput_) emulatorOutput_ << "c/pT and phi0" << std::endl;
  estimatedPars_ = linearFitTransverse->trackParameters(correctedVarsPhiInt_, deltaPhi_, registerBits_, bitsX_, bitsA_, emulatorOutput_, saveOutput_);
  // Must be called after the trackParameters call.
  estimatedParsInt_ = linearFitTransverse->getParametersInt();
  // Parameter 1 must be phi0 for the rotation.
  if (estimatedPars_.size() > 1) estimatedPars_.at(1) += rotationFactor_;
  if (saveOutput_) emulatorOutput_ << "cot(theta) and z0" << std::endl;
  auto tempPars = linearFitLongitudinal->trackParameters(correctedVarsZInt_, deltaZ_, registerBits_, bitsX_, bitsA_, emulatorOutput_, saveOutput_);
//  estimatedPars_ = linearFitTransverse->trackParameters(varsPhiInt_, deltaPhi_, registerBits_, bitsX_, bitsA_);
//  auto tempPars = linearFitLongitudinal->trackParameters(varsZInt_, deltaZ_, registerBits_, bitsX_, bitsA_);
  estimatedPars_.insert(estimatedPars_.end(), tempPars.begin(), tempPars.end());
  // Must be called after the trackParameters call.
  auto tempParsInt = linearFitLongitudinal->getParametersInt();
  estimatedParsInt_.insert(estimatedParsInt_.end(), tempParsInt.begin(), tempParsInt.end());

//  std::cout << "c/pT = " << estimatedPars_[0] << std::setprecision(12) << std::endl;
//  std::cout << "parameters = " << std::endl;
//  for (auto v : estimatedPars_) std::cout << v << std::setprecision(12) << std::endl;


//  // This is for printing the deltaA values for all cases
//  // ----------------------------------------------------
//
//  std::vector<double> minDeltaD(2, 9999999.);
//  std::vector<double> minFactorD(2, 9999999.);
//  std::vector<double> maxDeltaD(2, 0.);
//  std::vector<double> maxFactorD(2, 0.);
//
//  std::vector<double> minDeltaV(6, 9999999.);
//  std::vector<double> minFactorV(6, 9999999.);
//  std::vector<double> maxDeltaV(6, 0.);
//  std::vector<double> maxFactorV(6, 0.);
//
//  // for (auto it = linearFitLongitudinal_.begin(); it != linearFitLongitudinal_.end(); ++it) {
//  for (auto it = linearFitLowPt_.begin(); it != linearFitLowPt_.end(); ++it) {
//  // for (auto it = linearFitHighPt_.begin(); it != linearFitHighPt_.end(); ++it) {
//    size_t vSize = it->second.getDeltaVSize();
//    int shift = 0;
//    if (vSize == 5) shift = 1;
//    for (size_t i=0; i<vSize; ++i) {
//      if(minDeltaV[i+shift] > it->second.getDeltaV(i)) {
//        minDeltaV[i+shift] = it->second.getDeltaV(i);
//        minFactorV[i+shift] = it->second.getFactorV(i);
//      }
//      if(maxDeltaV[i+shift] < it->second.getDeltaV(i)) {
//        maxDeltaV[i+shift] = it->second.getDeltaV(i);
//        maxFactorV[i+shift] = it->second.getFactorV(i);
//      }
//    }
////    for (size_t i=0; i<it->second.getDeltaDSize(); ++i) {
//    for (size_t i=0; i<2; ++i) {
//      if(minDeltaD[i] > it->second.getDeltaD(i)) {
//        minDeltaD[i] = it->second.getDeltaD(i);
//        minFactorD[i] = it->second.getFactorD(i);
//      }
//      if(maxDeltaD[i] < it->second.getDeltaD(i)) {
//        maxDeltaD[i] = it->second.getDeltaD(i);
//        maxFactorD[i] = it->second.getFactorD(i);
//      }
//    }
//  }
//  for (int i=0; i<2; ++i) {
//    std::cout << "Maximum variation of bit-shift factor for parameter "<< i << " = " <<
//        maxFactorD[i] - minFactorD[i] << std::endl;
//    std::cout << "minDeltaD["<<i<<"] = " << minDeltaD[i] << std::endl;
//    std::cout << "maxDeltaD["<<i<<"] = " << maxDeltaD[i] << std::endl;
//    std::cout << "minFactorD["<<i<<"] = " << minFactorD[i] << std::endl;
//    std::cout << "maxFactorD["<<i<<"] = " << maxFactorD[i] << std::endl;
//  }
//  for (int i=0; i<6; ++i) {
//    std::cout << "Maximum variation of bit-shift factor for chi2 term "<< i << " = " <<
//    maxFactorV[i] - minFactorV[i] << std::endl;
//    std::cout << "minDeltaV["<<i<<"] = " << minDeltaV[i] << std::endl;
//    std::cout << "maxDeltaV["<<i<<"] = " << maxDeltaV[i] << std::endl;
//    std::cout << "minFactorV["<<i<<"] = " << minFactorV[i] << std::endl;
//    std::cout << "maxFactorV["<<i<<"] = " << maxFactorV[i] << std::endl;
//  }
//
//  double minDeltaA = 9999999.;
//  double maxDeltaA = 0.;
//  int minFactorA = 9999999;
//  int maxFactorA = 0;
//  // for (auto it = cotThetaEstimator_.begin(); it != cotThetaEstimator_.end(); ++it) {
//  for (auto it = chargeOverPtEstimator_.begin(); it != chargeOverPtEstimator_.end(); ++it) {
//    if(minDeltaA > it->second.deltaA()) {
//      minDeltaA = it->second.deltaA();
//      minFactorA = it->second.reduction();
//    }
//    if(maxDeltaA < it->second.deltaA()) {
//      maxDeltaA = it->second.deltaA();
//      maxFactorA = it->second.reduction();
//    }
//  }
//  std::cout << "minDeltaA = " << minDeltaA << std::endl;
//  std::cout << "maxDeltaA = " << maxDeltaA << std::endl;
//  std::cout << "minFactorA = " << minFactorA << std::endl;
//  std::cout << "maxFactorA = " << maxFactorA << std::endl;
//  std::cout << "Maximum variation of bit-shift factor for pre-estimate = " << maxFactorA - minFactorA << std::endl;


  return (chi2Transverse_+chi2Longitudinal_)/(ndofTransverse_+ndofLongitudinal_);
}


std::vector<double> LinearizedTrackFitterEmulator::normalizedPrincipalComponents()
{
  normalizedPrincipalComponents_.clear();
  MatrixReaderEmulator * linearFitTransverse = nullptr;
  if ((preChargeOverTwoRho_ < -oneOverTwoRhoSplitValueInt_) || (preChargeOverTwoRho_ >= oneOverTwoRhoSplitValueInt_)) {
    linearFitTransverse = &(linearFitLowPt_.find(combinationIndex_)->second);
  }
  else linearFitTransverse = &(linearFitHighPt_.find(combinationIndex_)->second);

  normalizedPrincipalComponents_ = linearFitTransverse->normalizedPrincipalComponents(correctedVarsPhiInt_, deltaPhi_,
                                                                                      registerBits_, bitsX_, bitsA_,
                                                                                      emulatorOutput_, saveOutput_);
  if (alignPrincipals_) alignPrincipals(normalizedPrincipalComponents_, linearFitTransverse->nDof());

//  auto tempPrincipalFromZ = linearFitLongitudinal_.find(combinationIndex_)->second.normalizedPrincipalComponents(correctedVarsZInt_, deltaZ_,
//                                                                                                                 registerBits_, bitsX_, bitsA_,
//                                                                                                                 emulatorOutput_, saveOutput_);
  MatrixReaderEmulator * linearFitLongitudinal = &(linearFitLongitudinal_.find(combinationIndex_)->second);
  auto tempPrincipalFromZ = linearFitLongitudinal->normalizedPrincipalComponents(correctedVarsZInt_, deltaZ_,
                                                                                 registerBits_, bitsX_, bitsA_,
                                                                                 emulatorOutput_, saveOutput_);
  if (alignPrincipals_) alignPrincipals(tempPrincipalFromZ, linearFitLongitudinal->nDof());
  normalizedPrincipalComponents_.insert(normalizedPrincipalComponents_.end(), tempPrincipalFromZ.begin(), tempPrincipalFromZ.end());
  return normalizedPrincipalComponents_;
}

