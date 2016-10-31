#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearizedTrackFitter.h"

LinearizedTrackFitter::LinearizedTrackFitter(const std::string & baseDir, const bool inputExtrapolateR,
                                             const int inputExtrapolatedRPrecision,
                                             const bool inputCorrectNonRadialStrips, const int regionsNumber,
                                             const std::string & preEstimatePtDirName,
                                             const std::string & preEstimateCotThetaDirName,
                                             const std::string & linearFitLowPtDirName,
                                             const std::string & linearFitHighPtDirName,
                                             const std::string & linearFitLongitudinalDirName,
                                             const bool alignPrincipals) :
    LinearizedTrackFitterBase(baseDir, inputExtrapolateR, inputExtrapolatedRPrecision,
                              inputCorrectNonRadialStrips, regionsNumber,
                              preEstimatePtDirName, preEstimateCotThetaDirName,
                              linearFitLowPtDirName, linearFitHighPtDirName, linearFitLongitudinalDirName,
                              alignPrincipals),
    preEstimatedPt_(0.)
{
  // Fill all pre-estimates
  fillMatrices(preEstimatePtDirName_, "matrixVD_0_pre_chargeOverPt.txt", &chargeOverPtEstimator_);
  // R and z are assumed to have the same number of layers. If not the estimator needs to be modified.
  fillMatrices(preEstimateCotThetaDirName_, "matrixVD_0_pre_cotTheta.txt", &cotThetaEstimator_);
  fillMatrices(preEstimateCotThetaDirName_, "matrixVD_0_pre_tgTheta.txt", &tgThetaEstimator_);

  // Fill all PCA coefficients for parameters and chi2 estimates
  fillMatrices(linearFitLowPtDirName_, "matrixVD_0.txt", &linearFitLowPt_);
  fillMatrices(linearFitHighPtDirName_, "matrixVD_0.txt", &linearFitHighPt_);
  fillMatrices(linearFitLongitudinalDirName_, "matrixVD_0.txt", &linearFitLongitudinal_);
}


double LinearizedTrackFitter::fit(const std::vector<double> & vars, const int bits)
{
  std::vector<double> cleanedVars;
  std::vector<int> layers;
  fillVariablesFromBits(vars, bits, cleanedVars, layers);
  return fit(cleanedVars, layers);
}


double LinearizedTrackFitter::fit(const std::vector<double> & vars, const std::vector<int> & layers)
{
//  std::cout << "vars =" << std::endl;
//  for (auto it : vars) std::cout << it << std::endl;

  consistencyCheck(vars, layers);
  varsNum_ = vars.size()/3;
  varsR_.clear();
  varsR_.reserve(varsNum_);
  correctedVarsPhi_ = Matrix<long double, Dynamic, 1>(varsNum_);
  correctedVarsZ_ = Matrix<long double, Dynamic, 1>(varsNum_);

  // Compute the phi alignment value
  rotationFactor_ = computeRotationFactor(vars);

  // Rotate phi accounting for the discontinuity at +/-pi.
  // The rotation factor is computed for the innermost stub. If any of the outer stubs goes across the +/-pi edge it
  // appears as a discontinuity. Given the size of a tower it cannot happen that the innermost stub is above/below
  // +/-pi/2 (rotation factor of at least +/-1.2) and any other stub has opposite sign below the 0 since they cannot be
  // together in the same trigger tower. If rotationFactor < -1.2 then the only possible change of sign within a tower
  // is for a stub to be across the -pi border (becomes positive sign) and for rotationFactor > 1.2 the only possible
  // change of sign is for a stub across the +pi border (becomes negative sign). Therefore, we rotate those cases by
  // 2pi to make all the stubs have consistent sign and avoid the discontinuity.
  for (unsigned int i=0; i<varsNum_; ++i) {
    if (rotationFactor_ < -1.2 && vars[i*3] > 0.) correctedVarsPhi_(i) = vars[i*3]-2*M_PI - rotationFactor_;
    else if (rotationFactor_ > 1.2 && vars[i*3] < 0.) correctedVarsPhi_(i) = vars[i*3]+2*M_PI - rotationFactor_;
    else correctedVarsPhi_(i) = vars[i*3] - rotationFactor_;
  }
  for (unsigned int i=0; i<varsNum_; ++i) { varsR_.push_back(vars[i*3+1]); }
  for (unsigned int i=0; i<varsNum_; ++i) { correctedVarsZ_(i) = vars[i*3+2]; }
  extrapolatedR_ = varsR_;

  uniqueLayers_ = layers;
  std::sort(uniqueLayers_.begin(), uniqueLayers_.end());
  uniqueLayers_.erase(std::unique(uniqueLayers_.begin(), uniqueLayers_.end()), uniqueLayers_.end());
  combinationIndex_ = combinationIndex(uniqueLayers_, varsR_, regionsNumber_);

  auto iterPt = chargeOverPtEstimator_.find(combinationIndex_);
  if (iterPt == chargeOverPtEstimator_.end()) {
    return -1.;
  }

  if (!readMean(preEstimateCotThetaDirName_, "MeanRadius_", combinationIndex_, meanRadius_)) {
    std::cout << "Error: mean radii not found for combination = " << combinationIndex_ << std::endl;
    throw;
  }

  // Correct the input variables and split them between phi and z vectors
  EstimatorSimple & chargeOverPtEstimator = iterPt->second;
  double preEstimatedChargeOverPt = chargeOverPtEstimator.estimate(correctedVarsPhi_);
  auto iterCotTheta = cotThetaEstimator_.find(combinationIndex_);
  EstimatorSimple & cotThetaEstimator = iterCotTheta->second;

  preEstimatedPt_ = 1./fabs(preEstimatedChargeOverPt);
  // Retake it here because we need it with the charge
  double chargeOverTwoRho = (3.8114*0.003)*preEstimatedChargeOverPt/2.;
  double cotTheta = cotThetaEstimator.estimate(varsR_, correctedVarsZ_);

  // Extrapolate R if required
  // Warning: do not put in the following loop or the correctedVarsZ will be modified for all the elements after the
  // first one and the results will be incorrect.
  double tgTheta = 0.;
  if (extrapolateR_) {
    auto iterTgTheta = tgThetaEstimator_.find(combinationIndex_);
    EstimatorSimple &tgThetaEstimator = iterTgTheta->second;
    tgTheta = tgThetaEstimator.estimate(varsR_, correctedVarsZ_);
  }

  return fit(chargeOverTwoRho, cotTheta, tgTheta);
}


double LinearizedTrackFitter::fit(const double & chargeOverTwoRho, const double & cotTheta, const double & tgTheta)
{
  // Extrapolate R if required
  // Warning: do not put in the following loop or the correctedVarsZ will be modified for all the elements after the
  // first one and the results will be incorrect.
  if (extrapolateR_) {
    firstOrderTerm_.clear();
    secondOrderTerm1_.clear();
    secondOrderTerm2_.clear();
    secondOrderTerm3_.clear();

    // Force first order R extrapolation
    // extrapolatedRPrecision_ = 0;

    for (unsigned int i=0; i<varsNum_; ++i) {
      if (extrapolatedRPrecision_ == 0) {
        extrapolatedR_[i] = extrapolateRFirstOrder(varsR_[i], correctedVarsZ_[i], uniqueLayers_[i], tgTheta,
                                                   chargeOverTwoRho, uniqueLayers_, varsR_, correctedVarsZ_,
                                                   firstOrderTerm_);
      }
      else if (extrapolatedRPrecision_ == 1) {
        extrapolatedR_[i] = extrapolateRSecondOrderFirstTermOnly(varsR_[i], correctedVarsZ_[i], uniqueLayers_[i], tgTheta,
                                                                 chargeOverTwoRho, uniqueLayers_, varsR_, correctedVarsZ_,
                                                                 firstOrderTerm_, secondOrderTerm1_);
      }
      else if (extrapolatedRPrecision_ == 2) {
        extrapolatedR_[i] = extrapolateRSecondOrderFirstTwoTermsOnly(varsR_[i], correctedVarsZ_[i], uniqueLayers_[i], tgTheta,
                                                                     chargeOverTwoRho, uniqueLayers_, varsR_, correctedVarsZ_,
                                                                     firstOrderTerm_, secondOrderTerm1_, secondOrderTerm2_);
      }
      else if (extrapolatedRPrecision_ == 3) {
        extrapolatedR_[i] = extrapolateRSecondOrder(varsR_[i], correctedVarsZ_[i], uniqueLayers_[i], tgTheta,
                                                    chargeOverTwoRho, uniqueLayers_, varsR_, correctedVarsZ_,
                                                    firstOrderTerm_, secondOrderTerm1_, secondOrderTerm2_,
                                                    secondOrderTerm3_);
      }
      else if (extrapolatedRPrecision_ == 4) {
        extrapolatedR_[i] = extrapolateRExact(varsR_[i], correctedVarsZ_[i], uniqueLayers_[i], tgTheta,
                                              chargeOverTwoRho, uniqueLayers_, varsR_, correctedVarsZ_);
      }
      if (correctNonRadialStrips_) {
        // We correct for the rotation factor to allow for the lookup table to work properly.
        correctedVarsPhi_[i] = correctPhiForNonRadialStripsLookup_.correctPhiForNonRadialStrips(correctedVarsPhi_[i]+rotationFactor_, 0.009, extrapolatedR_[i],
                                                                                                varsR_[i], correctedVarsZ_[i], uniqueLayers_[i]) - rotationFactor_;
      }
    }
  }

  for (unsigned int i=0; i<varsNum_; ++i) {
    double DeltaR = varsR_[i] - meanRadius_[combinationIndex_][i];
    double RCube = std::pow(varsR_[i], 3);
    // Note: the extrapolatedR = R and is only extrapolated for 2S modules in the disks if requested
    double DeltaExtrapolatedR = extrapolatedR_[i] - meanRadius_[combinationIndex_][i];
    double extrapolatedRCube = std::pow(extrapolatedR_[i], 3);
    correctedVarsPhi_[i] += chargeOverTwoRho * DeltaExtrapolatedR + extrapolatedRCube * std::pow(chargeOverTwoRho, 3) / 6.;
    // We use the regular R for the R-z plane. We could recompute constants with the extrapolated R (likely negligible difference).
    correctedVarsZ_[i] -= (DeltaR + 1/6.*RCube*(chargeOverTwoRho*chargeOverTwoRho))*cotTheta;
  }

  // Evaluate the chi2/ndof
  MatrixReader * linearFitLongitudinal = &(linearFitLongitudinal_.find(combinationIndex_)->second);
  ndofLongitudinal_ = linearFitLongitudinal->nDof();
  chi2Longitudinal_ = linearFitLongitudinal->normChi2(correctedVarsZ_)*ndofLongitudinal_;
  MatrixReader * linearFitTransverse = nullptr;
  if (preEstimatedPt_ < ptSplitValue_) linearFitTransverse = &(linearFitLowPt_.find(combinationIndex_)->second);
  else linearFitTransverse = &(linearFitHighPt_.find(combinationIndex_)->second);
  ndofTransverse_ = linearFitTransverse->nDof();
  chi2Transverse_ = linearFitTransverse->normChi2(correctedVarsPhi_)*ndofTransverse_;

  // Estimate the track parameters
  estimatedPars_.clear();
  estimatedPars_ = linearFitTransverse->trackParameters(correctedVarsPhi_);
  // Parameter 1 must be phi0 for the rotation.
  if (estimatedPars_.size() > 1) estimatedPars_.at(1) += rotationFactor_;
  auto tempPars = linearFitLongitudinal->trackParameters(correctedVarsZ_);
  estimatedPars_.insert(estimatedPars_.end(), tempPars.begin(), tempPars.end());

  return (chi2Transverse_+chi2Longitudinal_)/(ndofTransverse_+ndofLongitudinal_);
}


std::vector<double> LinearizedTrackFitter::principalComponents()
{
  principalComponents_.clear();
  MatrixReader * linearFitTransverse = nullptr;
  if (preEstimatedPt_ < ptSplitValue_) linearFitTransverse = &(linearFitLowPt_.find(combinationIndex_)->second);
  else linearFitTransverse = &(linearFitHighPt_.find(combinationIndex_)->second);
  principalComponents_ = linearFitTransverse->principalComponents(correctedVarsPhi_);
  if (alignPrincipals_) alignPrincipals(principalComponents_, linearFitTransverse->nDof());
//  auto tempPrincipalFromZ = linearFitLongitudinal_.find(combinationIndex_)->second.principalComponents(correctedVarsZ_);
  MatrixReader * linearFitLongitudinal = &(linearFitLongitudinal_.find(combinationIndex_)->second);
  auto tempPrincipalFromZ = linearFitLongitudinal->principalComponents(correctedVarsZ_);
  if (alignPrincipals_) alignPrincipals(tempPrincipalFromZ, linearFitLongitudinal->nDof());
  principalComponents_.insert(principalComponents_.end(), tempPrincipalFromZ.begin(), tempPrincipalFromZ.end());
  return principalComponents_;
}


std::vector<double> LinearizedTrackFitter::normalizedPrincipalComponents()
{
  normalizedPrincipalComponents_.clear();
  MatrixReader * linearFitTransverse = nullptr;
  if (preEstimatedPt_ < ptSplitValue_) linearFitTransverse = &(linearFitLowPt_.find(combinationIndex_)->second);
  else linearFitTransverse = &(linearFitHighPt_.find(combinationIndex_)->second);
  normalizedPrincipalComponents_ = linearFitTransverse->normalizedPrincipalComponents(correctedVarsPhi_);
  if (alignPrincipals_) alignPrincipals(normalizedPrincipalComponents_, linearFitTransverse->nDof());
//  auto tempPrincipalFromZ = linearFitLongitudinal_.find(combinationIndex_)->second.normalizedPrincipalComponents(correctedVarsZ_);
  MatrixReader * linearFitLongitudinal = &(linearFitLongitudinal_.find(combinationIndex_)->second);
  auto tempPrincipalFromZ = linearFitLongitudinal->normalizedPrincipalComponents(correctedVarsZ_);
  if (alignPrincipals_) alignPrincipals(tempPrincipalFromZ, linearFitLongitudinal->nDof());
  normalizedPrincipalComponents_.insert(normalizedPrincipalComponents_.end(), tempPrincipalFromZ.begin(), tempPrincipalFromZ.end());
  return normalizedPrincipalComponents_;
}
