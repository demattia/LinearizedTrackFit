//
// Created by Marco De Mattia on 9/7/15.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/EmulatorTools.h"


bigInt encode(const double & x, const double & deltaX, const int bits)
{
  if (x >= deltaX/2.) {
    std::cout << "exceeding upper edge for x = " << x << " with upper edge = " << deltaX/2. << std::endl;
//    throw;
    return ((bigInt(1)<<(bits-1))-1);
  }
  if (x < -deltaX/2.) {
    std::cout << "exceeding lower edge for x = " << x << " with lower edge = " << -deltaX/2. << std::endl;
//    throw;
    return -(bigInt(1)<<(bits-1));
  }
  double temp = x;
  for (int i=0; i<bits; ++i) {
    temp*=2;
  }
  return bigInt(floor(temp/deltaX + 0.5));
}

template <>
void encodeVector(const std::vector<double> & vec, std::vector<bigInt> & vecInt, const double & deltaX, const int bits)
{
  for (size_t i=0; i<vec.size(); ++i) {
    vecInt.push_back(encode(vec.at(i), deltaX, bits));
  }
}


double decode(const bigInt & xInt, const double & deltaX, const int bits)
{
  // double temp = xInt.convert_to<double>();
  double temp = double(xInt);
  // std::cout << "temp = " << temp << std::endl;
  for (int i=0; i<bits; ++i) {
    temp/=2.;
  }
  return temp*deltaX;
}


void truncate(bigInt & value, const int truncateBits)
{
  // // Use this with boost::multiprecision
  // bool minusOne = false;
  // if (value < 0 && (value % bigInt(std::pow(2, truncateBits)) != 0)) minusOne = true;
  // value >>= truncateBits;
  // // If there is a rounding of a negative number, add -1 to emulate the rounding by truncation.
  // if (minusOne) value -= 1;

  // Use this with int64_t
  value >>= truncateBits;
}


void checkSizeConsistency(const std::vector<bigInt> & a, const std::vector<bigInt> & b)
{
  if (a.size() != b.size()) {
    std::cout << "Inconsistent array lengths" << std::endl;
    throw;
  }
}


void checkBufferOverflow(const bigInt & value, const int bits, const std::string & errorMessage)
{
  if ((value > ((bigInt(1)<<(bits-1))-1)) || (value < -(bigInt(1)<<(bits-1)))) {
    std::cout << errorMessage;
    std::cout << "value = " << value << std::endl;
    std::cout << "min = " << -(bigInt(1)<<(bits-1)) << std::endl;
    std::cout << "max = " << ((bigInt(1)<<(bits-1))-1) << std::endl;
    throw;
  }
}


void preDiffer(const std::vector<bigInt> & a, const std::vector<bigInt> & b, std::vector<bigInt> & c, const int bits)
{
  checkSizeConsistency(a, b);
  for (size_t i=0; i<a.size(); ++i) {
    c.push_back(a[i]-b[i]);
  }
  if (bits > 0) {
    for (auto cc : c) {
      checkBufferOverflow(cc, bits);
    }
  }
}


void alignBits(bigInt & x, const int bitShift)
{
  if (bitShift < 0) x = (x << (-bitShift));
  else truncate(x, bitShift);
}


bigInt multiplyUnit(const bigInt & x, const bigInt & y, const int registerBits, const int bitsY,
                    const int reduction, const int maxBitsX, const int maxBitsY, const bool doBitShift)
{
  bigInt tempX = x;
  bigInt tempY = y;
  if (bitsY > maxBitsY) {
    // Would be too big to fit, reduce it before the product and take this
    // additional scale into account when scaling the end result.
    truncate(tempY, bitsY-maxBitsY);
  }
  checkBufferOverflow(tempX, maxBitsX, "multiply unit overflow x");
  checkBufferOverflow(tempY, maxBitsY, "multiply unit overflow y");
  bigInt xy = tempX*tempY;
  checkBufferOverflow(xy, registerBits, "multiply unit overflow");
  if (doBitShift) {
    int bitShift = std::min(bitsY, maxBitsY) - reduction;
    alignBits(xy, bitShift);
  }
  return xy;
}


bigInt macUnit(const std::vector<bigInt> & coordinates, const std::vector<bigInt> & coefficients,
               const bigInt & parMean, const int registerBits, const int bitsCoeff, const int bitsMultiplyUnitX,
               const int bitsMultiplyUnitY)
{
  checkSizeConsistency(coordinates, coefficients);
  bigInt s = parMean;
  for (size_t i=0; i<coordinates.size(); ++i) {
    bigInt product = multiplyUnit(coordinates[i], coefficients[i], registerBits, bitsCoeff,
                                  0, bitsMultiplyUnitX, bitsMultiplyUnitY, false);
    s += product;
    checkBufferOverflow(s, registerBits, "accumulate overflow");
  }
  if (bitsCoeff != 0) {
    truncate(s, bitsCoeff);
  }
  checkBufferOverflow(s, registerBits, "mean addition overflow");
  return s;
}


bigInt adderUnit(const bigInt & x, const bigInt & y, const int bits)
{
  bigInt s = x+y;
  checkBufferOverflow(s, bits, "adder unit overflow");
  return s;
}


std::pair<int, double> computeFactor(const std::vector<double> & coefficients, const double & deltaX,
                                     const double & scaleFactor, const bool powerTwoRange)
{
  auto minA = std::min_element(coefficients.begin(), coefficients.end());
  if (minA == coefficients.end()) {
    std::cout << "computeFactor error: minimum not found." << std::endl;
    throw;
  }
  auto maxA = std::max_element(coefficients.begin(), coefficients.end());
  if (maxA == coefficients.end()) {
    std::cout << "computeFactor error: maximum not found." << std::endl;
    throw;
  }
  double deltaA = 2.*scaleFactor*std::max(-(*minA), *maxA);
  int factor = 1;
  while ((std::pow(2, factor))/deltaX < deltaA) {
    factor += 1;
  }
  if (powerTwoRange) return std::make_pair(factor, (std::pow(2, factor))/deltaX);
  return std::make_pair(factor, deltaA);
}


void phiCorrectionUnit(const std::vector<bigInt> & phiInt, const double & deltaPhi, const int bitsPhi,
                       const std::vector<bigInt> & radiusInt, std::vector<bigInt> & deltaRadiusInt,
                       const bigInt & parInt, const int registerBits, const int reductionPt, const int reductionPhi,
                       const int reductionZ, const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                       const std::vector<bigInt> & stripCorrectionTermsInt,
                       std::vector<bigInt> & correctedPhiInt)
{
  bigInt oneOverSix = encode(1/6., deltaPhi, bitsPhi);
  for (size_t i=0; i<phiInt.size(); ++i) {




//    // Only for this term, since DeltaR is much smaller than R, we can preserve 5 bits without risking an overflow.
//    bigInt term1 = multiplyUnit(deltaRadiusInt[i], parInt, registerBits, bitsPhi, reductionPt,
//                                maxBitsMultiplyUnitX+5, maxBitsMultiplyUnitY);
//    // Second order term
//    bigInt temp = multiplyUnit(parInt, radiusInt[i], registerBits, bitsPhi, reductionPt + reductionPhi,
//                               maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
//    bigInt term2 = multiplyUnit(temp, temp, registerBits, bitsPhi, 0, maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
//    term2 = multiplyUnit(temp, term2, registerBits, bitsPhi, 0, maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
//    term2 = multiplyUnit(term2, oneOverSix, registerBits, bitsPhi, 0, maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
//    // Full correction term
//    bigInt fullTerm = adderUnit(term1, term2, registerBits);
//    // bigInt stripCorrection = stripCorrectionTermsInt[i];
//    // alignBits(stripCorrection, 16);
//    // bigInt correctedPhi = adderUnit(phiInt[i], stripCorrection, registerBits);
//    bigInt cPhiInt = adderUnit(phiInt[i], fullTerm, registerBits);
//    correctedPhiInt.push_back(cPhiInt);


    // First order corrections
    bigInt firstOrderTerm = multiplyUnit(deltaRadiusInt[i], parInt, registerBits, bitsPhi, 0,
                                         maxBitsMultiplyUnitX, maxBitsMultiplyUnitY, false);
    bigInt alignedPhi = phiInt[i];
    // Internally in the multiplier it was truncated to 18 bits without any re-alignment (the last input is false)
    alignBits(alignedPhi, -(18-reductionPt));
    bigInt firstOrderCorrectedPhi = adderUnit(firstOrderTerm, alignedPhi, registerBits);
    alignBits(firstOrderCorrectedPhi, (18-reductionPt));

    // Strip index corrections
    // std::cout << "stripCorrectionTermsInt["<<i<<"] = " << stripCorrectionTermsInt[i] << std::endl;
    alignBits(firstOrderCorrectedPhi, -(18-(-8+reductionZ-reductionPhi)));
    firstOrderCorrectedPhi = adderUnit(firstOrderCorrectedPhi, stripCorrectionTermsInt[i], registerBits);
    truncate(firstOrderCorrectedPhi, (18-(-8+reductionZ-reductionPhi)));

    // Second order corrections
    bigInt RcOverTwoRho = multiplyUnit(parInt, radiusInt[i], registerBits, bitsPhi, reductionPt + reductionPhi,
                                       maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
    bigInt RcOverTwoRhoSquared = multiplyUnit(RcOverTwoRho, RcOverTwoRho, registerBits, bitsPhi, 0, maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
    bigInt oneOverSixRcOverTwoRho = multiplyUnit(RcOverTwoRho, oneOverSix, registerBits, bitsPhi, 0, maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
    bigInt oneOverSixRcOverTwoRhoCube = multiplyUnit(oneOverSixRcOverTwoRho, RcOverTwoRhoSquared, registerBits, bitsPhi, 0, maxBitsMultiplyUnitX, maxBitsMultiplyUnitY, false);
    alignBits(firstOrderCorrectedPhi, -18);
    bigInt newCPhiInt = adderUnit(firstOrderCorrectedPhi, oneOverSixRcOverTwoRhoCube, registerBits);
    alignBits(newCPhiInt, 18);
    // bigInt newCPhiInt = adderUnit(phiInt[i], newFullTerm, registerBits);

    // Corrected phi
//    bigInt cPhiInt = adderUnit(phiInt[i], fullTerm, registerBits);
//    correctedPhiInt.push_back(cPhiInt);
    correctedPhiInt.push_back(newCPhiInt);

  }
}


//void zCorrectionUnit(const std::vector<bigInt> & zInt, const double & deltaZ, const int bitsZ,
//                     const std::vector<bigInt> & radiusInt, std::vector<bigInt> & deltaRadiusInt,
//                     const bigInt & parPtInt, const bigInt & parCotThetaInt, const int registerBits,
//                     const int reductionPt, const int reductionCotTheta, const int reductionPhi,
//                     const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
//                     std::vector<bigInt> & correctedZInt)
//{
//  bigInt oneOverSix = encode(1 / 6., 1., bitsZ);
//  for (size_t i = 0; i < zInt.size(); ++i) {
//    // Set to the default bits of the multiplier lowest number + 5 since the deltaR is much smaller than R.
//    bigInt term1 = multiplyUnit(deltaRadiusInt[i], parCotThetaInt, registerBits, bitsZ, reductionCotTheta,
//                                maxBitsMultiplyUnitX+5, maxBitsMultiplyUnitY);
//    // Second order term
//    // The (R*c/(2rho))^2 part is the same as for the phi corrections.
//    bigInt temp = multiplyUnit(parPtInt, radiusInt[i], registerBits, bitsZ, reductionPt + reductionPhi,
//                               maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
//    bigInt term2 = multiplyUnit(temp, temp, registerBits, bitsZ, 0,
//                                maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
//    bigInt RCotTheta = multiplyUnit(parCotThetaInt, radiusInt[i], registerBits, bitsZ, reductionCotTheta,
//                                    maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
//    term2 = multiplyUnit(RCotTheta, term2, registerBits, bitsZ, 0, maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
//    term2 = multiplyUnit(term2, oneOverSix, registerBits, bitsZ, 0, maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
//    // Full correction term
//    // The terms have the negative sign from the pre-estimate being -cot(theta)
//    bigInt fullTerm = adderUnit(term1, term2, registerBits);
//    // Corrected z
//    bigInt cZInt = adderUnit(zInt[i], fullTerm, registerBits);
//    correctedZInt.push_back(cZInt);
//  }
//}


//void zCorrectionUnit(const std::vector<bigInt> & zInt, const double & deltaZ, const int bitsZ,
//                     const std::vector<bigInt> & radiusInt, std::vector<bigInt> & deltaRadiusInt,
//                     const bigInt & parPtInt, const bigInt & parCotThetaInt, const int registerBits,
//                     const int reductionPt, const int reductionCotTheta, const int reductionPhi,
//                     const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
//                     std::vector<bigInt> & correctedZInt)
//{
//  for (size_t i = 0; i < zInt.size(); ++i) {
//    // Set to the default bits of the multiplier lowest number + 5 since the deltaR is much smaller than R.
//    bigInt term1 = multiplyUnit(deltaRadiusInt[i], parCotThetaInt, registerBits, bitsZ, reductionCotTheta,
//                                maxBitsMultiplyUnitX+5, maxBitsMultiplyUnitY);
//    // Corrected z
//    bigInt cZInt = adderUnit(zInt[i], term1, registerBits);
//    correctedZInt.push_back(cZInt);
//  }
//}


void zCorrectionUnit(const std::vector<bigInt> & zInt, const double & deltaZ, const int bitsZ,
                     const std::vector<bigInt> & radiusInt, std::vector<bigInt> & deltaRadiusInt,
                     const bigInt & parPtInt, const bigInt & parCotThetaInt, const int registerBits,
                     const int reductionPt, const int reductionCotTheta, const int reductionPhi,
                     const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                     std::vector<bigInt> & correctedZInt)
{
  bigInt oneOverSix = encode(1 / 6., 1., bitsZ);
  for (size_t i = 0; i < zInt.size(); ++i) {
//    // Set to the default bits of the multiplier lowest number + 5 since the deltaR is much smaller than R.
//    bigInt term1 = multiplyUnit(deltaRadiusInt[i], parCotThetaInt, registerBits, bitsZ, reductionCotTheta,
//                                maxBitsMultiplyUnitX+5, maxBitsMultiplyUnitY);

    // New parallel method
    bigInt firstOrderTerm = multiplyUnit(deltaRadiusInt[i], parCotThetaInt, registerBits, bitsZ, 0,
                                         maxBitsMultiplyUnitX, maxBitsMultiplyUnitY, false);
    bigInt alignedZ = zInt[i];
    // Internally in the multiplier it was truncated to 18 bits without any re-alignment (the last input is false)
    alignBits(alignedZ, -(18-reductionCotTheta));
    bigInt firstOrderCorrectedZ = adderUnit(firstOrderTerm, alignedZ, registerBits);
    alignBits(firstOrderCorrectedZ, (18-reductionCotTheta));

    bigInt RCotTheta = multiplyUnit(parCotThetaInt, radiusInt[i], registerBits, bitsZ, reductionCotTheta,
                                    maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
    bigInt RcOverTwoRho = multiplyUnit(parPtInt, radiusInt[i], registerBits, bitsZ, reductionPt + reductionPhi,
                                       maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
    bigInt RcOverTwoRhoSquared = multiplyUnit(RcOverTwoRho, RcOverTwoRho, registerBits, bitsZ, 0, maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
    bigInt oneOverSixRCotTheta = multiplyUnit(RCotTheta, oneOverSix, registerBits, bitsZ, 0, maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);

    bigInt secondOrderTerm = multiplyUnit(oneOverSixRCotTheta, RcOverTwoRhoSquared, registerBits, bitsZ, 0, maxBitsMultiplyUnitX, maxBitsMultiplyUnitY, false);
    alignBits(firstOrderCorrectedZ, -18);
    bigInt cZInt = adderUnit(firstOrderCorrectedZ, secondOrderTerm, registerBits);
    alignBits(cZInt, 18);

//    // Second order term
//    // The (R*c/(2rho))^2 part is the same as for the phi corrections.
//    bigInt temp = multiplyUnit(parPtInt, radiusInt[i], registerBits, bitsZ, reductionPt + reductionPhi,
//                               maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
//    bigInt term2 = multiplyUnit(temp, temp, registerBits, bitsZ, 0,
//                                maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
//    bigInt RCotTheta = multiplyUnit(parCotThetaInt, radiusInt[i], registerBits, bitsZ, reductionCotTheta,
//                                    maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
//    term2 = multiplyUnit(RCotTheta, term2, registerBits, bitsZ, 0, maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
//    term2 = multiplyUnit(term2, oneOverSix, registerBits, bitsZ, 0, maxBitsMultiplyUnitX, maxBitsMultiplyUnitY);
//    // Full correction term
//    // The terms have the negative sign from the pre-estimate being -cot(theta)
//    bigInt fullTerm = adderUnit(term1, term2, registerBits);
//    // Corrected z
//    bigInt cZInt = adderUnit(zInt[i], fullTerm, registerBits);
//    // Corrected z
//    bigInt cZInt = adderUnit(zInt[i], term1, registerBits);
    correctedZInt.push_back(cZInt);
  }
}


bigInt extrapolateRFirstOrderEmulator(const bigInt & R, const bigInt & z, const int layer, const bigInt & tgTheta,
                                      const bigInt & chargeOverTwoRho, const std::vector<int> & uniqueLayers,
                                      const std::vector<bigInt> & originalRInt,
                                      const std::vector<bigInt> & originalZInt,
                                      const bigInt & Rcut,
                                      const int registerBits, const int bitsVar, const double & scaleZ,
                                      const int reductionPhi, const int reductionZ,
                                      const int reductionPt, const int reductionTgTheta,
                                      const int maxBitsX, const int maxBitsY)
{
  if (layer > 10 && R > Rcut) {
    // Extrapolate R from the outermost PS module in this stubs combination
    for (int i = int(uniqueLayers.size()) - 1; i >= 0; --i) {
      if (uniqueLayers[i] < 8 || (uniqueLayers[i] > 10 && originalRInt[i] < Rcut)) {
        bigInt deltaZ = adderUnit(z, -originalZInt[i], registerBits);
        // bigInt deltaZTgTheta = multiplyUnit(tgTheta, deltaZ, registerBits, bitsVar, reductionTgTheta, maxBitsX, maxBitsY, true);
        // bigInt deltaZTgTheta_noShift = multiplyUnit(tgTheta, deltaZ, registerBits, bitsVar, reductionTgTheta, maxBitsX, maxBitsY, false);
        bigInt deltaZTgTheta_noShift_reversed = multiplyUnit(deltaZ, tgTheta, registerBits, bitsVar, reductionTgTheta, maxBitsX, maxBitsY, false);
        // bigInt extrapolatedRInt = adderUnit(originalRInt[i], deltaZTgTheta, registerBits);
        bigInt alignedOriginalR = originalRInt[i];
        alignBits(alignedOriginalR, -11);
        bigInt extrapolatedRInt_noShift_reversed = adderUnit(alignedOriginalR, deltaZTgTheta_noShift_reversed, registerBits);
        alignBits(extrapolatedRInt_noShift_reversed, 11);
        return extrapolatedRInt_noShift_reversed;
        // return extrapolatedRInt;
      }
    }
  }
  return R;
}


bigInt extrapolateRSecondOrderEmulator(const bigInt & R, const bigInt & z, const int layer, const bigInt & tgTheta,
                                       const bigInt & chargeOverTwoRho, const std::vector<int> & uniqueLayers,
                                       const std::vector<bigInt> & originalRInt,
                                       const std::vector<bigInt> & originalZInt,
                                       const bigInt & Rcut,
                                       const int registerBits, const int bitsVar, const double & scaleZ,
                                       const int reductionPhi, const int reductionZ,
                                       const int reductionPt, const int reductionTgTheta,
                                       const int maxBitsX, const int maxBitsY)
{
  if (layer > 10 && R > Rcut) {
    // Extrapolate R from the outermost PS module in this stubs combination
    for (int i = int(uniqueLayers.size()) - 1; i >= 0; --i) {
      if (uniqueLayers[i] < 8 || (uniqueLayers[i] > 10 && originalRInt[i] < Rcut)) {
        bigInt deltaZ = adderUnit(z, -originalZInt[i], registerBits);
        bigInt deltaZTgTheta = multiplyUnit(tgTheta, deltaZ, registerBits, bitsVar, reductionTgTheta, maxBitsX, maxBitsY, true);
        bigInt deltaZTgThetaChargeOverTwoRho = multiplyUnit(chargeOverTwoRho, deltaZTgTheta, registerBits, bitsVar,
                                                            reductionPt + reductionPhi, maxBitsX, maxBitsY, true);
        bigInt originalRChargeOverTwoRho = multiplyUnit(chargeOverTwoRho, originalRInt[i], registerBits, bitsVar,
                                                        reductionPt + reductionPhi, maxBitsX, maxBitsY, true);
        bigInt partialTerm1 = multiplyUnit(deltaZTgThetaChargeOverTwoRho, deltaZTgThetaChargeOverTwoRho, registerBits, bitsVar, 0, maxBitsX, maxBitsY, true);
        bigInt partialTerm2 = multiplyUnit(originalRChargeOverTwoRho, originalRChargeOverTwoRho, registerBits, bitsVar, 0, maxBitsX, maxBitsY, true);
        // Note that we use bitsVar+1 instead of bitsVar to also do the division by 2
        bigInt term1 = - multiplyUnit(partialTerm2, deltaZTgTheta, registerBits, bitsVar+1, 0, maxBitsX, maxBitsY, true);
        bigInt term2 = - multiplyUnit(partialTerm1, originalRInt[i], registerBits, bitsVar+1, 0, maxBitsX, maxBitsY, true);
        bigInt oneOverThree = encode(1/3., scaleZ, bitsVar);
        bigInt deltaZTgThetaOverThree = multiplyUnit(deltaZTgTheta, oneOverThree, registerBits, bitsVar, reductionZ, maxBitsX, maxBitsY, true);
        bigInt term3 = - multiplyUnit(deltaZTgThetaOverThree, partialTerm1, registerBits, bitsVar+1, 0, maxBitsX, maxBitsY, true);
        bigInt secondOrderTerm = adderUnit(adderUnit(term1, term2, registerBits), term3, registerBits);
        bigInt extrapolatedRInt = adderUnit(adderUnit(originalRInt[i], deltaZTgTheta, registerBits), secondOrderTerm, registerBits);
        return extrapolatedRInt;
      }
    }
  }
  return R;
}


bigInt extrapolateRSecondOrderFirstTwoTermsOnlyEmulator(const bigInt & R, const bigInt & z, const int layer, const bigInt & tgTheta,
                                                        const bigInt & chargeOverTwoRho, const std::vector<int> & uniqueLayers,
                                                        const std::vector<bigInt> & originalRInt,
                                                        const std::vector<bigInt> & originalZInt,
                                                        const bigInt & Rcut,
                                                        const int registerBits, const int bitsVar, const double & scaleZ,
                                                        const int reductionPhi, const int reductionZ,
                                                        const int reductionPt, const int reductionTgTheta,
                                                        const int maxBitsX, const int maxBitsY)
{
  if (layer > 10 && R > Rcut) {
    // Extrapolate R from the outermost PS module in this stubs combination
    for (int i = int(uniqueLayers.size()) - 1; i >= 0; --i) {
      if (uniqueLayers[i] < 8 || (uniqueLayers[i] > 10 && originalRInt[i] < Rcut)) {
        bigInt deltaZ = adderUnit(z, -originalZInt[i], registerBits);
        bigInt deltaZTgTheta = multiplyUnit(tgTheta, deltaZ, registerBits, bitsVar, reductionTgTheta, maxBitsX, maxBitsY, true);
        bigInt deltaZTgThetaChargeOverTwoRho = multiplyUnit(chargeOverTwoRho, deltaZTgTheta, registerBits, bitsVar,
                                                            reductionPt + reductionPhi, maxBitsX, maxBitsY, true);
        bigInt originalRChargeOverTwoRho = multiplyUnit(chargeOverTwoRho, originalRInt[i], registerBits, bitsVar,
                                                        reductionPt + reductionPhi, maxBitsX, maxBitsY, true);
        bigInt partialTerm1 = multiplyUnit(deltaZTgThetaChargeOverTwoRho, deltaZTgThetaChargeOverTwoRho, registerBits, bitsVar, 0, maxBitsX, maxBitsY, true);
        bigInt partialTerm2 = multiplyUnit(originalRChargeOverTwoRho, originalRChargeOverTwoRho, registerBits, bitsVar, 0, maxBitsX, maxBitsY, true);
        // Note that we use bitsVar+1 instead of bitsVar to also do the division by 2
        bigInt term1 = - multiplyUnit(partialTerm2, deltaZTgTheta, registerBits, bitsVar+1, 0, maxBitsX, maxBitsY, true);
        bigInt term2 = - multiplyUnit(partialTerm1, originalRInt[i], registerBits, bitsVar+1, 0, maxBitsX, maxBitsY, true);
//        bigInt oneOverThree = encode(1/3., scaleZ, bitsVar);
        bigInt secondOrderTerm = adderUnit(term1, term2, registerBits);
        bigInt extrapolatedRInt = adderUnit(adderUnit(originalRInt[i], deltaZTgTheta, registerBits), secondOrderTerm, registerBits);
        return extrapolatedRInt;
      }
    }
  }
  return R;
}


bigInt extrapolateRSecondOrderFirstTermOnlyEmulator(const bigInt & R, const bigInt & z, const int layer, const bigInt & tgTheta,
                                                    const bigInt & chargeOverTwoRho, const std::vector<int> & uniqueLayers,
                                                    const std::vector<bigInt> & originalRInt,
                                                    const std::vector<bigInt> & originalZInt,
                                                    const bigInt & Rcut,
                                                    const int registerBits, const int bitsVar, const double & scaleZ,
                                                    const int reductionPhi, const int reductionZ,
                                                    const int reductionPt, const int reductionTgTheta,
                                                    const int maxBitsX, const int maxBitsY)
{
  if (layer > 10 && R > Rcut) {
    // Extrapolate R from the outermost PS module in this stubs combination
    for (int i = int(uniqueLayers.size()) - 1; i >= 0; --i) {
      if (uniqueLayers[i] < 8 || (uniqueLayers[i] > 10 && originalRInt[i] < Rcut)) {
        bigInt deltaZ = adderUnit(z, -originalZInt[i], registerBits);
        bigInt deltaZTgTheta = multiplyUnit(tgTheta, deltaZ, registerBits, bitsVar, reductionTgTheta, maxBitsX, maxBitsY, true);
        bigInt originalRChargeOverTwoRho = multiplyUnit(chargeOverTwoRho, originalRInt[i], registerBits, bitsVar,
                                                        reductionPt + reductionPhi, maxBitsX, maxBitsY, true);
        bigInt partialTerm2 = multiplyUnit(originalRChargeOverTwoRho, originalRChargeOverTwoRho, registerBits, bitsVar, 0, maxBitsX, maxBitsY, true);
        // Note that we use bitsVar+1 instead of bitsVar to also do the division by 2
        bigInt term1 = - multiplyUnit(partialTerm2, deltaZTgTheta, registerBits, bitsVar+1, 0, maxBitsX, maxBitsY, true);
        bigInt extrapolatedRInt = adderUnit(adderUnit(originalRInt[i], deltaZTgTheta, registerBits), term1, registerBits);
        return extrapolatedRInt;
      }
    }
  }
  return R;
}


CorrectPhiForNonRadialStripsEmulator::CorrectPhiForNonRadialStripsEmulator() : stripMean_(508)
{
  // Encode the pitch as a 18 bits number. The range used is 2^(-6) = 0.015625.
  stripPitch_ = encode(0.009, std::pow(2, -5), 18);
  // Lookup table for 0.009/R^2
  lookupOneOverRSquared_.insert(std::make_pair(31, 154));
  lookupOneOverRSquared_.insert(std::make_pair(51, 57));
  lookupOneOverRSquared_.insert(std::make_pair(53, 52));
  lookupOneOverRSquared_.insert(std::make_pair(33, 132));
  lookupOneOverRSquared_.insert(std::make_pair(37, 109));
  lookupOneOverRSquared_.insert(std::make_pair(48, 64));
  lookupOneOverRSquared_.insert(std::make_pair(34, 125));
  lookupOneOverRSquared_.insert(std::make_pair(42, 83));
  lookupOneOverRSquared_.insert(std::make_pair(45, 72));
  lookupOneOverRSquared_.insert(std::make_pair(43, 81));
  lookupOneOverRSquared_.insert(std::make_pair(40, 94));
  lookupOneOverRSquared_.insert(std::make_pair(50, 58));
  // lookupOneOverRSquared_.insert(std::make_pair(51, 58));
}


//bigInt CorrectPhiForNonRadialStripsEmulator::correctPhiForNonRadialStrips(const bigInt & phi, const bigInt & stripIndex,
//                                                                          const bigInt & extrapolatedR, const bigInt & R,
//                                                                          const bigInt & RCut, const bigInt & oneOverR,
//                                                                          const int layer, const int registerBits, const int bitsVar,
//                                                                          const int maxBitsX, const int maxBitsY,
//                                                                          const int reductionZ, const int reductionPhi)
//{
//  if (layer > 10 && R > RCut) {
//    bigInt deltaStripIndex = adderUnit(stripIndex, -stripMean_, registerBits);
//    bigInt deltaR = adderUnit(extrapolatedR, -R, registerBits);
//    // TODO: The oneOverR should be taken from a lookup table
//    // Scale it by reductionZ ( = reductionR) to bring it to scaleR (it is scaleR^2).
//    bigInt deltaROverR = multiplyUnit(oneOverR, deltaR, registerBits, bitsVar, reductionZ, maxBitsX, maxBitsY, true);
//    // The strip index is a number from 0 to 1015. No encoding and no shift required
//    // (equivalent to a scale factor equal to the bitshift).
//    bigInt deltaStripIndexOverR = multiplyUnit(oneOverR, deltaStripIndex, registerBits, 0, 0, maxBitsX, maxBitsY, true);
//    // The order of this operation is critical for the performance of the 27x18 multiplication.
//    bigInt correctionTerm = multiplyUnit(deltaROverR, deltaStripIndexOverR, registerBits, bitsVar, reductionZ, maxBitsX, maxBitsY, true);
//    correctionTerm = multiplyUnit(correctionTerm, stripPitch_, registerBits, 18, -5+reductionZ-reductionPhi, maxBitsX, maxBitsY, true);
//    bigInt correctedPhi = adderUnit(phi, -correctionTerm, registerBits);
//    return correctedPhi;
//  }
//  return phi;
//}


//bigInt CorrectPhiForNonRadialStripsEmulator::correctPhiForNonRadialStrips(const bigInt & phi, const bigInt & stripIndex,
//                                                                          const bigInt & extrapolatedR, const bigInt & R,
//                                                                          const bigInt & RCut, const bigInt & oneOverR,
//                                                                          const bigInt & oneOverRSquared,
//                                                                          const int layer, const int registerBits, const int bitsVar,
//                                                                          const int maxBitsX, const int maxBitsY,
//                                                                          const int reductionZ, const int reductionPhi)
//{
//  if (layer > 10 && R > RCut) {
//    if (stripIndex > 1015) {
//      std::cout << "Error: Strip index = " << stripIndex << std::endl;
//      std::cout << "Maximum strip index should be 1015 for 2S modules in the disks." << std::endl;
//      throw;
//    }
//    bigInt deltaStripIndex = adderUnit(stripIndex, -stripMean_, registerBits);
//    bigInt term1 = multiplyUnit(oneOverRSquared, deltaStripIndex, registerBits, 0, 0, maxBitsX, maxBitsY, true);
//    bigInt deltaR = adderUnit(extrapolatedR, -R, registerBits);
//    bigInt correctionTerm = multiplyUnit(deltaR, term1, registerBits, 18, -8+reductionZ-reductionPhi, maxBitsX, maxBitsY, true);
//    double correctionTermFloat = decode(correctionTerm, 2, 27);
//    // TODO: The oneOverR should be taken from a lookup table
//    // Scale it by reductionZ ( = reductionR) to bring it to scaleR (it is scaleR^2).
//    bigInt deltaROverR = multiplyUnit(oneOverR, deltaR, registerBits, bitsVar, reductionZ, maxBitsX, maxBitsY, true);
//    // The strip index is a number from 0 to 1015. No encoding and no shift required
//    // (equivalent to a scale factor equal to the bitshift).
//    bigInt deltaStripIndexOverR = multiplyUnit(oneOverR, deltaStripIndex, registerBits, 0, 0, maxBitsX, maxBitsY, true);
//    // The order of this operation is critical for the performance of the 27x18 multiplication.
//    bigInt correctionTermOld = multiplyUnit(deltaROverR, deltaStripIndexOverR, registerBits, bitsVar, reductionZ, maxBitsX, maxBitsY, true);
//    correctionTermOld = multiplyUnit(correctionTermOld, stripPitch_, registerBits, 18, -5+reductionZ-reductionPhi, maxBitsX, maxBitsY, true);
//    double correctionTermOldFloat = decode(correctionTermOld, 2, 27);
//    bigInt correctedPhi = adderUnit(phi, -correctionTerm, registerBits);
//    return correctedPhi;
//  }
//  return phi;
//}


bigInt CorrectPhiForNonRadialStripsEmulator::correctPhiForNonRadialStrips(const bigInt & phi, const bigInt & stripIndex,
                                                                          const bigInt & extrapolatedR, const bigInt & R,
                                                                          const bigInt & RCut,
                                                                          const int layer, const int registerBits, const int bitsVar,
                                                                          const int maxBitsX, const int maxBitsY,
                                                                          const int reductionZ, const int reductionPhi)
{
  if (layer > 10 && R > RCut) {
    if (stripIndex > 1015) {
      std::cout << "Error: Strip index = " << stripIndex << std::endl;
      std::cout << "Maximum strip index should be 1015 for 2S modules in the disks." << std::endl;
      throw;
    }
    bigInt deltaStripIndex = adderUnit(stripIndex, -stripMean_, registerBits);
    bigInt lookupR = R;
    truncate(lookupR, 18);
    auto it = lookupOneOverRSquared_.find(lookupR);
    if (it == lookupOneOverRSquared_.end()) {
      std::cout << "Error: lookupR = " << lookupR << " not found for R = " << R << std::endl;
      throw;
    }
//    bigInt term1 = multiplyUnit(oneOverRSquared, deltaStripIndex, registerBits, 0, 0, maxBitsX, maxBitsY, true);
    bigInt term1 = multiplyUnit(it->second, deltaStripIndex, registerBits, 0, 0, maxBitsX, maxBitsY, true);
    bigInt negativeDeltaR = adderUnit(R, -extrapolatedR, registerBits);
    // 0.009/R^2 is encoded with 18 bits and a range of 2^-8.
//    bigInt correctionTerm = multiplyUnit(deltaR, term1, registerBits, 18, -8+reductionZ-reductionPhi, maxBitsX, maxBitsY, true);
    // bigInt correctionTermAligned = multiplyUnit(deltaR, term1, registerBits, 18, -8+reductionZ-reductionPhi, maxBitsX, maxBitsY, true);
//    bigInt correctionTerm = multiplyUnit(negativeDeltaR, term1, registerBits, 18, 0, maxBitsX, maxBitsY, false);

    return multiplyUnit(negativeDeltaR, term1, registerBits, 18, 0, maxBitsX, maxBitsY, false);

    // double decodedCorrectionTerm = decode(correctionTermAligned, 1., 18);
    // double decodedCorrectionTermAlignedNew = decode(correctionTerm, 4., 27+18);

//    bigInt correctedPhi = phi;
//    alignBits(correctedPhi, -16);
//    correctedPhi = adderUnit(correctedPhi, correctionTerm, registerBits);
//    truncate(correctedPhi, 16);

//    bigInt correctedPhi = adderUnit(phi, -correctionTerm, registerBits);
    // double decodedCorrectedPhi = decode(correctedPhi, 2., 27);
//    return correctedPhi;
//    return correctionTerm;
  }
//  return phi;
  return 0;
}


int computeShiftBits(const int bitsX, const int reducedAccuracyBitsX)
{
  int shiftBitsX = 0;
  if (reducedAccuracyBitsX != 0) {
    shiftBitsX = bitsX - reducedAccuracyBitsX;
    if (shiftBitsX < 0) {
      std::cout <<
      "Error: reducedAccuracyBitsX is bigger than bitsX. The allowed range is 0 <= reducedAccuracyBitsX < bitsX" <<
      std::endl;
      std::cout << "bitsX = " << bitsX << std::endl;
      std::cout << "reducedAccuracyBitsX = " << reducedAccuracyBitsX << std::endl;
      throw;
    }
  }
  return shiftBitsX;
}


std::string convertToBitString(const bigInt & value, const int bits)
{
  std::string bitString;
  int sign = value > 0 ? 1 : -1;
  boost::dynamic_bitset<> bitsValue(bits, abs(int(value)));
  // Compute the negative under two's complement notation if the original sign is negative
  // We do it explicitly to make sure it works with input types different from signed long
  if (sign == -1) bitsValue = boost::dynamic_bitset<>(bits, bitsValue.flip().to_ulong() + 1);
  boost::to_string(bitsValue, bitString);
  return bitString;
}
