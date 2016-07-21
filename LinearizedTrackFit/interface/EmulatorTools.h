//
// Created by Marco De Mattia on 9/7/15.
//

#ifndef REMOTEPROJECTS_EMULATORTOOLS_H_H
#define REMOTEPROJECTS_EMULATORTOOLS_H_H

// #include <boost/multiprecision/cpp_int.hpp>
#include <boost/algorithm/string/replace.hpp>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GetVariables.h"
#include <boost/dynamic_bitset.hpp>

// typedef boost::multiprecision::int128_t bigInt;
// typedef long long int bigInt;
typedef int64_t bigInt;


struct MaxDeltaAndFactors
{
  MaxDeltaAndFactors() : fixValues(false)
  {}

  MaxDeltaAndFactors(std::vector<double> & inputMaxDeltaD, std::vector<int> & inputMaxFactorD,
                     std::vector<double> & inputMaxDeltaV, std::vector<int> & inputMaxFactorV) :
      maxDeltaD(std::move(inputMaxDeltaD)), maxFactorD(std::move(inputMaxFactorD)),
      maxDeltaV(std::move(inputMaxDeltaV)), maxFactorV(std::move(inputMaxFactorV)),
      fixValues(true)
  {}

  MaxDeltaAndFactors(const double & inputMaxDeltaD, const int inputMaxFactorD) :
      fixValues(true)
  {
    maxDeltaD.push_back(inputMaxDeltaD);
    maxFactorD.push_back(inputMaxFactorD);
  }

  std::vector<double> maxDeltaD;
  std::vector<int> maxFactorD;
  std::vector<double> maxDeltaV;
  std::vector<int> maxFactorV;
  bool fixValues;
};


bigInt encode(const double & x, const double & deltaX, const int bits);


template <class T>
void encodeVector(const T & vec, std::vector<bigInt> & vecInt, const double & deltaX, const int bits)
{
  for (int row=0; row<vec.rows(); ++row) {
    vecInt.push_back(encode(vec(row, 0), deltaX, bits));
  }
}


/// Specialized version for std::vector<double>
template <>
void encodeVector(const std::vector<double> & vec, std::vector<bigInt> & vecInt, const double & deltaX, const int bits);


template <class T>
void encodeVector(const T & vec, std::vector<bigInt> & vecInt, const double & deltaX,
                  const std::vector<double> & deltaA, const int bits)
{
  if (vec.rows() != int(deltaA.size())) {
    std::cout << "encodeVector error: vec and deltaA have different size:" << std::endl;
    std::cout << "vec size = " << vec.size() << ", deltaA size = " << deltaA.size() << std::endl;
    throw;
  }
  for (int i=0; i<vec.rows(); ++i) {
//    std::cout << "encoding coefficient " << vec.at(i) << " to int = " << encode(vec.at(i), deltaX, bits) << std::endl;
    vecInt.push_back(encode(vec(i, 0), deltaX*deltaA.at(i), bits));
  }
}


void truncate(bigInt & value, const int truncateBits);


void alignBits(bigInt & x, const int bitShift);


std::pair<int, double> computeFactor(const std::vector<double> & coefficients, const double & deltaX,
                                     const double & scaleFactor=1., const bool powerTwoRange=true);


template <class T>
void encodeMatrix(const T & matrix, std::vector<std::vector<bigInt> > & matrixInt, const double & deltaX,
                  std::vector<double> & deltaA, std::vector<int> & factors, const int bits,
                  const std::vector<double> & maxDeltas, const std::vector<int> & maxFactors, const bool Vcoeff,
                  const std::vector<bool> & powerTwoRanges)
{
  // If we fix the ranges used to encode the coefficients to be always the same for same rows
  // across different combinations we shift the chi2 terms by one. This is because the smallest
  // chi2 term is the one with the least amount of information and considering this as the extra
  // term in the 6/6 case minimizes the variation in the coefficient ranges for the chi2 terms.
  // We therefore align the coefficient terms as:
  // 6/6 [0], [1], [2], [3]
  // 5/6      [9], [1], [3]
  bool maxRanges = false;
  if (maxDeltas.size() != 0) {
    maxRanges = true;
    int shift = 0;
    if (Vcoeff && matrix.rows() == 5) shift = 1;
    std::copy(maxDeltas.begin()+shift, maxDeltas.end(), std::back_inserter(deltaA));
    std::copy(maxFactors.begin()+shift, maxFactors.end(), std::back_inserter(factors));
  }

  // Align the chi2 coefficients for the 5/6 cases to the 6/6
  for (int row=0; row<matrix.rows(); ++row) {
    std::vector<double> coefficients;
    for (int col=0; col<matrix.cols(); ++col) {
      // Find the value of deltaA and the corresponding factor for this column.
      coefficients.push_back(matrix(row, col));
    }
    if (!maxRanges) {
      auto result(computeFactor(coefficients, deltaX, 1., powerTwoRanges.at(row)));
      factors.push_back(result.first);
      deltaA.push_back(result.second);
    }

    std::vector<bigInt> coefficientsInt;
    encodeVector(coefficients, coefficientsInt, deltaA.at(row), bits);
    matrixInt.push_back(coefficientsInt);
  }
}


double decode(const bigInt & i, const double & deltaX, const int bits);

void checkSizeConsistency(const std::vector<bigInt> & a, const std::vector<bigInt> & b);

void checkBufferOverflow(const bigInt & value, const int bits, const std::string & errorMessage="Buffer overflow");

void preDiffer(const std::vector<bigInt> & a, const std::vector<bigInt> & b, std::vector<bigInt> & c, const int bits=0);


/**
 * Multiplies two numbers and checks for buffer overflow.
 * The maximum number of bits of each number are constrained to be no more
 * than the limits. The default values are 17 bits for the first number and
 * 27 bits for the second.
 */
bigInt multiplyUnit(const bigInt & x, const bigInt & y, const int registerBits, const int bitsX,
                    const int reduction, const int maxBitsX, const int maxBitsY, const bool doBitShift=true);

bigInt macUnit(const std::vector<bigInt> & coordinates, const std::vector<bigInt> & coefficients,
               const bigInt & parMean, const int registerBits, const int bitsCoeff, const int bitsMultiplyUnitX,
               const int bitsMultiplyUnitY);

bigInt adderUnit(const bigInt & x, const bigInt & y, const int bits);


//void phiCorrectionUnit(const std::vector<bigInt> & phiInt, const double & deltaPhi, const int bitsPhi,
//                       const std::vector<bigInt> & radiusInt, std::vector<bigInt> & deltaRadiusInt,
//                       const bigInt & parInt, const int registerBits, const int reduction, const int reductionPhi,
//                       const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
//                       std::vector<bigInt> & correctedPhiInt);
void phiCorrectionUnit(const std::vector<bigInt> & phiInt, const double & deltaPhi, const int bitsPhi,
                       const std::vector<bigInt> & radiusInt, std::vector<bigInt> & deltaRadiusInt,
                       const bigInt & parInt, const int registerBits, const int reductionPt, const int reductionPhi,
                       const int reductionZ, const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                       const std::vector<bigInt> & stripCorrectionTermsInt,
                       std::vector<bigInt> & correctedPhiInt);


void zCorrectionUnit(const std::vector<bigInt> & zInt, const double & deltaZ, const int bitsZ,
                     const std::vector<bigInt> & radiusInt, std::vector<bigInt> & deltaRadiusInt,
                     const bigInt & parPtInt, const bigInt & parCotThetaInt, const int registerBits,
                     const int reductionPt, const int reductionCotTheta, const int reductionPhi,
                     const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                     std::vector<bigInt> & correctedZInt);


/// Emulator version of the EstimatorSimple utilized to compute the pre-estimates
class EstimatorSimpleEmulator : public EstimatorSimple
{
 public:
  EstimatorSimpleEmulator(const TString & inputFileName, const double & deltaX, const double & deltaR,
                          const int registerBits, const int bitsX, const int bitsA,
                          const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                          const std::vector<double> & deltaA, const std::vector<int> & factorA,
                          const double scaleFactor = 1., const double & ptSplitValue = 10.) :
      EstimatorSimple(inputFileName, scaleFactor), // deltaX_(deltaX),
      registerBits_(registerBits), bitsX_(bitsX), bitsA_(bitsA),
      maxBitsMultiplyUnitX_(maxBitsMultiplyUnitX), maxBitsMultiplyUnitY_(maxBitsMultiplyUnitY),
      inputFileName_(inputFileName)
  {
    // Compute the deltaA and the reduction "factor" given the range of A and R such that
    // deltaA*deltaR is a power of 2^factor.
    if (deltaA.size() > 0) {
      deltaA_ = deltaA.at(0);
      reduction_ = factorA.at(0);
    }
    else {
      auto result(computeFactor(coeff_, deltaR));
      reduction_ = result.first;
      deltaA_ = result.second;
    }

    if (ptSplitValue == 0.) {
      std::cout << "EstimatorSimpleEmulator error: ptSplitValue = " << ptSplitValue << std::endl;
      throw;
    }

    chargeOverTwoRhoSplitValueInt_ = encode(scaleFactor/ptSplitValue, deltaX*deltaA_, bitsX_);

    parameterMeanInt_ = encode(parameterMean_, deltaX*deltaA_, bitsX_+bitsA_);
    if (means_.size() < 10) {
      encodeVector(means_, meansInt_, deltaX, bitsX_);
      encodeVector(coeff_, coeffInt_, deltaA_, bitsA_);
    }
    else {
      std::vector<double> means1;
      std::vector<double> means2;
      std::vector<double> coeff1;
      std::vector<double> coeff2;
      for (size_t i=0; i<means_.size()/2; ++i) {
        means1.push_back(means_[i*2]);
        means2.push_back(means_[i*2+1]);
        coeff1.push_back(coeff_[i*2]);
        coeff2.push_back(coeff_[i*2+1]);
      }
      encodeVector(means1, meansInt_, deltaX, bitsX_);
      encodeVector(means2, meansInt2_, deltaX, bitsX_);
      encodeVector(coeff1, coeffInt_, deltaA_, bitsA_);
      encodeVector(coeff2, coeffInt_, deltaA_, bitsA_);
      encodeVector(coeff1, coeffInt1_, deltaA_, bitsA_);
      encodeVector(coeff2, coeffInt2_, deltaA_, bitsA_);
    }
    bool writeCoefficients = false;
    if (writeCoefficients) write();
  }

  int reduction() const { return reduction_; }
  bigInt chargeOverTwoRhoSplitValue() const { return chargeOverTwoRhoSplitValueInt_; }
  double deltaA() const { return deltaA_; }

  template <class T>
  bigInt estimate(const T & var)
  {
    std::vector<bigInt> varDeltas;
    preDiffer(var, meansInt_, varDeltas, registerBits_);
    return macUnit(varDeltas, coeffInt_, parameterMeanInt_, registerBits_, bitsA_,
                   maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_);
  }

  template <class T, class U>
  bigInt estimate(const T & var1, const U & var2)
  {
    std::vector<bigInt> varDeltas;
    preDiffer(var1, meansInt_, varDeltas, registerBits_);
    preDiffer(var2, meansInt2_, varDeltas, registerBits_);
    std::vector<bigInt> varDeltas1;
    preDiffer(var1, meansInt_, varDeltas1, registerBits_);
    bigInt first = macUnit(varDeltas1, coeffInt1_, parameterMeanInt_, registerBits_, 0,
                           maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_);
    std::vector<bigInt> varDeltas2;
    preDiffer(var2, meansInt2_, varDeltas2, registerBits_);
    bigInt second = macUnit(varDeltas2, coeffInt2_, 0, registerBits_, 0,
                            maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_);
    bigInt outputInt = adderUnit(first, second, registerBits_);
    truncate(outputInt, bitsA_);
    return outputInt;
  }

  void write()
  {
    std::string fileName = inputFileName_;
    // Sanitize the file name.
    boost::replace_all(fileName, "//", "/");
    std::size_t foundOne = fileName.rfind("/");
    if (foundOne != std::string::npos) {
      std::size_t found = fileName.rfind("/", foundOne-2);
      if (found != std::string::npos) {
        fileName = fileName.replace(foundOne, 1, "_").substr(found + 1);
      }
      else {
        std::cout << "MatrixReaderEmulator: Error building output file name from input file name = " << inputFileName_ << std::endl;
        throw;
      }
    }

    std::ofstream outfile;
    outfile.open("Emulator_"+fileName);
    if(!outfile) {
      std::cout << "error opening Emulator"+fileName << std::endl;
      return;
    }

    outfile << "Mean values:" << std::endl;
    for (const auto & meanValue : meansInt_) {
      outfile << meanValue << std::endl;
    }
    outfile << "Mean values 2:" << std::endl;
    for (const auto & meanValue : meansInt2_) {
      outfile << meanValue << std::endl;
    }
    outfile << "Coefficients:" << std::endl;
    for (const auto & coeff : coeffInt_) {
      outfile << coeff << std::endl;
    }
    outfile << "Parameter mean:" << std::endl;
    outfile << parameterMeanInt_ << std::endl;
    outfile << "DeltaA:" << std::endl;
    outfile << deltaA_ << std::endl;
    outfile << "Reduction:" << std::endl;
    // Output in a format already usable by the firmware.
    outfile << bitsA_ - reduction_ << std::endl;
    outfile << "chargeOverTwoRhoSplit:" << std::endl;
    outfile << chargeOverTwoRhoSplitValueInt_ << std::endl;
    outfile << std::endl;
    outfile.close();
  }

 private:
  std::vector<bigInt> meansInt_;
  std::vector<bigInt> meansInt2_;
  std::vector<bigInt> coeffInt_;
  std::vector<bigInt> coeffInt1_;
  std::vector<bigInt> coeffInt2_;
  double deltaA_;
  bigInt parameterMeanInt_;
  int registerBits_;
  int bitsX_;
  int bitsA_;
  int reduction_;
  bigInt chargeOverTwoRhoSplitValueInt_;
  const int maxBitsMultiplyUnitX_;
  const int maxBitsMultiplyUnitY_;
  std::string inputFileName_;
};


// Extrapolate R to from the outermost PS module to the given z position using the given tgTheta.
// Use a first order approximation.
bigInt extrapolateRFirstOrderEmulator(const bigInt & R, const bigInt & z, const int layer, const bigInt & tgTheta,
                                      const bigInt & chargeOverTwoRho, const std::vector<int> & uniqueLayers,
                                      const std::vector<bigInt> & originalRInt,
                                      const std::vector<bigInt> & originalZInt,
                                      const bigInt & Rcut,
                                      const int registerBits, const int bitsVar, const double & scaleZ,
                                      const int reductionPhi, const int reductionZ,
                                      const int reductionPt, const int reductionTgTheta,
                                      const int maxBitsX, const int maxBitsY);

// Extrapolate R to from the outermost PS module to the given z position using the given tgTheta.
// Use a second order approximation with only one term.
bigInt extrapolateRSecondOrderFirstTermOnlyEmulator(const bigInt & R, const bigInt & z, const int layer, const bigInt & tgTheta,
                                                    const bigInt & chargeOverTwoRho, const std::vector<int> & uniqueLayers,
                                                    const std::vector<bigInt> & originalRInt,
                                                    const std::vector<bigInt> & originalZInt,
                                                    const bigInt & Rcut,
                                                    const int registerBits, const int bitsVar, const double & scaleZ,
                                                    const int reductionPhi, const int reductionZ,
                                                    const int reductionPt, const int reductionCotTheta,
                                                    const int maxBitsX, const int maxBitsY);

// Extrapolate R to from the outermost PS module to the given z position using the given tgTheta.
// Use a second order approximation with only the first two terms.
bigInt extrapolateRSecondOrderFirstTwoTermsOnlyEmulator(const bigInt & R, const bigInt & z, const int layer, const bigInt & tgTheta,
                                                        const bigInt & chargeOverTwoRho, const std::vector<int> & uniqueLayers,
                                                        const std::vector<bigInt> & originalRInt,
                                                        const std::vector<bigInt> & originalZInt,
                                                        const bigInt & Rcut,
                                                        const int registerBits, const int bitsVar, const double & scaleZ,
                                                        const int reductionPhi, const int reductionZ,
                                                        const int reductionPt, const int reductionCotTheta,
                                                        const int maxBitsX, const int maxBitsY);

// Extrapolate R to from the outermost PS module to the given z position using the given tgTheta.
// Use a second order approximation.
bigInt extrapolateRSecondOrderEmulator(const bigInt & R, const bigInt & z, const int layer, const bigInt & tgTheta,
                                       const bigInt & chargeOverTwoRho, const std::vector<int> & uniqueLayers,
                                       const std::vector<bigInt> & originalRInt,
                                       const std::vector<bigInt> & originalZInt,
                                       const bigInt & Rcut,
                                       const int registerBits, const int bitsVar, const double & scaleZ,
                                       const int reductionPhi, const int reductionZ,
                                       const int reductionPt, const int reductionCotTheta,
                                       const int maxBitsX, const int maxBitsY);

class CorrectPhiForNonRadialStripsEmulator
{
 public:
  CorrectPhiForNonRadialStripsEmulator();

//  bigInt correctPhiForNonRadialStrips(const bigInt & phi, const bigInt & stripIndex,
//                                      const bigInt & extrapolatedR, const bigInt & R,
//                                      const bigInt & RCut, const bigInt & oneOverR,
//                                      const int layer, const int registerBits, const int bitsVar,
//                                      const int maxBitsX, const int maxBitsY,
//                                      const int reductionZ, const int reductionPhi);
  bigInt correctPhiForNonRadialStrips(const bigInt & phi, const bigInt & stripIndex,
                                      const bigInt & extrapolatedR, const bigInt & R,
                                      const bigInt & RCut, // const bigInt & oneOverR, const bigInt & oneOverRSquared,
                                      const int layer, const int registerBits, const int bitsVar,
                                      const int maxBitsX, const int maxBitsY,
                                      const int reductionZ, const int reductionPhi);
 private:
  bigInt stripPitch_;
  bigInt stripMean_;
  std::map<bigInt, bigInt> lookupOneOverRSquared_;
};


int computeShiftBits(const int bitsX, const int reducedAccuracyBitsX);


std::string convertToBitString(const bigInt & value, const int bits);


#endif //REMOTEPROJECTS_EMULATORTOOLS_H_H
