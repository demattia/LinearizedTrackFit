//
// Created by Marco De Mattia on 9/7/15.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixReaderEmulator.h"

MatrixReaderEmulator::MatrixReaderEmulator(const std::string & inputFileName, const double & deltaX,
                                           const int registerBits, const int bitsX, const int bitsA,
                                           const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                                           const MaxDeltaAndFactors & maxDeltaAndFactors,
                                           const std::vector<bool> & powerTwoRanges) :
    MatrixReader(inputFileName), maxBitsMultiplyUnitX_(maxBitsMultiplyUnitX), maxBitsMultiplyUnitY_(maxBitsMultiplyUnitY),
    inputFileName_(inputFileName), bitsA_(bitsA)
{
  encodeVector(meanValues_, meanValuesInt_, deltaX, bitsX);
//  if (maxDeltaAndFactors.fixValues == true) {
  bool Vcoeff = true;
//  if (powerTwoRanges.size() == 0) {
//    std::cout << "Error: empty powerTwoRanges received" << std::endl;
//  }
  encodeMatrix(V_, VInt_, deltaX, deltaV_, factorsV_, bitsA, maxDeltaAndFactors.maxDeltaV, maxDeltaAndFactors.maxFactorV, Vcoeff, std::vector<bool>(V_.rows(), true));
  if (powerTwoRanges.size() == 0) {
    encodeMatrix(D_, DInt_, deltaX, deltaD_, factorsD_, bitsA, maxDeltaAndFactors.maxDeltaD,
                 maxDeltaAndFactors.maxFactorD, false, std::vector<bool>(V_.rows(), true));
  }
  else {
    encodeMatrix(D_, DInt_, deltaX, deltaD_, factorsD_, bitsA, maxDeltaAndFactors.maxDeltaD,
                 maxDeltaAndFactors.maxFactorD, false, powerTwoRanges);
  }
//  }
//  else {
//    encodeMatrix(V_, VInt_, deltaX, deltaV_, factorsV_, bitsA);
//    encodeMatrix(D_, DInt_, deltaX, deltaD_, factorsD_, bitsA);
//  }
  encodeVector(meanPars_, meanParsInt_, deltaX, deltaD_, bitsX+bitsA);
}


double MatrixReaderEmulator::normChi2(const std::vector<bigInt> & vars, const double & deltaX,
                                      const int registerBits, const int bitsX, const int bitsA,
                                      std::ofstream & emulatorOutput, const bool saveOutput) const
{
  // Compute the deltas
  std::vector<bigInt> deltaVars;
  preDiffer(vars, meanValuesInt_, deltaVars, bitsX);

  chiTermsInt_.clear();

  // Evaluate the chi2
  bigInt chi2Int = 0;
  // If the factorsV_ are too large, truncate so that the output fits in 48 bits
  std::vector<int>::const_iterator it = std::max_element(factorsV_.begin(), factorsV_.begin()+nDof_);
  int maxBits = bitsX+2*(*it)-registerBits;
  for (int i=0; i<nDof_; ++i) {
    bigInt chiTerm = macUnit(deltaVars, VInt_.at(i), 0, registerBits, bitsA,
                             maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_);
    bigInt tempChiTerm = macUnit(deltaVars, VInt_.at(i), 0, registerBits, 0,
                                 maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_);
    chiTermsInt_.push_back(tempChiTerm);
    if (saveOutput) emulatorOutput << "normChiTerm["<<i<<"] = " << tempChiTerm << std::endl;

    if (maxBits <= 0) maxBits = 0;
    chi2Int = adderUnit(chi2Int, multiplyUnit(chiTerm, chiTerm, registerBits, bitsX, 2 * factorsV_.at(i) - maxBits,
                                              maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_), registerBits);
//    chi2Int = adderUnit(chi2Int, multiplyUnit(tempChiTerm, tempChiTerm, registerBits, bitsX+bitsA, 2*factorsV_.at(i),
//                                              maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_), registerBits);
  }
  // The 2*factorsV_ above already shifts for deltaX*deltaA
  return decode(chi2Int, 1., bitsX-maxBits)/nDof_;
//  return decode(chi2Int, 1., bitsX+bitsA)/nDof_;
}


std::vector<double> MatrixReaderEmulator::trackParameters(const std::vector<bigInt> & vars, const double & deltaX,
                                                          const int registerBits, const int bitsX, const int bitsA,
                                                          std::ofstream & emulatorOutput, const bool saveOutput) const
{
  // Compute the deltas
  std::vector<bigInt> deltaVars;
  preDiffer(vars, meanValuesInt_, deltaVars, bitsX);

  parametersInt_.clear();

  std::vector<double> pars;
  // Estimate track parameters
  for (int i=0; i<nTrackParameters_; ++i) {
    // Saving the values before the truncation
    bigInt tempParameter = macUnit(deltaVars, DInt_.at(i), meanParsInt_.at(i), registerBits,
                                   0, maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_);
    parametersInt_.push_back(tempParameter);
    if (saveOutput) emulatorOutput << "parInt = " << tempParameter << std::endl;
    pars.push_back(decode(macUnit(deltaVars, DInt_.at(i), meanParsInt_.at(i), registerBits, bitsA,
                                  maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_), deltaX*deltaD_.at(i), bitsX));
//    pars.push_back(decode(macUnit(deltaVars, DInt_.at(i), meanParsInt_.at(i), registerBits, 0,
//                                  maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_), deltaX*deltaD_.at(i), bitsA+bitsX));
  }

  return pars;
}


std::vector<double> MatrixReaderEmulator::principalComponents(const std::vector<bigInt> & vars) const
{
  return std::vector<double>(nVars_, 0.);
}


std::vector<double> MatrixReaderEmulator::normalizedPrincipalComponents(const std::vector<bigInt> & vars,
                                                                        const double & deltaX, const int registerBits,
                                                                        const int bitsX, const int bitsA,
                                                                        std::ofstream & emulatorOutput,
                                                                        const bool saveOutput) const
{
  std::vector<double> nPC;
  // Compute the deltas
  std::vector<bigInt> deltaVars;
  preDiffer(vars, meanValuesInt_, deltaVars, bitsX);

  for (int i=0; i<nDof_; ++i) {
    bigInt chiTerm = macUnit(deltaVars, VInt_.at(i), 0, registerBits, bitsA,
                             maxBitsMultiplyUnitX_, maxBitsMultiplyUnitY_);
    nPC.push_back(decode(chiTerm, deltaX*deltaV_.at(i), bitsX));
  }
  return nPC;
}


void MatrixReaderEmulator::write() const
{
  std::string fileName = inputFileName_;
  std::size_t foundOne = fileName.rfind("/");
  if (foundOne != std::string::npos) {
    std::size_t found = fileName.rfind("/", foundOne-2);
    if (found != std::string::npos) {
      fileName = fileName.replace(foundOne-1, 2, "_").substr(found + 1);
    }
    else {
      std::cout << "MatrixReaderEmulator: Error building output file name from input file name = " << inputFileName_ << std::endl;
      throw;
    }
  }
  // std::cout << inputFileName_ << std::endl;
  std::ofstream outfile;
  outfile.open("Emulator_"+fileName);
  if(!outfile) {
    std::cout << "error opening Emulator_"+fileName << std::endl;
    return;
  }

  outfile << "Mean values:" << std::endl;
  for (const auto & meanValue : meanValuesInt_) {
    outfile << meanValue << std::endl;
  }

  outfile << "Mean pars:" << std::endl;
  for (const auto & meanPar : meanParsInt_) {
    outfile << meanPar << std::endl;
  }

  outfile << "V:" << std::endl;
  for (const auto & matrixV : VInt_) {
    for (const auto & V : matrixV) {
      outfile << V << std::endl;
    }
  }

  outfile << "D:" << std::endl;
  for (const auto & matrixD : DInt_) {
    for (const auto & D : matrixD) {
      outfile << D << std::endl;
    }
  }

  outfile << "deltaV:" << std::endl;
  for (const auto & deltaV : deltaV_) {
    outfile << deltaV << std::endl;
  }
  outfile << "factorsV:" << std::endl;
  for (const auto & factorsV : factorsV_) {
    // Write it already in a format usable by the firmware.
    outfile << (bitsA_ - factorsV) << std::endl;
  }

  outfile << "deltaD:" << std::endl;
  for (const auto & deltaD : deltaD_) {
    outfile << deltaD << std::endl;
  }
  outfile << "factorsD:" << std::endl;
  for (const auto & factorsD : factorsD_) {
    // Write it already in a format usable by the firmware.
    outfile << (bitsA_ - factorsD) << std::endl;
  }

  outfile << std::endl;
  outfile.close();
}