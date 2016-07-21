//
// Created by Marco De Mattia on 9/7/15.
//

#ifndef REMOTEPROJECTS_MATRIXREADEREMULATOR_H
#define REMOTEPROJECTS_MATRIXREADEREMULATOR_H

#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixReader.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/EmulatorTools.h"
// #include <boost/multiprecision/cpp_int.hpp>

class MatrixReaderEmulator : public MatrixReader
{
 public:
  MatrixReaderEmulator(const std::string & inputFileName, const double & deltaX,
                       const int registerBits, const int bitsX, const int bitsA,
                       const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                       const MaxDeltaAndFactors & maxDeltaAndFactors,
                       const std::vector<bool> & powerTwoRanges);
  double normChi2(const std::vector<bigInt> & vars, const double & deltaX,
                  const int registerBits, const int bitsX, const int bitsA,
                  std::ofstream & emulatorOutput, const bool saveOutput) const;
  std::vector<double> trackParameters(const std::vector<bigInt> & vars, const double & deltaX,
                                      const int registerBits, const int bitsX, const int bitsA,
                                      std::ofstream & emulatorOutput, const bool saveOutput) const;
  std::vector<double> principalComponents(const std::vector<bigInt> & vars) const;
  std::vector<double> normalizedPrincipalComponents(const std::vector<bigInt> & vars,
                                                    const double & deltaX, const int registerBits,
                                                    const int bitsX, const int bitsA,
                                                    std::ofstream & emulatorOutput, const bool saveOutput) const;
  size_t getDeltaVSize() const { return deltaV_.size(); }
  size_t getDeltaDSize() const { return deltaD_.size(); }
  double getDeltaV(const int i) const { return deltaV_.at(i); }
  double getDeltaD(const int i) const { return deltaD_.at(i); }
  int getFactorV(const int i) const { return factorsV_.at(i); }
  int getFactorD(const int i) const { return factorsD_.at(i); }
  std::vector<bigInt> getChiTermsInt() const { return chiTermsInt_; }
  std::vector<bigInt> getParametersInt() const { return parametersInt_; }

  void write() const;

 private:
  // typedef boost::multiprecision::int128_t bigInt;
  // typedef long long int bigInt;
  // typedef int64_t bigInt;
  std::vector<bigInt> meanValuesInt_;
  std::vector<bigInt> meanParsInt_;
  std::vector<std::vector<bigInt> > VInt_;
  std::vector<std::vector<bigInt> > DInt_;
  std::vector<double> deltaV_;
  std::vector<int> factorsV_;
  std::vector<double> deltaD_;
  std::vector<int> factorsD_;
  const int maxBitsMultiplyUnitX_;
  const int maxBitsMultiplyUnitY_;
  std::string inputFileName_;
  int bitsA_;
  mutable std::vector<bigInt> chiTermsInt_;
  mutable std::vector<bigInt> parametersInt_;
};

#endif //REMOTEPROJECTS_MATRIXREADEREMULATOR_H
