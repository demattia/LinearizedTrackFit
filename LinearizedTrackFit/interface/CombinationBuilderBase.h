//
// Created by Marco De Mattia on 7/5/16.
//

#ifndef LINEARIZEDTRACKFIT_COMBINATIONBUILDERBASE_H_H
#define LINEARIZEDTRACKFIT_COMBINATIONBUILDERBASE_H_H

#include "LinearizedTrackFit/LinearizedTrackFit/interface/StubsCombination.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/Road.h"

/**
 * Defines the interface of a combination builder and implements some of the common data members and methods.
 * The fillEmptyStubs argument in the constructor specifies if the combination builder should return also
 * empty stubs in the SubsCombination vector. This is used in the integration with Fermilab's version of the
 * AM simulation.
 */

class CombinationBuilderBase
{
 public:
  CombinationBuilderBase(const bool fillEmptyStubs = false);
  virtual void initialize(const Road & road) = 0;
  virtual StubsCombination nextCombination() = 0;
  int totalCombinations() const { return totalCombinations_; }
  void write(std::ofstream & outputFile, const StubsCombination & stubsCombination,
             const Road & road, const int bitsX, const int reducedAccuracyBitsX,
             const double & deltaPhi, const double & deltaR, const double & deltaZ);

 protected:
  void fillGenInfo(StubsCombination & stubsCombination);
  std::vector<int> carryOver_;
  const Road * road_;
  int totalCombinations_;
  bool fillEmptyStubs_;
};

#endif //LINEARIZEDTRACKFIT_COMBINATIONBUILDERBASE_H_H
