//
// Created by Marco De Mattia on 2/28/16.
//

#ifndef LINEARIZEDTRACKFIT_ROAD_H
#define LINEARIZEDTRACKFIT_ROAD_H

#include "LinearizedTrackFit/LinearizedTrackFit/interface/RoadsTree.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/Stub.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GetVariables.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/EmulatorTools.h"
// #include <boost/multiprecision/cpp_int.hpp>
// #include <bitset>
// #include <boost/dynamic_bitset.hpp>

// typedef boost::multiprecision::int128_t bigInt;
// typedef long long int bigInt;
typedef int64_t bigInt;

/**
 * A Road contains a vector of 6 elements, one per layer. Each vector contains the stubs in that layer.
 * For simplicity of access and scalability in case we need more variables, we store a pointer to the tree
 * and provide methods to access the different vectors. A check is made in the RoadsTreeReader to ensure
 * that there are exactly 6 vectors of stubs.
 */

class Road
{
 public:
  // Road(const std::vector<Hit*> & stubs);
  Road(const std::shared_ptr<RoadsTree> & tree, const size_t & roadIndex);
  Road(const std::shared_ptr<RoadsTree> & tree, const std::map<unsigned long, unsigned long> & stubsIndexes);
  Road(const std::vector<std::vector<unsigned> > & stubRefs);

  std::vector<std::vector<Stub> >::const_iterator beginLayer() const { return stubs_.begin(); }
  std::vector<std::vector<Stub> >::const_iterator endLayer() const { return stubs_.end(); }
  size_t stubsNum(const size_t & layer) const { return stubs_.at(layer).size(); }
  Stub getStub(const size_t & layer, const int count) const { return stubs_.at(layer).at(count); }
  int pdgId(const unsigned int stubRef) const;
  int tpId(const unsigned int stubRef) const { return tree_->TTStubs_tpId->at(stubRef); }
  double genChargeOverPt(const unsigned int stubRef) const;
  double genPhi0(const unsigned int stubRef) const;
  double genD0(const unsigned int stubRef) const;
  double genCotTheta(const unsigned int stubRef) const;
  double genZ0(const unsigned int stubRef) const;
  void write(std::vector<std::ofstream> & outputFile, const int bitsX, const int reducedAccuracyBitsPhi,
             const int reducedAccuracyBitsR, const int reducedAccuracyBitsZ,
             const double & deltaPhi, const double & deltaR, const double & deltaZ) const;

 private:
  std::shared_ptr<RoadsTree> tree_;
  int maxStubsPerLayer_;
  std::vector<std::vector<Stub> > stubs_;
};

#endif //LINEARIZEDTRACKFIT_ROAD_H
