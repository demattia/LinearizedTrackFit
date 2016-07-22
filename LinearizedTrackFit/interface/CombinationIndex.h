//
// Created by Marco De Mattia on 7/11/15.
//

#ifndef REMOTEPROJECTS_COMBINATIONINDEX_H
#define REMOTEPROJECTS_COMBINATIONINDEX_H

#include <vector>
#include <bitset>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/Stub.h"
typedef int64_t bigInt;


inline void combinationIndex(const std::vector<int> & layers, std::bitset<32> & bits) { for (auto l : layers) bits.set(l, 1); }


unsigned long combinationIndex(const std::vector<int> & layers, const int region);


unsigned long combinationIndex(const std::vector<int> & layers, const std::vector<double> & radius, const int regionsNumber);


unsigned long combinationIndex(const std::vector<int> & layers, const std::vector<bigInt> & radius,
                               const int regionsNumber, const bigInt & radiusCutPS2S,
                               const bigInt & radiusCut2S1, const bigInt & radiusCut2S2);


unsigned long combinationIndex(const std::vector<Stub> & stubs, const int regionsNumber);


void allCombinationIndexes(const std::vector<int> & layers, const std::vector<double> & radius,
                           std::vector<unsigned long> & combinationIndexList, const bool fiveOutOfSix);


template <class T>
inline void setLayerRadiusBits(const int layer, const double & radius, T & bits, const int regionsNumber)
{
  if (layer > 10) {
    if (radius < 61.) {
      bits.set(layer + 5, 1);
    }
    // Split the 2S modules part of the disks
    else if (regionsNumber == 14 &&
             (((layer == 11 || layer == 12) && radius < 82.5) ||
              ((layer == 13) && radius < 82.5) ||
              ((layer == 14) && radius < 77.) ||
              ((layer == 15) && radius < 77.))) {
      bits.set(layer + 10, 1);
    }
//    else if (((layer == 11 || layer == 12) && radius < 82.5) ||
//             ((layer == 13) && radius < 77.) ||
//             ((layer == 14) && radius < 72.) ||
//             ((layer == 15) && radius < 77.)) {
//      bits.set(layer + 10, 1);
//    }
  }
}


template <class T>
inline void setLayerRadiusBits(const int layer, const bigInt & radius, T & bits, const int regionsNumber,
                               const bigInt & radiusCutPS2S, const bigInt & radiusCut2S1, const bigInt & radiusCut2S2)
{
  if (layer > 10) {
    if (radius < radiusCutPS2S) {
      bits.set(layer + 5, 1);
    }
    // Split the 2S modules part of the disks
    else if (regionsNumber == 14 &&
             (((layer == 11 || layer == 12) && radius < radiusCut2S1) ||
              ((layer == 13) && radius < radiusCut2S1) ||
              ((layer == 14) && radius < radiusCut2S2) ||
              ((layer == 15) && radius < radiusCut2S2))) {
      bits.set(layer + 10, 1);
    }
  }
}


double radiusRange(const int layer, const double & radius, const int regionsNumber);


double radiusRange(const int layer, const bigInt & radius, const int regionsNumber,
                   const bigInt & radiusCut, const bigInt & radiusCut1, const bigInt & radiusCut2);


bool combIndexIsBarrel(const int combIndex);


#endif //REMOTEPROJECTS_COMBINATIONINDEX_H
