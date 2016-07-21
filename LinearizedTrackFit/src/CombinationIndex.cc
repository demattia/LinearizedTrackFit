//
// Created by Marco De Mattia on 7/11/15.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/CombinationIndex.h"
#include <iostream>
typedef int64_t bigInt;


unsigned long combinationIndex(const std::vector<int> & layers, const std::vector<double> & radius,
                               const int regionsNumber)
{
  std::bitset<32> bits;
  // Set the bits for the layers
  combinationIndex(layers, bits);

  // Set bits to determine the type of modules in the disks (2S vs PS) and upper or lower 2S radius.
  // Their positions in the bitset are the disk index + 5, since there
  // are 5 disks per side and they never appear together.
  for (unsigned int i=0; i<layers.size(); ++i) {
    setLayerRadiusBits(layers[i], radius[i], bits, regionsNumber);
  }
  return bits.to_ulong();
}


unsigned long combinationIndex(const std::vector<int> & layers, const std::vector<bigInt> & radius,
                               const int regionsNumber, const bigInt & radiusCutPS2S,
                               const bigInt & radiusCut2S1, const bigInt & radiusCut2S2)
{
  std::bitset<32> bits;
  // Set the bits for the layers
  combinationIndex(layers, bits);

  // Set bits to determine the type of modules in the disks (2S vs PS) and upper or lower 2S radius.
  // Their positions in the bitset are the disk index + 5, since there
  // are 5 disks per side and they never appear together.
  for (unsigned int i=0; i<layers.size(); ++i) {
    setLayerRadiusBits(layers[i], radius[i], bits, regionsNumber, radiusCutPS2S, radiusCut2S1, radiusCut2S2);
  }
  return bits.to_ulong();
}


unsigned long combinationIndex(const std::vector<Stub> & stubs, const int regionsNumber)
{
  std::bitset<32> bits;
  for (const Stub & stub : stubs) {
    int layer = stub.layer();
    bits.set(layer, 1);
    // Set bits to determine the type of modules in the disks (2S vs PS).
    // Their positions in the bitset are the disk index + 10, since there are 10 disks in total.
    setLayerRadiusBits(layer, stub.R(), bits, regionsNumber);
  }
  return bits.to_ulong();
}


double radiusRange(const int layer, const double & radius, const int regionsNumber)
{
  if (layer > 10) {
    if (radius < 61.) {
      return 0.;
    }
    // Split the 2S modules part of the disks
    else if (regionsNumber == 14 && (((layer == 11 || layer == 12 || layer == 13) && radius < 82.5) ||
                                     ((layer == 14 || layer == 15) && radius < 77.))) {
      return 70.;
    }
    return 110.;
  }
  return 0.;
}


double radiusRange(const int layer, const bigInt & radius, const int regionsNumber,
                   const bigInt & radiusCut, const bigInt & radiusCut1, const bigInt & radiusCut2)
{
  if (layer > 10) {
    if (radius < radiusCut) {
      return 0.;
    }
    // Split the 2S modules part of the disks
    else if (regionsNumber == 14 && (((layer == 11 || layer == 12 || layer == 13) && radius < radiusCut1) ||
                                     ((layer == 14 || layer == 15) && radius < radiusCut2))) {
      return 70.;
    }
    return 110.;
  }
  return 0.;
}
