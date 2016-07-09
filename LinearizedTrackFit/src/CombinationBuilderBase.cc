//
// Created by Marco De Mattia on 7/5/16.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/CombinationBuilderBase.h"


CombinationBuilderBase::CombinationBuilderBase(const bool fillEmptyStubs) :
  road_(nullptr), totalCombinations_(1), fillEmptyStubs_(fillEmptyStubs)
{}


void CombinationBuilderBase::write(std::ofstream & outputFile, const StubsCombination & stubsCombination,
                                   const Road & road, const int bitsX, const int reducedAccuracyBitsX,
                                   const double & deltaPhi, const double & deltaR, const double & deltaZ)
{
//   int shiftBitsX = computeShiftBits(bitsX, reducedAccuracyBitsX);
//   double rotationFactor = computeRotationFactor(std::vector<double>(1, stubsCombination.phi(0)));
//   int stubsCombinationIndex = 0;
//   for (size_t i = 0; i < 6; ++i) {
//     outputFile << "layer " << i << std::endl;
//     if (road.stubsNum(i) == 0) {
//       for (int j=0; j<3; ++j) outputFile << 0 << std::endl;
//       continue;
//     }
//     else {
//       bigInt tempVarPhiInt = encode(stubsCombination.phi(stubsCombinationIndex) - rotationFactor, deltaPhi, bitsX);
//       if (shiftBitsX != 0) tempVarPhiInt = (tempVarPhiInt >> shiftBitsX);
//       outputFile << tempVarPhiInt << std::endl;
//       bigInt tempVarRInt = encode(stubsCombination.R(stubsCombinationIndex), deltaR, bitsX);
//       if (shiftBitsX != 0) tempVarRInt = (tempVarRInt >> shiftBitsX);
//       outputFile << tempVarRInt << std::endl;
//       bigInt tempVarZInt = encode(stubsCombination.z(stubsCombinationIndex), deltaZ, bitsX);
//       if (shiftBitsX != 0) tempVarZInt = (tempVarZInt >> shiftBitsX);
//       outputFile << tempVarZInt << std::endl;
//       ++stubsCombinationIndex;
//     }
//   }
}


void CombinationBuilderBase::fillGenInfo(StubsCombination & stubsCombination)
{
//   bool first = true;
//   int tpId = -1;
//   for (auto stub : stubsCombination) {
//     if (first) {
//       tpId = road_->tpId(stub.stubRef());
//       first = false;
//     }
//     else if (tpId != road_->tpId(stub.stubRef())) {
//       tpId = -1;
//       break;
//     }
//   }
//   if (tpId != -1) {
//     unsigned int stubRef = stubsCombination.begin()->stubRef();
//     stubsCombination.setGenTrack(road_->genChargeOverPt(stubRef), road_->genPhi0(stubRef), road_->genD0(stubRef),
//                                  road_->genCotTheta(stubRef), road_->genZ0(stubRef));
//   }
}
