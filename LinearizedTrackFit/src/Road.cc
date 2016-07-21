//
// Created by Marco De Mattia on 2/29/16.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/Road.h"

//// This constructor is to integrate with Seb's framework
//// -----------------------------------------------------
//Road::Road(const std::vector<Hit*> & stubs) :
//    maxStubsPerLayer_(4)
//{
//  // Loop over all stubs, arrange them by layer and store their information, including the index.
//  std::map<int, std::vector<Hit*> > stubMap;
//  for (Hit * s : stubs) {
//    int layer = (int)s->getLayer();
//    if (stubMap.find(layer) == stubMap.end()) {
//      stubMap.insert(std::make_pair(layer, std::vector<Hit*>()));
//    }
//    stubMap[layer].push_back(s);
//  }
//
//  // // If this road covers more than 6 layers take the innermost 6
//  int num = 0;
//  for (auto s : stubMap) {
//    if (num >= 6) {
//      std::cout << "road with more than 6 layers: ";
//      for (auto sss : stubMap) std::cout << sss.first << " ";
//      std::cout << std::endl;
//      // std::cout << "road with more than 6 layers. Stopping at layer " << s.first << std::endl;
//      break;
//    }
//    std::vector<Stub> stubs;
//    int maxStubIndex = int(s.second.size())-1;
//    for (int i = maxStubIndex; i >= 0; --i) {
//      // Only consider the last 4 stubs
//      if (maxStubIndex-i >= maxStubsPerLayer_) break;
//      // Loop over the stubs in this layer
//      Hit * h = s.second.at(i);
//      double phi = std::atan2(h->getY(), h->getX());
//      double R = std::sqrt(h->getX()*h->getX() + h->getY()*h->getY());
//      stubs.push_back(Stub(phi, R, h->getZ(), s.first, h->getStripNumber(), h->getID()));
//    }
//    stubs_.push_back(stubs);
//    // Only take the innermost 6 layers
//    ++num;
//  }
//
//  // If we inserted only 4 or 5 vectors add one more since the CB always expects 6 (some may be empty)
//  for (size_t i=stubs_.size(); i<6; ++i) {
//    stubs_.push_back(std::vector<Stub>());
//  }
//}


// This constructor is to integrate with Sergo's framework
// -------------------------------------------------------
Road::Road(const std::vector<std::vector<unsigned> > & stubRefs) :
    maxStubsPerLayer_(4)
{
  for (const auto & layerStubs : stubRefs) {
    std::vector<Stub> stubs;
    int maxStubIndex = int(layerStubs.size())-1;
    for (int i = maxStubIndex; i >= 0; --i) {
      // Only consider the last 4 stubs
      if (maxStubIndex-i >= maxStubsPerLayer_) break;
      // Loop over the stubs in this layer. For this constructor we only need to store the stubRef.
      stubs.push_back(Stub(0., 0., 0., 0, 0, layerStubs.at(i)));
    }
    stubs_.push_back(stubs);
  }
}


// This constructor is to integrate with the standalone fitter framework
// ---------------------------------------------------------------------
Road::Road(const std::shared_ptr<RoadsTree> & tree, const size_t & roadIndex) :
    tree_(tree), maxStubsPerLayer_(4)
{
  for (size_t layer = 0; layer < tree_->AMTTRoads_stubRefs->at(roadIndex).size(); ++layer) {
    std::vector<Stub> stubs;
    int maxStubIndex = int(tree_->AMTTRoads_stubRefs->at(roadIndex).at(layer).size())-1;
    for (int i = maxStubIndex; i >= 0; --i) {
      // Only consider the last 4 stubs
      if (maxStubIndex-i >= maxStubsPerLayer_) break;
      // Loop over the stubs in this layer
      unsigned int stubRef = tree_->AMTTRoads_stubRefs->at(roadIndex).at(layer).at(i);
      stubs.push_back(Stub(tree_->TTStubs_phi->at(stubRef), tree_->TTStubs_r->at(stubRef),
                           tree_->TTStubs_z->at(stubRef), layer+5, 0, stubRef));
    }
    stubs_.push_back(stubs);
  }
}


Road::Road(const std::shared_ptr<RoadsTree> & tree, const std::map<unsigned long, unsigned long> & stubsIndexes) :
    tree_(tree)
{
  stubs_ = std::vector<std::vector<Stub> >(6, std::vector<Stub>());
  for (auto it : stubsIndexes) {
    unsigned long stubRef = it.second;
    stubs_.at(it.first).push_back(Stub(tree_->TTStubs_phi->at(stubRef), tree_->TTStubs_r->at(stubRef),
                                       tree_->TTStubs_z->at(stubRef), it.first+5, 0, stubRef));
  }
}


int Road::pdgId(const unsigned int stubRef) const
{
  int tpId = tree_->TTStubs_tpId->at(stubRef);
  return tpId == -1 ? 0 : tree_->trkParts_pdgId->at(tpId);
}


double Road::genChargeOverPt(const unsigned int stubRef) const
{
  int tpId = tree_->TTStubs_tpId->at(stubRef);
  return tpId == -1 ? 0. : (tree_->trkParts_pt->at(tpId) != 0. ? tree_->trkParts_charge->at(tpId) /
                                                                 tree_->trkParts_pt->at(tpId) : 0.);
}


double Road::genPhi0(const unsigned int stubRef) const
{
  int tpId = tree_->TTStubs_tpId->at(stubRef);
  return tpId == -1 ? 0. : tree_->trkParts_phi->at(tpId);
}


double Road::genD0(const unsigned int stubRef) const
{
  int tpId = tree_->TTStubs_tpId->at(stubRef);
  return tpId == -1 ? 0. : 0.;
}


double Road::genCotTheta(const unsigned int stubRef) const
{
  int tpId = tree_->TTStubs_tpId->at(stubRef);
  return tpId == -1 ? 0. : 1./tan(2*atan(exp(-tree_->trkParts_eta->at(tpId))));
}


double Road::genZ0(const unsigned int stubRef) const
{
  int tpId = tree_->TTStubs_tpId->at(stubRef);
  return tpId == -1 ? 0. : tree_->trkParts_vz->at(tpId);
}


//void Road::write(std::vector<std::ofstream> & outputFile, const int bitsX, const int reducedAccuracyBitsX,
//                 const double & deltaPhi, const double & deltaR, const double & deltaZ) const
//{
//  assert(reducedAccuracyBitsX == 18);
//  double rotationFactor = 0.;
//  int shiftBitsX = computeShiftBits(bitsX, reducedAccuracyBitsX);
//  if (stubs_.at(0).size() > 0) rotationFactor = computeRotationFactor(std::vector<double>(1, stubs_.at(0).at(0).phi()));
//  else if (stubs_.at(1).size() > 0) rotationFactor = computeRotationFactor(std::vector<double>(1, stubs_.at(1).at(0).phi()));
//  else {
//    std::cout << "Road::write error: missing stubs in both the first and second layer" << std::endl;
//    throw;
//  }
//  // std::vector<std::string> output(6, std::string(240, '0'));
//  size_t i = 0;
//  for (auto layer = stubs_.begin(); layer != stubs_.end(); ++layer) {
//    std::string output(240, '0');
//    size_t s = 0;
//    for (auto stub = layer->begin(); stub != layer->end(); ++stub) {
//
//      // Encode the stubs with the same settings used in the emulator.
//      output.replace(s, 1, "1");
//
//      bigInt tempVarPhiInt = encode(stub->phi()-rotationFactor, deltaPhi, bitsX);
//      if (shiftBitsX != 0) tempVarPhiInt = (tempVarPhiInt >> shiftBitsX);
//      bigInt tempVarRInt = encode(stub->R(), deltaR, bitsX);
//      if (shiftBitsX != 0) tempVarRInt = (tempVarRInt >> shiftBitsX);
//      bigInt tempVarZInt = encode(stub->z(), deltaZ, bitsX);
//      if (shiftBitsX != 0) tempVarZInt = (tempVarZInt >> shiftBitsX);
//
//      int signPhi = tempVarPhiInt > 0 ? 1 : -1;
//      int signR = tempVarRInt > 0 ? 1 : -1;
//      int signZ = tempVarZInt > 0 ? 1 : -1;
//
//      // Note that this format is likely to change in the future.
//      boost::dynamic_bitset<> bitsPhi(reducedAccuracyBitsX, abs(int(tempVarPhiInt)));
//      boost::dynamic_bitset<> bitsR(reducedAccuracyBitsX, abs(int(tempVarRInt)));
//      boost::dynamic_bitset<> bitsZ(reducedAccuracyBitsX, abs(int(tempVarZInt)));
//
//      // Compute the negative under two's complement notation if the original sign is negative
//      // Note that this format is likely to change in the future.
//      if (signPhi == -1) bitsPhi = boost::dynamic_bitset<>(reducedAccuracyBitsX, bitsPhi.flip().to_ulong() + 1);
//      if (signR == -1) bitsR = boost::dynamic_bitset<>(reducedAccuracyBitsX, bitsR.flip().to_ulong() + 1);
//      if (signZ == -1) bitsZ = boost::dynamic_bitset<>(reducedAccuracyBitsX, bitsZ.flip().to_ulong() + 1);
//
//      // The input is still in bitsX bits while the output is truncated to the 18 bits in the L2G output format.
//      // Note that this format is likely to change in the future.
//      std::string stringBitsPhi;
//      boost::to_string(bitsPhi, stringBitsPhi);
//      std::string stringBitsR;
//      boost::to_string(bitsR, stringBitsR);
//      std::string stringBitsZ;
//      boost::to_string(bitsZ, stringBitsZ);
//      output.replace(s+1, reducedAccuracyBitsX, stringBitsPhi);
//      output.replace(s+19, reducedAccuracyBitsX, stringBitsR);
//      output.replace(s+37, reducedAccuracyBitsX, stringBitsZ);
//      s += 60;
//    }
////    std::cout << "output["<<i<<"] = " << output << std::endl;
//    outputFile[i] << output << std::endl;
//    ++i;
//  }
//}


void Road::write(std::vector<std::ofstream> & outputFile, const int bitsX, const int reducedAccuracyBitsPhi,
                 const int reducedAccuracyBitsR, const int reducedAccuracyBitsZ,
                 const double & deltaPhi, const double & deltaR, const double & deltaZ) const
{
  // assert(reducedAccuracyBitsX == 18);
  double rotationFactor = 0.;
  int shiftBitsPhi = computeShiftBits(bitsX, reducedAccuracyBitsPhi);
  int shiftBitsR = computeShiftBits(bitsX, reducedAccuracyBitsR);
  int shiftBitsZ = computeShiftBits(bitsX, reducedAccuracyBitsZ);
  if (stubs_.at(0).size() > 0) rotationFactor = computeRotationFactor(std::vector<double>(1, stubs_.at(0).at(0).phi()));
  else if (stubs_.at(1).size() > 0) rotationFactor = computeRotationFactor(std::vector<double>(1, stubs_.at(1).at(0).phi()));
  else {
    std::cout << "Road::write error: missing stubs in both the first and second layer" << std::endl;
    throw;
  }
  size_t i = 0;
  for (auto layer = stubs_.begin(); layer != stubs_.end(); ++layer) {
    std::string output(6*(1+reducedAccuracyBitsPhi+reducedAccuracyBitsR+reducedAccuracyBitsZ+5), '0');
    size_t s = 0;
    for (auto stub = layer->begin(); stub != layer->end(); ++stub) {

      // Encode the stubs with the same settings used in the emulator.
      output.replace(s, 1, "1");

      bigInt tempVarPhiInt = encode(stub->phi()-rotationFactor, deltaPhi, bitsX);
      if (shiftBitsPhi != 0) tempVarPhiInt = (tempVarPhiInt >> shiftBitsPhi);
      bigInt tempVarRInt = encode(stub->R(), deltaR, bitsX);
      if (shiftBitsR != 0) tempVarRInt = (tempVarRInt >> shiftBitsR);
      bigInt tempVarZInt = encode(stub->z(), deltaZ, bitsX);
      if (shiftBitsZ != 0) tempVarZInt = (tempVarZInt >> shiftBitsZ);

      std::string stringBitsPhi(convertToBitString(tempVarPhiInt, reducedAccuracyBitsPhi));
      std::string stringBitsR(convertToBitString(tempVarRInt, reducedAccuracyBitsR));
      std::string stringBitsZ(convertToBitString(tempVarZInt, reducedAccuracyBitsZ));
      output.replace(s+1, reducedAccuracyBitsPhi, stringBitsPhi);
      output.replace(s+1+reducedAccuracyBitsPhi, reducedAccuracyBitsR, stringBitsR);
      output.replace(s+1+reducedAccuracyBitsPhi+reducedAccuracyBitsR, reducedAccuracyBitsZ, stringBitsZ);
      s += 1+reducedAccuracyBitsPhi+reducedAccuracyBitsR+reducedAccuracyBitsZ+5;
    }
//    std::cout << "output["<<i<<"] = " << output << std::endl;
    outputFile[i] << output << std::endl;
    ++i;
  }
}