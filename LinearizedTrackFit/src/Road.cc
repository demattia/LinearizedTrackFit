//
// Created by Marco De Mattia on 2/29/16.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/Road.h"

// // This constructor is to integrate with Seb's framework
// // -----------------------------------------------------
// Road::Road(const std::vector<Hit*> & stubs) :
//   maxStubsPerLayer_(4)
// {
//   // Loop over all stubs, arrange them by layer and store their information, including the index.
//   std::map<int, std::vector<Hit*> > stubMap;
//   for (Hit * s : stubs) {
//     int layer = (int)s->getLayer();
//     if (stubMap.find(layer) == stubMap.end()) {
//       stubMap.insert(std::make_pair(layer, std::vector<Hit*>()));
//     }
//     stubMap[layer].push_back(s);
//   }

//   // // If this road covers more than 6 layers take the innermost 6
//   int num = 0;
//   for (auto s : stubMap) {
//     if (num >= 6) {
//       std::cout << "road with more than 6 layers: ";
//       for (auto sss : stubMap) std::cout << sss.first << " ";
//       std::cout << std::endl;
//       // std::cout << "road with more than 6 layers. Stopping at layer " << s.first << std::endl;
//       break;
//     }
//     std::vector<Stub> stubs;
//     int maxStubIndex = int(s.second.size())-1;
//     for (int i = maxStubIndex; i >= 0; --i) {
//       // Only consider the last 4 stubs
//       if (maxStubIndex-i >= maxStubsPerLayer_) break;
//       // Loop over the stubs in this layer
//       Hit * h = s.second.at(i);
//       double phi = std::atan2(h->getY(), h->getX());
//       double R = std::sqrt(h->getX()*h->getX() + h->getY()*h->getY());
//       stubs.push_back(Stub(phi, R, h->getZ(), s.first, h->getStripNumber(), h->getID()));
//     }
//     stubs_.push_back(stubs);
//     // Only take the innermost 6 layers
//     ++num;
//   }

//   // If we inserted only 4 or 5 vectors add one more since the CB always expects 6 (some may be empty)
//   for (size_t i=stubs_.size(); i<6; ++i) {
//     stubs_.push_back(std::vector<Stub>());
//   }
// }


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


// // This constructor is to integrate with the standalone fitter framework
// // ---------------------------------------------------------------------
// Road::Road(const std::shared_ptr<RoadsTree> & tree, const size_t & roadIndex) :
//     tree_(tree), maxStubsPerLayer_(4)
// {
//   for (size_t layer = 0; layer < tree_->AMTTRoads_stubRefs->at(roadIndex).size(); ++layer) {
//     std::vector<Stub> stubs;
//     int maxStubIndex = int(tree_->AMTTRoads_stubRefs->at(roadIndex).at(layer).size())-1;
//     for (int i = maxStubIndex; i >= 0; --i) {
//       // Only consider the last 4 stubs
//       if (maxStubIndex-i >= maxStubsPerLayer_) break;
//       // Loop over the stubs in this layer
//       unsigned int stubRef = tree_->AMTTRoads_stubRefs->at(roadIndex).at(layer).at(i);
//       stubs.push_back(Stub(tree_->TTStubs_phi->at(stubRef), tree_->TTStubs_r->at(stubRef),
//                            tree_->TTStubs_z->at(stubRef), layer+5, 0, stubRef));
//     }
//     stubs_.push_back(stubs);
//   }
// }
// 
// 
// Road::Road(const std::shared_ptr<RoadsTree> & tree, const std::map<unsigned long, unsigned long> & stubsIndexes) :
//     tree_(tree)
// {
//   stubs_ = std::vector<std::vector<Stub> >(6, std::vector<Stub>());
//   for (auto it : stubsIndexes) {
//     unsigned long stubRef = it.second;
//     stubs_.at(it.first).push_back(Stub(tree_->TTStubs_phi->at(stubRef), tree_->TTStubs_r->at(stubRef),
//                                        tree_->TTStubs_z->at(stubRef), it.first+5, 0, stubRef));
//   }
// }
