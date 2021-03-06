//
// Created by Marco De Mattia on 2/29/16.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/SimpleCombinationBuilder.h"


SimpleCombinationBuilder::SimpleCombinationBuilder(const bool fillEmptyStubs) :
  CombinationBuilderBase(fillEmptyStubs)
{
  for (int layer = 0; layer < 6; ++layer) {
    combUnits_.push_back(CombUnit(0));
  }
}


void SimpleCombinationBuilder::initialize(const Road & road)
{
  // Set the ranges for the counters to this combination
  totalCombinations_ = 1;
  size_t layer = 0;
  for (auto l = road.beginLayer(); l != road.endLayer(); ++l) {
    int range = l->size() == 0 ? 0 : int(l->size())-1;
    combUnits_.at(layer).setRange(range);
    totalCombinations_ *= range+1;
    ++layer;
  }

  carryOver_ = std::vector<int>(7, 0);
  carryOver_[0] = 1;

  road_ = &road;
}


StubsCombination SimpleCombinationBuilder::nextCombination()
{
  StubsCombination stubsCombination;
  if (carryOver_[6] == 0) {
    for (size_t layer = 0; layer < 6; ++layer) {
      if (road_->stubsNum(layer) > 0) {
        stubsCombination.pushStub(road_->getStub(layer, combUnits_[layer].comb()));
      }
      else if (fillEmptyStubs_) {
        stubsCombination.pushStub(Stub(0., 0., 0., int(layer), 0, 999999999));
      }
    }
    for (size_t layer = 0; layer < 6; ++layer) {
      carryOver_[layer + 1] = combUnits_[layer].decrement(carryOver_[layer]);
    }
  }
  // Fill generator-level information
  fillGenInfo(stubsCombination);

  return stubsCombination;
}
