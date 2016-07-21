//
// Created by Marco De Mattia on 2/28/16.
//

#ifndef LINEARIZEDTRACKFIT_ROADSTREE_H
#define LINEARIZEDTRACKFIT_ROADSTREE_H

#include <vector>
#include <iostream>
#include "TString.h"
#include "TChain.h"
#include <cmath>

class RoadsTree
{
 public:
  RoadsTree(const TString & fileName);

  // Load the values for the event
  void getEntry(const int evt)
  {
    L1TT->GetEntry(evt);
  }

  void printInfo();

  TChain *L1TT;
  int n_entries;
  int m_stub;
//  std::vector<bool>    *trkParts_intime;
//  std::vector<bool>    *trkParts_primary;
//  std::vector<bool>    *trkParts_signal;
//  std::vector<float>   *genParts_eta;
//  std::vector<float>   *genParts_phi;
//  std::vector<float>   *genParts_pt;
//  std::vector<float>   *genParts_vx;
//  std::vector<float>   *genParts_vy;
//  std::vector<float>   *genParts_vz;
//  std::vector<float>   *TTStubs_coordx;
//  std::vector<float>   *TTStubs_coordy;
//  std::vector<float>   *TTStubs_eta;
  std::vector<float>   *TTStubs_phi;
  std::vector<float>   *TTStubs_r;
//  std::vector<float>   *TTStubs_roughPt;
//  std::vector<float>   *TTStubs_trigBend;
//  std::vector<float>   *TTStubs_x;
//  std::vector<float>   *TTStubs_y;
  std::vector<float>   *TTStubs_z;
  std::vector<float>   *trkParts_eta;
  std::vector<float>   *trkParts_phi;
  std::vector<float>   *trkParts_pt;
  std::vector<float>   *trkParts_vx;
  std::vector<float>   *trkParts_vy;
  std::vector<float>   *trkParts_vz;
//  std::vector<int>     *genParts_charge;
  std::vector<int>     *TTStubs_tpId;
  std::vector<int>     *trkParts_charge;
  std::vector<int>     *trkParts_pdgId;
//  std::vector<unsigned int> *TTStubs_clusWidth0;
//  std::vector<unsigned int> *TTStubs_clusWidth1;
  std::vector<unsigned int> *TTStubs_modId;
  std::vector<unsigned int> *AMTTRoads_patternRef;
  std::vector<unsigned int> *AMTTRoads_tower;
  std::vector<unsigned int> *AMTTRoads_nstubs;
  std::vector<float>   *AMTTRoads_patternInvPt;
  std::vector<std::vector<unsigned int> > *AMTTRoads_superstripIds;
  std::vector<std::vector<std::vector<unsigned int> > > *AMTTRoads_stubRefs;
};

#endif //LINEARIZEDTRACKFIT_ROADSTREE_H
