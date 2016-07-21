//
// Created by Marco De Mattia on 2/28/16.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/RoadsTree.h"

RoadsTree::RoadsTree(const TString & fileName)
{
  L1TT = new TChain("ntupler/tree");
  L1TT->Add(fileName);

  L1TT->SetBranchAddress("TTStubs_phi",            &TTStubs_phi);
  L1TT->SetBranchAddress("TTStubs_r",              &TTStubs_r);
  L1TT->SetBranchAddress("TTStubs_z",              &TTStubs_z);

  L1TT->SetBranchAddress("trkParts_eta",           &trkParts_eta);
  L1TT->SetBranchAddress("trkParts_phi",           &trkParts_phi);
  L1TT->SetBranchAddress("trkParts_pt",            &trkParts_pt);
  L1TT->SetBranchAddress("trkParts_vx",            &trkParts_vx);
  L1TT->SetBranchAddress("trkParts_vy",            &trkParts_vy);
  L1TT->SetBranchAddress("trkParts_vz",            &trkParts_vz);
  L1TT->SetBranchAddress("TTStubs_tpId",           &TTStubs_tpId);
  L1TT->SetBranchAddress("trkParts_charge",        &trkParts_charge);
  L1TT->SetBranchAddress("trkParts_pdgId",         &trkParts_pdgId);
  L1TT->SetBranchAddress("TTStubs_modId",          &TTStubs_modId);
  L1TT->SetBranchAddress("AMTTRoads_patternRef",   &AMTTRoads_patternRef);
  L1TT->SetBranchAddress("AMTTRoads_tower",        &AMTTRoads_tower);
  L1TT->SetBranchAddress("AMTTRoads_nstubs",       &AMTTRoads_nstubs);
  L1TT->SetBranchAddress("AMTTRoads_patternInvPt", &AMTTRoads_patternInvPt);
  L1TT->SetBranchAddress("AMTTRoads_superstripId", &AMTTRoads_superstripIds);
  L1TT->SetBranchAddress("AMTTRoads_stubRefs",     &AMTTRoads_stubRefs);

  n_entries = L1TT->GetEntries();
}


void RoadsTree::printInfo()
{
  // std::cout << "where " << m_stub << " stub(s) where produced" << std::endl;
}
