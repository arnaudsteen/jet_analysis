#include "JetCombination.h"

JetCombination::JetCombination(std::vector<fastjet::PseudoJet> &jets, float E, bool flag=true)
{
  _inputJets=jets;
  _cmEnergy=E;
  _jetPairList.clear();
  SAME_PARTICLE_MASS=flag;

  CreateAllPotentialJetPairs();
  CombineJetPairs();
  FindBestJetCombiation();
}

void JetCombination::CreateAllPotentialJetPairs()
{
  for(std::vector<fastjet::PseudoJet>::iterator it=_inputJets.begin(); it!=_inputJets.end(); ++it){
    for(std::vector<fastjet::PseudoJet>::iterator jt=it+1; jt!=_inputJets.end(); ++jt){
      JetPair jetPair( (*it),(*jt) );
      _jetPairList.push_back(jetPair);
    }
  }
}

void JetCombination::CombineJetPairs()
{
  for(JetPairList::iterator it=_jetPairList.begin(); it!=_jetPairList.end(); ++it){
    for(JetPairList::iterator jt=it+1; jt!=_jetPairList.end(); ++jt){
      if( (*it).first==(*jt).first  ||
	  (*it).first==(*jt).second ||
	  (*it).second==(*jt).first ||
	  (*it).second==(*jt).second )
	continue;
      else{
	PairOfJetPair pJP( (*it),(*jt) );
	_pairOfJetPairList.push_back(pJP);
      }
    }
  }
}

void JetCombination::FindBestJetCombiation()
{
  _bestChi2 = std::numeric_limits<float>::max();
  for( PairOfJetPairList::iterator it=_pairOfJetPairList.begin(); it!=_pairOfJetPairList.end(); ++it ){
    fastjet::PseudoJet pJ1= (*it).first.first + (*it).first.second;
    fastjet::PseudoJet pJ2= (*it).second.first + (*it).second.second;
    float chi2 = ( ( pJ1.e() + pJ2.e() - _cmEnergy)*( pJ1.e() + pJ2.e() - _cmEnergy) +
		   ( pJ1.m() - pJ2.m() )*( pJ1.m() - pJ2.m() ) );
    if( chi2 < _bestChi2 ){
      _bestPairOfPair=(*it);
      _bestChi2=chi2;
    }
  }
  _outputJets.push_back( _bestPairOfPair.first.first+_bestPairOfPair.first.second );
  _outputJets.push_back( _bestPairOfPair.second.first+_bestPairOfPair.second.second );
  _massDiff = fabs( (_bestPairOfPair.first.first+_bestPairOfPair.first.second).m() -
		    (_bestPairOfPair.second.first+_bestPairOfPair.second.second).m() );
  _energyDiff = fabs( (_bestPairOfPair.first.first+_bestPairOfPair.first.second).e() +
		      (_bestPairOfPair.second.first+_bestPairOfPair.second.second).e() -
		      _cmEnergy );
}

