#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <EVENT/ReconstructedParticle.h>

class JetCombination;

typedef std::pair<fastjet::PseudoJet,fastjet::PseudoJet> JetPair;
typedef std::vector<JetPair> JetPairList;
typedef std::pair<JetPair,JetPair> PairOfJetPair;
typedef std::vector<PairOfJetPair> PairOfJetPairList;

class JetCombination{
 public:
  JetCombination(std::vector<fastjet::PseudoJet> &jets, float E, bool flag);
  ~JetCombination();
  void CreateAllPotentialJetPairs();
  void CombineJetPairs();
  void FindBestJetCombiation();

  std::vector<fastjet::PseudoJet> &getFinalJets(){return _outputJets;}
  float getChi2(){return _bestChi2;}
 protected:
  std::vector<fastjet::PseudoJet> _inputJets;
  bool SAME_PARTICLE_MASS;
  float _cmEnergy;

  JetPairList _jetPairList;
  PairOfJetPairList _pairOfJetPairList;

  PairOfJetPair _bestPairOfPair;
  float _bestChi2;
  float _massDiff;
  float _energyDiff;

  //after jet combination it's easier to keep particles inside pseudo jet
  std::vector<fastjet::PseudoJet> _outputJets;

};

bool operator==(fastjet::PseudoJet jet1,fastjet::PseudoJet jet2)
{
  return ( jet1.e()-jet2.e()<std::numeric_limits<float>::epsilon() &&
	   jet1.px()-jet2.px()<std::numeric_limits<float>::epsilon() &&
	   jet1.py()-jet2.py()<std::numeric_limits<float>::epsilon() &&
	   jet1.pz()-jet2.pz()<std::numeric_limits<float>::epsilon() )  ? true : false ;
 }
