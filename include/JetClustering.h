#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <EVENT/ReconstructedParticle.h>


class JetClustering{
 public:
  JetClustering(fastjet::JetAlgorithm jetAlgo, 
		double par, 
		fastjet::RecombinationScheme recomb_scheme, 
		fastjet::Strategy strategy);
  
  JetClustering(fastjet::JetAlgorithm jetAlgo,
		fastjet::RecombinationScheme recomb_scheme,
		fastjet::Strategy strategy);
  
  JetClustering(fastjet::JetAlgorithm jetAlgo,
		double par, double extra_param,
		fastjet::RecombinationScheme recomb_scheme, 
		fastjet::Strategy strategy);

  ~JetClustering();
  
  void Init(std::vector<EVENT::ReconstructedParticle*> &vec);
  void Init(std::vector<fastjet::PseudoJet> &vec){_inputJets=vec;}
  void Clustering(int nJet);

  inline std::vector<fastjet::PseudoJet> &getFinalJets(){return _reconstructedJets;}
  inline float getYcutValueNM1_N(){return _ycutNminus1_N;}
  inline float getYcutValueN_NP1(){return _ycutN_Nplus1;}
  inline float getDcutValueNM1_N(){return _dcutNminus1_N;}
  inline float getDcutValueN_NP1(){return _dcutN_Nplus1;}

 protected:
  std::vector<fastjet::PseudoJet> _inputJets;
  fastjet::JetDefinition _jetDefinition;

  std::vector<fastjet::PseudoJet> _reconstructedJets;
  float _ycutN_Nplus1;
  float _ycutNminus1_N;
  float _dcutN_Nplus1;
  float _dcutNminus1_N;
};
