#include "JetClustering.h"

JetClustering::JetClustering(fastjet::JetAlgorithm jetAlgo=fastjet::kt_algorithm, 
			     double par=0.7, 
			     fastjet::RecombinationScheme recomb_scheme=fastjet::E_scheme, 
			     fastjet::Strategy strategy=fastjet::Best)
{
  _jetDefinition = fastjet::JetDefinition(jetAlgo,par,recomb_scheme,strategy);
}

JetClustering::JetClustering(fastjet::JetAlgorithm jetAlgo=fastjet::ee_kt_algorithm,
			     fastjet::RecombinationScheme recomb_scheme=fastjet::E_scheme,
			     fastjet::Strategy strategy=fastjet::Best)
{
  _jetDefinition = fastjet::JetDefinition(jetAlgo,recomb_scheme,strategy);
}

JetClustering::JetClustering(fastjet::JetAlgorithm jetAlgo=fastjet::genkt_algorithm,
			     double par=0.7, double extra_param=0.0,
			     fastjet::RecombinationScheme recomb_scheme=fastjet::E_scheme, 
			     fastjet::Strategy strategy=fastjet::Best)
{
  _jetDefinition = fastjet::JetDefinition(jetAlgo,par,extra_param,recomb_scheme,strategy);
}

JetClustering::~JetClustering()
{
  _inputJets.clear();
  _reconstructedJets.clear();
}

void JetClustering::Init(std::vector<EVENT::ReconstructedParticle*> &vec)
{
  for(std::vector<EVENT::ReconstructedParticle*>::iterator it=vec.begin(); it!=vec.end(); ++it){
    fastjet::PseudoJet jet( (*it)->getMomentum()[0],
			    (*it)->getMomentum()[1],
			    (*it)->getMomentum()[2],
			    (*it)->getEnergy() );
    _inputJets.push_back(jet);
  }
}

void JetClustering::Clustering(int nJet)
{
  fastjet::ClusterSequence clusterSequence(_inputJets,_jetDefinition);
  //_reconstructedJets=clusterSequence.inclusive_jets();
  _reconstructedJets=clusterSequence.exclusive_jets(nJet);
  //std::cout << "rparam after jet algo = " << _jetDefinition.R() << std::endl;
  _ycutN_Nplus1=(float)clusterSequence.exclusive_ymerge(nJet);
  _ycutNminus1_N=(float)clusterSequence.exclusive_ymerge(nJet-1);
  _dcutN_Nplus1=(float)clusterSequence.exclusive_dmerge(nJet);
  _dcutNminus1_N=(float)clusterSequence.exclusive_dmerge(nJet-1);
}
