#include "MyFastJetProcessor.h"
#include <iostream>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include "cmath"
#include <UTIL/LCTOOLS.h>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

MyFastJetProcessor aMyFastJetProcessor ;

MyFastJetProcessor::MyFastJetProcessor() : Processor("MyFastJetProcessor") {
  
  // Processor description
  _description = "FastJet Clustering ..." ;
 

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			   "InputCollection",
			   "Collection of reconstructed particles",
			   _inputCollection,
			   std::string("Unset") );
    
  registerProcessorParameter("Algorithm",
			     "FastJet algorithm",
			     sAlgorithm,
			     std::string("antikt_algorithm")); 

  registerProcessorParameter("R",
			     "R Parameter",
			     _RPar,
			     double(0.7)); 

  registerProcessorParameter("EjetMin",
			     "Ejet",
			     _eJet,
			     double(10.0)); 

  registerProcessorParameter("NJets",
			     "max nb of jets",
			     _nJetMax,
			     int(25)); 

}

void MyFastJetProcessor::init() { 
  
  printParameters() ;

  _rootfile = new TFile("MyFastJetProcessor.root","RECREATE");
  _tree = (TTree*)_rootfile->Get("tree");
  if(!_tree){
    std::cout << "tree creation" << std::endl; 
    _tree = new TTree("tree","Jet output variables");
  }
  
  _tree->Branch("nInputRecoParticles",&_nInputRecoParticles,"nInputRecoParticles/I");
  _tree->Branch("nHEJet",&_nHEJet,"nHEJet/I");
  _tree->Branch("jetEnergy","std::vector<double>",&_jetEnergy);
  _tree->Branch("jetMass","std::vector<double>",&_jetMass);
  _tree->Branch("jetPhi","std::vector<double>",&_jetPhi);
  _tree->Branch("jetEta","std::vector<double>",&_jetEta);
  _tree->Branch("jetTheta","std::vector<double>",&_jetTheta);
  _tree->Branch("jetPx","std::vector<double>",&_jetPx);
  _tree->Branch("jetPy","std::vector<double>",&_jetPy);
  _tree->Branch("jetPz","std::vector<double>",&_jetPz);
  _tree->Branch("y_nm1_n",&_y_NM1_N,"y_nm1_n/D");
  _tree->Branch("y_n_np1",&_y_N_NP1,"y_n_np1/D");
  _tree->Branch("d_nm1_n",&_d_NM1_N,"d_nm1_n/D");
  _tree->Branch("d_n_np1",&_d_N_NP1,"d_n_np1/D");
  _tree->Branch("W_mass","std::vector<double>",&_WMass);
    

  if(sAlgorithm=="kt_algorithm"){
    SetAlgorithm(fastjet::kt_algorithm);
  }else if(sAlgorithm=="cambridge_algorithm"){ 
    SetAlgorithm(fastjet::cambridge_algorithm);
  }else if(sAlgorithm=="antikt_algorithm"){ 
    SetAlgorithm(fastjet::antikt_algorithm);
  }else if(sAlgorithm=="ee_kt_algorithm"){ 
    SetAlgorithm(fastjet::ee_kt_algorithm);
  }else if(sAlgorithm=="genkt_algorithm"){ 
    SetAlgorithm(fastjet::genkt_algorithm);
  }else if(sAlgorithm=="ee_genkt_algorithm"){ 
    SetAlgorithm(fastjet::ee_genkt_algorithm);
  }
  else{
    streamlog_out( ERROR ) << sAlgorithm << " is not an algorithm available option" << std::endl;
    throw;
  }
}

void MyFastJetProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun = run->getRunNumber() ;

} 

void MyFastJetProcessor::WriteOutputJetInfo()
{
  std::vector<fastjet::PseudoJet> finalJet=_jetClustering->getFinalJets();
  std::sort(finalJet.begin(), finalJet.end(), MyFastJetProcessor::SortJetWithEnergy);
  for(std::vector<fastjet::PseudoJet>::iterator it=finalJet.begin(); it!=finalJet.end(); ++it){
    streamlog_out( DEBUG ) << (*it).e() << std::endl;
    if( (*it).e() > _eJet ){
      _nHEJet++;
      streamlog_out( DEBUG ) << "jet energy= " << (*it).e() << " jet mass= " << (*it).m() << "   p= " << sqrt( (*it).px()*(*it).px() + (*it).py()*(*it).py() +(*it).pz()*(*it).pz() )
			       << "   eta= " << (*it).eta() << "   phi= " << (*it).phi() << " " 
			       << "   p= (" << (*it).px() << ";" << (*it).py() << ";" << (*it).pz() << ")" << std::endl;
      _jetEnergy.push_back((*it).e());
      _jetMass.push_back((*it).m());
      _jetPx.push_back((*it).px());
      _jetPy.push_back((*it).py());
      _jetPz.push_back((*it).pz());
      _jetPhi.push_back((*it).phi());
      _jetEta.push_back((*it).eta());
      _jetTheta.push_back(2*atan(exp(-(*it).eta())));
    }
  }
  _nHEJet=finalJet.size();
  //try jet pairing
  if(finalJet.size()>=4){
    double ew1=finalJet.at(0).e() + finalJet.at(3).e();    
    double ew2=finalJet.at(1).e() + finalJet.at(2).e();
    double pw1=sqrt( (finalJet.at(0).px() + finalJet.at(3).px())*(finalJet.at(0).px() + finalJet.at(3).px()) +
		     (finalJet.at(0).py() + finalJet.at(3).py())*(finalJet.at(0).py() + finalJet.at(3).py()) +
		     (finalJet.at(0).pz() + finalJet.at(3).pz())*(finalJet.at(0).pz() + finalJet.at(3).pz()) );
    double pw2=sqrt( (finalJet.at(1).px() + finalJet.at(2).px())*(finalJet.at(1).px() + finalJet.at(2).px()) +
		     (finalJet.at(1).py() + finalJet.at(2).py())*(finalJet.at(1).py() + finalJet.at(2).py()) +
		     (finalJet.at(1).pz() + finalJet.at(2).pz())*(finalJet.at(1).pz() + finalJet.at(2).pz()) );
    _WMass.push_back(sqrt(ew1*ew1-pw1*pw1));
    _WMass.push_back(sqrt(ew2*ew2-pw2*pw2));
    //std::cout << "W1 mass = " << sqrt(ew1*ew1-pw1*pw1) << "\t W2 mass = " << sqrt(ew2*ew2-pw2*pw2) << std::endl;
  }
}

void MyFastJetProcessor::doClustering()
{
  if(GetAlgorithm()==fastjet::kt_algorithm){
    _jetClustering=new JetClustering(GetAlgorithm(),_RPar,fastjet::E_scheme,fastjet::Best);
  }else if(GetAlgorithm()==fastjet::cambridge_algorithm){ 
    _jetClustering=new JetClustering(GetAlgorithm(),_RPar,fastjet::E_scheme,fastjet::Best);
  }else if(GetAlgorithm()==fastjet::antikt_algorithm){ 
    _jetClustering=new JetClustering(GetAlgorithm(),_RPar,fastjet::E_scheme,fastjet::Best);
  }else if(GetAlgorithm()==fastjet::ee_kt_algorithm){ 
    _jetClustering=new JetClustering(GetAlgorithm(),fastjet::E_scheme,fastjet::Best);
  }else if(GetAlgorithm()==fastjet::genkt_algorithm){ 
    _jetClustering=new JetClustering(GetAlgorithm(),_RPar,_xtra_param,fastjet::E_scheme,fastjet::Best);
  }else if(GetAlgorithm()==fastjet::ee_genkt_algorithm){ 
    _jetClustering=new JetClustering(GetAlgorithm(),_RPar,_xtra_param,fastjet::E_scheme,fastjet::Best);
  }
  else{
    streamlog_out( ERROR ) << GetAlgorithm() << " is not an algorithm available option" << std::endl;
    throw;
  }
  _jetClustering->Init(_reconstructedParticleVec);
  _jetClustering->Clustering(_nJetMax);
  _y_NM1_N=_jetClustering->getYcutValueNM1_N();
  _y_N_NP1=_jetClustering->getYcutValueN_NP1();
  _d_NM1_N=_jetClustering->getDcutValueNM1_N();
  _d_N_NP1=_jetClustering->getDcutValueN_NP1();
}

void MyFastJetProcessor::processEvent( LCEvent * evt ) { 

  //_reconstructedParticleCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

  _nRun = evt->getRunNumber();
  _nEvt = evt->getEventNumber();

  LCCollection* inputParticleFlowCol=evt->getCollection(_inputCollection);
  _nInputRecoParticles =  inputParticleFlowCol->getNumberOfElements(); 
  _nHEJet=0;
  _reconstructedParticleVec.clear();
  for ( int ienflow=0; ienflow<_nInputRecoParticles ; ienflow++){
    ReconstructedParticle* enflow = dynamic_cast<ReconstructedParticle*>(inputParticleFlowCol->getElementAt( ienflow ));
    _reconstructedParticleVec.push_back(enflow);
  }
  doClustering();
  //_jetClustering=new JetClustering(GetAlgorithm(),fastjet::E_scheme,fastjet::Best);
  //if(NULL==_jetClustering) std::cout << "ok" << std::endl;
  WriteOutputJetInfo();
  _rootfile->cd();
  _tree->Fill();
  _jetEnergy.clear();
  _jetMass.clear();
  _jetPhi.clear();
  _jetEta.clear();
  _jetTheta.clear();
  _jetPx.clear();_jetPy.clear();_jetPz.clear();
  _WMass.clear();
  //streamlog_out( MESSAGE ) << "event number " << _nEvt << " processed" << std::endl;
  //getchar();
}

void  MyFastJetProcessor::check( LCEvent * evt ) { 

}


void MyFastJetProcessor::end(){ 
  _rootfile->cd("");
  _rootfile->Write();
  _rootfile->Close();
}
