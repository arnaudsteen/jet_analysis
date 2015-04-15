
#ifndef MyFastJetProcessor_h
#define MyFastJetProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <cstring>
#include <vector>
#include <EVENT/LCCollection.h>
#include <EVENT/LCParameters.h>
#include <EVENT/ReconstructedParticle.h>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

#include "JetClustering.h"

using namespace lcio ;
using namespace marlin ;

class MyFastJetProcessor : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new MyFastJetProcessor ; }
  
  
  MyFastJetProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  void doClustering();
  void WriteOutputJetInfo();
  inline void SetAlgorithm(fastjet::JetAlgorithm algo){_jetAlgo=algo;}
  inline fastjet::JetAlgorithm &GetAlgorithm(){return _jetAlgo;}
  static bool SortJetWithEnergy(fastjet::PseudoJet p1, fastjet::PseudoJet p2){return p1.e()>p2.e();}
 protected:

  int _nRun ;
  int _nEvt ;
  /** Input collection name.
   */
  std::string _inputCollection;
  std::string sAlgorithm;
  double _RPar;
  double _xtra_param;
  double _eJet;
  int _nJetMax;

 private:
  //  int numElements;
  LCCollection * col;
  TFile *_rootfile;
  TTree *_tree;
  std::vector<EVENT::ReconstructedParticle*> _reconstructedParticleVec;
  JetClustering* _jetClustering;
  fastjet::JetAlgorithm _jetAlgo;

  int _nInputRecoParticles;
  int _nHEJet;
  std::vector<double> _jetEnergy;
  std::vector<double> _jetMass;
  std::vector<double> _jetPhi;
  std::vector<double> _jetEta;
  std::vector<double> _jetTheta;
  std::vector<double> _jetPx;
  std::vector<double> _jetPy;
  std::vector<double> _jetPz;
  double _y_NM1_N;
  double _y_N_NP1;
  double _d_NM1_N;
  double _d_N_NP1;
  std::vector<double> _WMass;

} ;

#endif
