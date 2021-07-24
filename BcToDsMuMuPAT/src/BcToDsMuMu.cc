

// -*- C++ -*-
//
// Package:    BcToDsMuMu
// Class:      BcToDsMuMu
// 
//=================================================
// Original by: Mohammad Mobassir Ameen           |
//<mohammad.moabassir.ameen@cern.ch>              |
//created:  Saturday Jul 3 2021                   |
//=================================================

// system include files
#include <memory>

// user include files
#include "bctodsmumu-analysis/BcToDsMuMuPAT/src/BcToDsMuMu.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/FWLite/interface/EventBase.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/CLHEP/interface/Migration.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"


#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"


#include <vector>
#include <utility>
#include <string>

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

//
BcToDsMuMu::BcToDsMuMu(const edm::ParameterSet& iConfig):

  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  //primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  genParticles_(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"))),
  packedGenToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter <edm::InputTag> ("packedGenParticles"))),
  //OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
  doMC_ ( iConfig.getUntrackedParameter<bool>("doMC",false) ),
  
  tree_(0), 

  nB(0),
  //B_mass(0), B_px(0), B_py(0), B_pz(0),
  
  //B_Phi_mass(0),
  B_Ds_mass(0),
  Ds_pt(0),
  pion_pt(0),
  //B_Ds_pt1(0), B_Ds_px1(0), B_Ds_py1(0), B_Ds_pz1(0), 
  //B_Ds_pt2(0), B_Ds_px2(0), B_Ds_py2(0), B_Ds_pz2(0), 
  //B_Ds_pt3(0), B_Ds_px3(0), B_Ds_py3(0), B_Ds_pz3(0),

  //B_Ds_px1_track(0), B_Ds_py1_track(0), B_Ds_pz1_track(0), 
  //B_Ds_px2_track(0), B_Ds_py2_track(0), B_Ds_pz2_track(0), 
  //B_Ds_px3_track(0), B_Ds_py3_track(0), B_Ds_pz3_track(0),

  //B_Ds_charge1(0), B_Ds_charge2(0),B_Ds_charge3(0)
  run(0), event(0),
  lumiblock(0)
  
  
{
   //now do what ever initialization is needed
}

BcToDsMuMu::~BcToDsMuMu()
{

}

// ------------ method called to for each event  ------------
void BcToDsMuMu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  // *********************************
  // Get event content information
  // *********************************  


  //Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 

  edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);

  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event();

  //we look for Ds candidates
  for(View<pat::PackedCandidate>::const_iterator iTrack1 = thePATTrackHandle->begin();
  iTrack1 != thePATTrackHandle->end(); ++iTrack1 )	
  {
	    
	    if(iTrack1->charge()==0) continue;		  
            //if(fabs(iTrack1->pdgId())!=321) continue;
	    //if(iTrack1->pt()<0.5) continue;
	    //if(iTrack1->pt()<0.8) continue;
	    if(iTrack1->pt()<2.0) continue; // this cut of tarck1 pt worked for official mc DsToPhi(MuMu)Pi
	    if(!(iTrack1->trackHighPurity())) continue;
	    
            for(View<pat::PackedCandidate>::const_iterator iTrack2 = iTrack1+1;
            iTrack2 != thePATTrackHandle->end(); ++iTrack2 ) 
              {

	      
	      if(iTrack1==iTrack2) continue;
	      if(iTrack2->charge()==0) continue;
              //if(fabs(iTrack1->pdgId())!=321) continue;
	      //if(iTrack2->pt()<0.5) continue;
	      if(iTrack2->pt()<2.0) continue; // this cut of track2 pt worked for official mc DsToPhi(MuMu)Pi
	      if(!(iTrack2->trackHighPurity())) continue;

	      if(iTrack1->charge() == iTrack2->charge()) continue;
	     
              // ********** Lets see first Phi Candidate.
              //TLorentzVector kaon14V,kaon24V,Phi4V;
              //Phi4V=kaon14V+kaon24V;
     
              
              
              //for(View<pat::PackedCandidate>::const_iterator iTrack3 = thePATTrackHandle->begin();
              //iTrack3 != thePATTrackHandle->end(); ++iTrack3 ) 
              
              for(View<pat::PackedCandidate>::const_iterator iTrack3 = iTrack2+1;
              iTrack3 != thePATTrackHandle->end(); ++iTrack3 )
                {
		  
		  if(iTrack3==iTrack2) continue;
		  if(iTrack3==iTrack1) continue; 
                  //if(fabs(iTrack1->pdgId())!=211) continue;
		  if(iTrack3->charge()==0) continue;
		  //if(iTrack3->pt()<0.2) continue;
		  if(iTrack3->pt()<1.0) continue; // this cut of track3 pt cut worked for official mc DsToPhi(MuMu)Pi
		  if(!(iTrack3->trackHighPurity())) continue;
		  
                  //Now let's see if these three tracks make a vertex		  
		  reco::TransientTrack pion1TT((*theB).build(iTrack1->pseudoTrack()));
		  reco::TransientTrack pion2TT((*theB).build(iTrack2->pseudoTrack()));
		  reco::TransientTrack pion3TT((*theB).build(iTrack3->pseudoTrack()));
		  

		
		ParticleMass pion_mass = 0.13957018;
		ParticleMass kaon_mass = 0.493677;
                //ParticleMass muon_mass = 0.10565837; //pdg mass
		float kaon_sigma = 1.6e-5;
		float pion_sigma = pion_mass*1.e-6;
		//float muon_sigma = muon_mass*1.e-6;
	// ------------ started Vertexing	
	        
                //Creating a KinematicParticleFactory
                KinematicParticleFactoryFromTransientTrack pFactory;

                //initial chi2 and ndf before kinematic fits.
                float chi = 0.;
                float ndf = 0.;
                vector<RefCountedKinematicParticle> pionParticles;
                try {
		  pionParticles.push_back(pFactory.particle(pion1TT,kaon_mass,chi,ndf,kaon_sigma)); // for Private mc BcToDs(KpKmPi)MuMu
		  //pionParticles.push_back(pFactory.particle(pion1TT,muon_mass,chi,ndf,muon_sigma)); //for official mc DsToPhi(MuMu)Pi
		  pionParticles.push_back(pFactory.particle(pion2TT,kaon_mass,chi,ndf,kaon_sigma)); //for Private mc BcToDs(KpKmPi)MuMu
		  //pionParticles.push_back(pFactory.particle(pion2TT,muon_mass,chi,ndf,muon_sigma)); //for official mc DsToPhi(MuMu)Pi
		  pionParticles.push_back(pFactory.particle(pion3TT,pion_mass,chi,ndf,pion_sigma)); // It is common for both
			
		}
                catch(...) {
		  std::cout<<" Exception caught ... continuing 3 "<<std::endl;
		  continue;
		}
                RefCountedKinematicTree DsVertexFitTree;
                KinematicParticleVertexFitter fitter; 
		try{
		  DsVertexFitTree = fitter.fit(pionParticles); 
		}
                catch(...) {
		  std::cout<<" Exception caught ... continuing 4 "<<std::endl;                   
		  continue;
		}
                if (!DsVertexFitTree->isValid()) 
		  {
		    //std::cout << "invalid vertex from the Ds vertex fit" << std::endl;
		    continue; 
		  }
	        DsVertexFitTree->movePointerToTheTop();	

                RefCountedKinematicParticle Ds_vFit_noMC = DsVertexFitTree->currentParticle();
                RefCountedKinematicVertex Ds_vFit_vertex_noMC = DsVertexFitTree->currentDecayVertex();
                
                if( Ds_vFit_vertex_noMC->chiSquared() < 0 )
		  { 
		    //std::cout << "negative chisq from ks fit" << endl;
		    continue;
		  }  

		//some loose cuts go here
		if(Ds_vFit_vertex_noMC->chiSquared()>50) continue;
		if(Ds_vFit_noMC->currentState().mass()<1.90 || Ds_vFit_noMC->currentState().mass()>2.04) continue;
		double Dspt= sqrt(pow(Ds_vFit_noMC->currentState().globalMomentum().x(),2) +pow(Ds_vFit_noMC->currentState().globalMomentum().y(),2));
                //if(Dspt<5) continue;
		double Ds_Prob_tmp  = TMath::Prob(Ds_vFit_vertex_noMC->chiSquared(),(int)Ds_vFit_vertex_noMC->degreesOfFreedom());
		if(Ds_Prob_tmp<0.01)
		  {
		    continue;
		  }


	// -------------------------------- End of Vertexing

                // ***************************
 	        // K+K-Pi+ invariant mass (before kinematic vertex fit)
 		// ***************************

 		//TLorentzVector pion4V, Ds4V; 
 		TLorentzVector kaon14V, kaon24V, pion4V, Ds4V; 
 		kaon14V.SetXYZM(iTrack1->px(),iTrack1->py(),iTrack1->pz(),kaon_mass); // for private mc
 		//kaon14V.SetXYZM(iTrack1->px(),iTrack1->py(),iTrack1->pz(),muon_mass); //Changed kaon_mass to muon_mass for official mc
 		kaon24V.SetXYZM(iTrack2->px(),iTrack2->py(),iTrack2->pz(),kaon_mass); // for private mc 
 		//kaon24V.SetXYZM(iTrack2->px(),iTrack2->py(),iTrack2->pz(),muon_mass); //Changed kaon_mass to muon_mass for official mc
 		pion4V.SetXYZM(iTrack3->px(),iTrack3->py(),iTrack3->pz(),pion_mass);



                //Phi4V=kaon14V+kaon24V;
	        //Ds4V=Phi4V+pion4V;
	        Ds4V=kaon14V+kaon24V+pion4V;
		//cout<< "Ds mass 1st : " << Ds4V.M() << endl;

                //if(Phi4V.M()<0.970 || Phi4V.M()>1.070) continue;
                if(Ds4V.M()<1.90 || Ds4V.M()>2.04) continue;
                //if(Ds4V.Pt()<0)continue;
	

                	
		// fill candidate variables now
		
		if(nB==0){
		   //nMu  = nMu_tmp;
		    //cout<< "Ds mass 2nd : " << Ds4V.M() << endl;
		} // end nB==0
		
		//B_Ds_mass->push_back(Ds4V.M());
		
		//B_Phi_mass->push_back(Phi4V.M());
		B_Ds_mass->push_back( Ds_vFit_noMC->currentState().mass() );
		//B_Phi_mass->push_back(Phi4V.M());
                Ds_pt->push_back(Ds4V.Pt());
	        pion_pt->push_back(iTrack3->pt());		
		//std::cout << " Ds mass Filled" << endl;
		
//		B_Ds_chi2->push_back(Ds_vFit_vertex_noMC->chiSquared());
		
		
		nB++;	
		
		
	    } //for iTrack3
	  }
}

   //fill the tree and clear the vectors
  //cout << "This is test1" << std::endl;
  if (nB > 0 ) 
    {
      std::cout << "filling tree " << nB<<std::endl;
      tree_->Fill();
    }
  //cout << "This is test2" << std::endl;
  // *********
  
  nB = 0;

  //B_Phi_mass->clear(); 
  B_Ds_mass->clear(); 
  Ds_pt->clear();
  pion_pt->clear(); 
  
     
  
}  
// ------------ method called once each job just before starting event loop  ------------

void 
BcToDsMuMu::beginJob()
{

  std::cout << "Beginning analyzer job with value of doMC = " << doMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","Ds+ -> K+k-Pi+ ntuple");

  tree_->Branch("nB",&nB,"nB/i");
  
  //tree_->Branch("B_Phi_mass", &B_Phi_mass);
  tree_->Branch("B_Ds_mass", &B_Ds_mass);
  tree_->Branch("Ds_pt", &Ds_pt);
  tree_->Branch("pion_pt", &pion_pt);

  tree_->Branch("run",        &run,       "run/I");
  tree_->Branch("event",        &event,     "event/I");
  tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");
}

// ------------ method called once each job just after ending the event loop  ------------
void BcToDsMuMu::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(BcToDsMuMu);
