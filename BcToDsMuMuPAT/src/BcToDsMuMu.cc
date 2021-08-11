

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

  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  //primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  genParticles_(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"))),
  packedGenToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter <edm::InputTag> ("packedGenParticles"))),

  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
  doMC_ ( iConfig.getUntrackedParameter<bool>("doMC",false) ),
  
  tree_(0), 

  nB(0), nMu(0),
  B_mass(0), B_px(0), B_py(0), B_pz(0),
  
  B_Ds_mass(0), B_Ds_px(0), B_Ds_py(0), B_Ds_pz(0),
  //B_Phi_mass(0),
  //B_mass(0),
  //B_Ds_mass(0),
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

  edm::Handle< View<pat::Muon> > thePATMuonHandle;
  iEvent.getByToken(dimuon_Label,thePATMuonHandle);

  edm::Handle<reco::GenParticleCollection> pruned;
  //edm::Handle<pat::PackedGenParticle> pruned; 
  //iEvent.getByToken(genCands_, pruned);
  iEvent.getByToken(genParticles_, pruned);

  edm::Handle<pat::PackedGenParticleCollection> packed;
  iEvent.getByToken(packedGenToken_,packed);

  //std::cout<<"size of the muon collection "<<thePATMuonHandle->size()<<" track size "<<theV0PtrHandle->size()<<" packed track "<<thePATTrackHandle->size()<<endl;

  //*********************************
  //  Get gen level information
  //*********************************  

  gen_b_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_ds_p4.SetPtEtaPhiM(0.,0.,0.,0.); 
  gen_kaon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_kaon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_pion_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_b_vtx.SetXYZ(0.,0.,0.);
  gen_jpsi_vtx.SetXYZ(0.,0.,0.);
  gen_b_ct = -9999.;

  if ( (isMC_ || OnlyGen_) && pruned.isValid() ) {
std::cout << "Output0: "<< std::endl;
    int foundit = 0;
    for (size_t i=0; i<pruned->size(); i++) {
       foundit = 0;
       const reco::Candidate *dau = &(*pruned)[i];

std::cout << "output1: "<< dau->pdgId() << std::endl;
       if ( (abs(dau->pdgId()) == 541) ) { //&& (dau->status() == 2) ) {          // Bc = 541
      
std::cout << "1: "<< dau->pdgId() << std::endl;

         foundit++;
         gen_b_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
	 gen_b_vtx.SetXYZ(dau->vx(),dau->vy(),dau->vz());
         for (size_t k=0; k<dau->numberOfDaughters(); k++) {
            const reco::Candidate *gdau = dau->daughter(k);
            if (gdau->pdgId()==443 ) { //&& gdau->status()==2) {                 // J/psi = 443  === J/psi -> mu+ mu-
            foundit++;
            gen_jpsi_vtx.SetXYZ(gdau->vx(),gdau->vy(),gdau->vz());
	    gen_b_ct = GetLifetime(gen_b_p4,gen_b_vtx,gen_jpsi_vtx);
            int nm=0;
            for (size_t l=0; l<gdau->numberOfDaughters(); l++) {
               const reco::Candidate *mm = gdau->daughter(l);
               if (mm->pdgId()==13) { foundit++;                               // mu- = 13
                 if (mm->status()!=1) {
                   for (size_t m=0; m<mm->numberOfDaughters(); m++) {
                      const reco::Candidate *mu = mm->daughter(m);
                      if (mu->pdgId()==13 ) { //&& mu->status()==1) {
                        nm++;
		        gen_muon1_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
		        break;
                      }
                    }
                  } else {
                    gen_muon1_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
		    nm++;
                  }
                }
                if (mm->pdgId()==-13) { foundit++;
		if (mm->status()!=1) {
		  for (size_t m=0; m<mm->numberOfDaughters(); m++) {
		    const reco::Candidate *mu = mm->daughter(m);
		    if (mu->pdgId()==-13 ) { //&& mu->status()==1) {
		      nm++;
		      gen_muon2_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
		      break;
		    }
		  }
		} else {
		  gen_muon2_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
		  nm++;
		}
	      }
	    }
	    if (nm==2) gen_jpsi_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
	    else foundit-=nm;
	  }

          for (size_t lk=0; lk<packed->size(); lk++) {
             const reco::Candidate * dauInPrunedColl = (*packed)[lk].mother(0);
             int stable_id = (*packed)[lk].pdgId();
             if (dauInPrunedColl != nullptr && isAncestor(gdau,dauInPrunedColl)) {
               if(stable_id == 321) {foundit++;                             // K+ = 321
                 gen_kaon1_p4.SetPtEtaPhiM((*packed)[lk].pt(),(*packed)[lk].eta(),(*packed)[lk].phi(),(*packed)[lk].mass());
               }
               if(stable_id == -321){ foundit++;                            // K- = -321
		 gen_kaon2_p4.SetPtEtaPhiM((*packed)[lk].pt(),(*packed)[lk].eta(),(*packed)[lk].phi(),(*packed)[lk].mass());
	       }
               if(stable_id == 211){ foundit++;                             // Pi+ = 211
		 gen_pion_p4.SetPtEtaPhiM((*packed)[lk].pt(),(*packed)[lk].eta(),(*packed)[lk].phi(),(*packed)[lk].mass());
	       }
	     }
	   }
	 
	} // for (size_t k
      } 
      if (foundit>=7) break;
    } // for i
    if (foundit!=7) {
      gen_b_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_ds_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_b_vtx.SetXYZ(0.,0.,0.);
      gen_jpsi_vtx.SetXYZ(0.,0.,0.);
      gen_b_ct = -9999.;
      std::cout << "Does not found the given decay " << run << "," << event << " foundit=" << foundit << std::endl; // sanity check
    }
  }

         


  

// ====================       yaha tak gen level onformation ka kaam karte hai/

  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event();

  //Let's begin by looking for mu+ and mu-
  unsigned int nMu_tmp = thePATMuonHandle->size();
  
  for(edm::View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1)
    {
      for(edm::View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2)
        {
      
          if(iMuon1==iMuon2) continue;
         
          //opposite charge 
          if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue;
          
          TrackRef glbTrackP;	  
	  TrackRef glbTrackM;	

          if(iMuon1->charge() == 1){glbTrackP = iMuon1->track();}
	  if(iMuon1->charge() == -1){glbTrackM = iMuon1->track();}

          if(iMuon2->charge() == 1) {glbTrackP = iMuon2->track();}
	  if(iMuon2->charge() == -1){glbTrackM = iMuon2->track();}

          if( glbTrackP.isNull() || glbTrackM.isNull() ) 
	    {
	      //std::cout << "continue due to no track ref" << endl;
	      continue;
	    } 

          if(iMuon1->track()->pt()<4.0) continue;
	  if(iMuon2->track()->pt()<4.0) continue;

          if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
	  if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;

          //Let's check the vertex and mass
          reco::TransientTrack muon1TT((*theB).build(glbTrackP));
	  reco::TransientTrack muon2TT((*theB).build(glbTrackM));

          // *****  Trajectory states to calculate DCA for the 2 muons *********************
          FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
	  FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();

          if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() ) continue;
          
          // Measure distance between tracks at their closest approach
          ClosestApproachInRPhi cApp;
	  cApp.calculate(mu1State, mu2State);
	  if( !cApp.status() ) continue;
	  float dca = fabs( cApp.distance() );	  
	  //if (dca < 0. || dca > 0.5) continue;
	  //cout<<" closest approach  "<<dca<<endl;
  
          // *****  end DCA for the 2 muons *********************


          //The mass of a muon and the insignificant mass sigma 
          //to avoid singularities in the covariance matrix.
          ParticleMass muon_mass = 0.10565837; //pdg mass
          //ParticleMass psi_mass = 3.096916;
          float muon_sigma = muon_mass*1.e-6;
          //float psi_sigma = psi_mass*1.e-6;

          //Creating a KinematicParticleFactory
          KinematicParticleFactoryFromTransientTrack pFactory;

          //initial chi2 and ndf before kinematic fits.
          float chi = 0.;
	  float ndf = 0.;
	  vector<RefCountedKinematicParticle> muonParticles;
	  try {
	    muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	    muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	  }
	  catch(...) { 
	    std::cout<<" Exception caught ... continuing 1 "<<std::endl; 
	    continue;
	  }

          KinematicParticleVertexFitter fitter;

          RefCountedKinematicTree psiVertexFitTree;
	  try {
	    psiVertexFitTree = fitter.fit(muonParticles); 
	  }
	  catch (...) { 
	    std::cout<<" Exception caught ... continuing 2 "<<std::endl; 
	    continue;
	  }

          if (!psiVertexFitTree->isValid()) 
	    {
	      //std::cout << "caught an exception in the psi vertex fit" << std::endl;
	      continue; 
	    }

          psiVertexFitTree->movePointerToTheTop();

          RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();
	  RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();

          if( psi_vFit_vertex_noMC->chiSquared() < 0 )
	    {
	      std::cout << "negative chisq from psi fit" << endl;
	      continue;
	    }

           double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
	  if(J_Prob_tmp<0.01)
	    {
	      continue;
	    }

  //we look for Ds candidates
  for(View<pat::PackedCandidate>::const_iterator iTrack1 = thePATTrackHandle->begin();
  iTrack1 != thePATTrackHandle->end(); ++iTrack1 )	
  {
	    
	    if(iTrack1->charge()==0) continue;		  
            //if(fabs(iTrack1->pdgId())!=321) continue;
	    //if(iTrack1->pt()<1.0) continue;  // changed for 2 to 1 for new MC
	    //if(iTrack1->pt()<0.8) continue;
	    //if(iTrack1->pt()<2.0) continue; // this cut of tarck1 pt worked for official mc DsToPhi(MuMu)Pi
	    if(!(iTrack1->trackHighPurity())) continue;
	    
            for(View<pat::PackedCandidate>::const_iterator iTrack2 = iTrack1+1;
            iTrack2 != thePATTrackHandle->end(); ++iTrack2 ) 
              {

	      
	      if(iTrack1==iTrack2) continue;
	      if(iTrack2->charge()==0) continue;
              //if(fabs(iTrack1->pdgId())!=321) continue;
	      //if(iTrack2->pt()<1.0) continue; // changed for 2 to 1 for new MC 
	      //if(iTrack2->pt()<2.0) continue; // this cut of track2 pt worked for official mc DsToPhi(MuMu)Pi
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
		  //if(iTrack3->pt()<1.0) continue; // this cut of track3 pt cut worked for official mc DsToPhi(MuMu)Pi
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
		//float kaon_sigma = kaon_mass*1.e-6;
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
                //VertexFitTree->movePointerToTheFirstChild();
                DsVertexFitTree->movePointerToTheFirstChild();
                RefCountedKinematicParticle T1CandMC = DsVertexFitTree->currentParticle();

                DsVertexFitTree->movePointerToTheNextChild();
                RefCountedKinematicParticle T2CandMC = DsVertexFitTree->currentParticle();

                DsVertexFitTree->movePointerToTheNextChild();
                RefCountedKinematicParticle T3CandMC = DsVertexFitTree->currentParticle();

                DsVertexFitTree->movePointerToTheTop();
                RefCountedKinematicParticle ks0_vFit_withMC = DsVertexFitTree->currentParticle();

                //Now we are ready to combine!

                //cout<<"Dspt "<<Dspt<<endl;

                vector<RefCountedKinematicParticle> vFitMCParticles;
                vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
                vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
                //vFitMCParticles.push_back(pFactory.particle(pion1TT,kaon_mass,chi,ndf,kaon_sigma));
                //vFitMCParticles.push_back(pFactory.particle(pion2TT,kaon_mass,chi,ndf,kaon_sigma));                

                vFitMCParticles.push_back(ks0_vFit_withMC);

                KinematicConstrainedVertexFitter kcvFitter;
                RefCountedKinematicTree vertexFitTree = kcvFitter.fit(vFitMCParticles);
                if (!vertexFitTree->isValid()) {
                  //std::cout << "B MC fit vertex is not valid" << endl;
                  continue;
                }

                vertexFitTree->movePointerToTheTop();

                RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
                RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
                if (!bDecayVertexMC->vertexIsValid()){
                  //std::cout << "B MC fit vertex is not valid" << endl;
                  continue;
                }

                if(bCandMC->currentState().mass()<5.7 || bCandMC->currentState().mass()>6.8) continue;

                if(bDecayVertexMC->chiSquared()<0 || bDecayVertexMC->chiSquared()>50 )
                {
                    //std::cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
                    continue;
                }

                double B_Prob_tmp  = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
                if(B_Prob_tmp<0.01)
                  {
                    continue;
                  }

                std::cout << "B mass "<<bCandMC->currentState().mass()<< std::endl;

                vertexFitTree->movePointerToTheFirstChild();
                RefCountedKinematicParticle mu1CandMC = vertexFitTree->currentParticle();
                vertexFitTree->movePointerToTheNextChild();
                RefCountedKinematicParticle mu2CandMC = vertexFitTree->currentParticle();

                vertexFitTree->movePointerToTheNextChild();
                RefCountedKinematicParticle DsCandMC = vertexFitTree->currentParticle();

                KinematicParameters psiMu1KP = mu1CandMC->currentState().kinematicParameters();
                KinematicParameters psiMu2KP = mu2CandMC->currentState().kinematicParameters();
                KinematicParameters psiMupKP;
                KinematicParameters psiMumKP;

                if ( mu1CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu1KP;
                if ( mu1CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu1KP;
                if ( mu2CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu2KP;
                if ( mu2CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu2KP;

                GlobalVector Jp1vec(mu1CandMC->currentState().globalMomentum().x(),
                                    mu1CandMC->currentState().globalMomentum().y(),
                                    mu1CandMC->currentState().globalMomentum().z());

                GlobalVector Jp2vec(mu2CandMC->currentState().globalMomentum().x(),
                                    mu2CandMC->currentState().globalMomentum().y(),
                                    mu2CandMC->currentState().globalMomentum().z());

                GlobalVector Dsp1vec(T1CandMC->currentState().globalMomentum().x(),
                                     T1CandMC->currentState().globalMomentum().y(),
                                     T1CandMC->currentState().globalMomentum().z());

                GlobalVector Dsp2vec(T2CandMC->currentState().globalMomentum().x(),
                                     T2CandMC->currentState().globalMomentum().y(),
                                     T2CandMC->currentState().globalMomentum().z());

                GlobalVector Dsp3vec(T3CandMC->currentState().globalMomentum().x(),
                                     T3CandMC->currentState().globalMomentum().y(),
                                     T3CandMC->currentState().globalMomentum().z());

                
                //KinematicParameters DsPi1KP = T1CandMC->currentState().kinematicParameters();
                //KinematicParameters DsPi2KP = T2CandMC->currentState().kinematicParameters();
                //KinematicParameters DsPi3KP = T3CandMC->currentState().kinematicParameters();




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
                //if(Ds4V.M()<1.90 || Ds4V.M()>2.04) continue;
                if(Ds4V.Pt()<0)continue; // uncomment for the private mc
	

                	
		// fill candidate variables now
		
		if(nB==0){
		   nMu  = nMu_tmp;
		    //cout<< "Ds mass 2nd : " << Ds4V.M() << endl;
		} // end nB==0

                B_mass->push_back(bCandMC->currentState().mass());           
                B_px->push_back(bCandMC->currentState().globalMomentum().x());
		B_py->push_back(bCandMC->currentState().globalMomentum().y());
		B_pz->push_back(bCandMC->currentState().globalMomentum().z());
		
		//B_Ds_mass->push_back(Ds4V.M());
		
		//B_Phi_mass->push_back(Phi4V.M());
		B_Ds_mass->push_back( Ds_vFit_noMC->currentState().mass() );
                B_Ds_px->push_back( Ds_vFit_noMC->currentState().globalMomentum().x() );
                B_Ds_py->push_back( Ds_vFit_noMC->currentState().globalMomentum().x() );
                B_Ds_pz->push_back( Ds_vFit_noMC->currentState().globalMomentum().x() );
		//B_Phi_mass->push_back(Phi4V.M());
                Ds_pt->push_back(Ds4V.Pt());
	        pion_pt->push_back(iTrack3->pt());		
		//std::cout << " Ds mass Filled" << endl;
		
//		B_Ds_chi2->push_back(Ds_vFit_vertex_noMC->chiSquared());
		
		
		nB++;	
       
                pionParticles.clear();
                muonParticles.clear();
                vFitMCParticles.clear();
		
	    } //for iTrack3
	  }
      }
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
  
  nB = 0; nMu = 0;

  B_mass->clear(); B_px->clear();    B_py->clear();    B_pz->clear();
  B_Ds_mass->clear(); B_Ds_px->clear(); B_Ds_py->clear(); B_Ds_pz->clear();

  //B_Phi_mass->clear(); 
  B_Ds_mass->clear(); 
  Ds_pt->clear();
  pion_pt->clear(); 
  
     
  
}  

bool BcToDsMuMu::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
    if (ancestor == particle ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(ancestor,particle->mother(i))) return true;
    }
    return false;
}

double BcToDsMuMu::GetLifetime(TLorentzVector b_p4, TVector3 production_vtx, TVector3 decay_vtx) {
   TVector3 pv_dv = decay_vtx - production_vtx;
   TVector3 b_p3  = b_p4.Vect();
   pv_dv.SetZ(0.);
   b_p3.SetZ(0.);
   Double_t lxy   = pv_dv.Dot(b_p3)/b_p3.Mag();
   return lxy*b_p4.M()/b_p3.Mag();
}
// ------------ method called once each job just before starting event loop  ------------

void 
BcToDsMuMu::beginJob()
{
  
  std::cout << "Beginning analyzer job with value of isMC_ = " << isMC_ << std::endl;  

  //std::cout << "Beginning analyzer job with value of doMC = " << doMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","Ds+ -> K+k-Pi+ ntuple");

  tree_->Branch("nB",&nB,"nB/i");
 
  tree_->Branch("B_mass", &B_mass); 
  tree_->Branch("B_px", &B_px);
  tree_->Branch("B_py", &B_py);
  tree_->Branch("B_pz", &B_pz);

  //tree_->Branch("B_Phi_mass", &B_Phi_mass);
  tree_->Branch("B_Ds_mass", &B_Ds_mass);
  tree_->Branch("B_Ds_px", &B_Ds_px);
  tree_->Branch("B_Ds_py", &B_Ds_py);
  tree_->Branch("B_Ds_pz", &B_Ds_pz);

  tree_->Branch("Ds_pt", &Ds_pt);
  tree_->Branch("pion_pt", &pion_pt);

  tree_->Branch("run",        &run,       "run/I");
  tree_->Branch("event",        &event,     "event/I");
  tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");

// gen
if (isMC_) {
     tree_->Branch("gen_b_p4",     "TLorentzVector",  &gen_b_p4);
     tree_->Branch("gen_ds_p4",   "TLorentzVector",  &gen_ds_p4);
     tree_->Branch("gen_kaon1_p4",  "TLorentzVector",  &gen_kaon1_p4);
     tree_->Branch("gen_kaon2_p4",  "TLorentzVector",  &gen_kaon2_p4);
     tree_->Branch("gen_pion_p4",  "TLorentzVector",  &gen_pion_p4);
     tree_->Branch("gen_jpsi_p4",   "TLorentzVector",  &gen_jpsi_p4);
     tree_->Branch("gen_muon1_p4",  "TLorentzVector",  &gen_muon1_p4);
     tree_->Branch("gen_muon2_p4",  "TLorentzVector",  &gen_muon2_p4);
     tree_->Branch("gen_b_vtx",    "TVector3",        &gen_b_vtx);
     tree_->Branch("gen_jpsi_vtx",  "TVector3",        &gen_jpsi_vtx);
     tree_->Branch("gen_b_ct",     &gen_b_ct,        "gen_b_ct/F");
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void BcToDsMuMu::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(BcToDsMuMu);
