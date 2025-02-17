// -*- C++ -*-
//
// Package:    L1Trigger/L1CaloTrigger
// Class:      Phase2L1TJetSeedReductionProducer
//
/**\class Phase2L1TJetSeedReductionProducer Phase2L1TJetSeedReductionProducer.cc L1Trigger/L1CaloTrigger/plugins/Phase2L1TJetSeedReductionProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Liam Robertshaw [SP/JB 2209]
//         Created:  Mon, 17 Feb 2025 15:00:41 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

//
// class declaration
//

class Phase2L1TJetSeedReductionProducer : public edm::stream::EDProducer<> {
public:
  explicit Phase2L1TJetSeedReductionProducer(const edm::ParameterSet&);
  ~Phase2L1TJetSeedReductionProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;


  // void endStream() override;

  //void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //void endRun(edm::Run const&, edm::EventSetup const&) override;
  //void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Phase2L1TJetSeedReductionProducer::Phase2L1TJetSeedReductionProducer(const edm::ParameterSet& cfg) 
  : seedsToken(consumes<l1t::PFCandidateCollection>(cfg.getParameter<edm::InputTag>("histoSeeds"))),
    outputCollectionName_(iConfig.getParameter<std::string>("outputCollectionName")) 
{
  //register your products
  /* Examples
  produces<ExampleData2>();

  //if do put with a label
  produces<ExampleData2>("label");
 
  //if you want to put into the Run
  produces<ExampleData2,InRun>();
  */
  //now do what ever other initialization is needed
  produces<l1t::PFCandidateCollection>(outputCollectionName_);
}


Phase2L1TJetSeedReductionProducer::~Phase2L1TJetSeedReductionProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void Phase2L1TJetSeedReductionProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // GET SEEDS FROM EDM INPUT
  edm::Handle<l1t::PFCandidateCollection> seedsHandle;
  std::vector<edm::Ptr<l1t::PFCandidate>> seeds;
  iEvent.getByToken(seedsToken, seedsHandle);
  for (unsigned i = 0; i < (*seedsHandle).size(); i++) {
    seeds.push_back(edm::Ptr<l1t::PFCandidate>(seedsHandle, i));
  }

  //code
  
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void Phase2L1TJetSeedReductionProducer::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void Phase2L1TJetSeedReductionProducer::endStream() {
  // please remove this method if not needed
}

/*
// ------------ method called when starting to processes a run  ------------
/*
void
Phase2L1TJetSeedReductionProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
Phase2L1TJetSeedReductionProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
Phase2L1TJetSeedReductionProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
Phase2L1TJetSeedReductionProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Phase2L1TJetSeedReductionProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("seeds", edm::InputTag('l1tHSC8PFL1PuppiSimSeedProducer', 'WIDEHSCEMUSEEDS'));
  
  descriptions.add(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Phase2L1TJetSeedReductionProducer);
