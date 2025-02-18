// -*- C++ -*-
//
// Package:    L1Trigger/L1CaloTrigger
// Class:      Phase2L1TJetSeedReductionProducer
//
/**\class Phase2L1TJetSeedReductionProducer Phase2L1TJetSeedReductionProducer.cc L1Trigger/L1CaloTrigger/plugins/Phase2L1TJetSeedReductionProducer.cc
*/

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"


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
  edm::EDGetTokenT<std::vector<l1t::PFCandidate>> seedsToken;
  std::string outputCollectionName_;

  l1t::PFCandidateCollection copyInputToOutput(l1t::PFCandidateCollection seeds);

};

Phase2L1TJetSeedReductionProducer::Phase2L1TJetSeedReductionProducer(const edm::ParameterSet& cfg) 
  : seedsToken(consumes<l1t::PFCandidateCollection>(cfg.getParameter<edm::InputTag>("seeds"))),
    outputCollectionName_(cfg.getParameter<std::string>("outputCollectionName")) 
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

void Phase2L1TJetSeedReductionProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // GET SEEDS FROM EDM INPUT
  edm::Handle<l1t::PFCandidateCollection> seedsHandle;
  // std::vector<edm::Ptr<l1t::PFCandidate>> seeds;
  l1t::PFCandidateCollection seeds;
  iEvent.getByToken(seedsToken, seedsHandle);
  for (unsigned i = 0; i < (*seedsHandle).size(); i++) {
    edm::Ptr<l1t::PFCandidate> seedPtr(seedsHandle, i);    // get pointer to seed i
    seeds.push_back(*seedPtr);    // dereference seed i pointer to get seed i and add to seeds vector
  }
  // WE NOW HAVE A VECTOR OF THE SEEDS PASSED TO THE MODULE

  // NOW WE WANT TO WRITE THE LOGIC WHICH WILL REDUCE THE SEEDS TO PRODUCE THE WIDE JET SEEDS
  // AVOID WRITING ALL THE CODE IN THE PRODUCE FUNCTION, INSTEAD BREAK INTO SMALLER MORE MODULAR FUNCTIONS (I.E. DELTA R MATCHING, ETC)
  // AS AN EXAMPLE, I HAVE WRITTEN A FUNCTION WHICH COPIES THE VECTOR OF SEEDS TO A NEW VECTOR OF SEEDS:
  l1t::PFCandidateCollection reducedSeeds = copyInputToOutput(seeds);

  int i = 0;
  for(auto& seed : reducedSeeds) {
    i += 1;
    std::cout << "Seed "<< i << ": " << seed.pt() << std::endl;
    // std::cout << "Seed "<< i << ": " << seed.pt() << std::endl;
  }

  std::unique_ptr<l1t::PFCandidateCollection> reducedSeedsPtr;//(new l1t::PFCandidateCollection);
  reducedSeedsPtr->swap(reducedSeeds);

  iEvent.put(std::move(reducedSeedsPtr), outputCollectionName_);    // NOTE CAN ONLY PUT POINTERS INTO THE EVENT, NOT THE OBJECTS THEMSELVES
}

/* DESTRUCTOR */
Phase2L1TJetSeedReductionProducer::~Phase2L1TJetSeedReductionProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

l1t::PFCandidateCollection Phase2L1TJetSeedReductionProducer::copyInputToOutput(l1t::PFCandidateCollection seeds){
  l1t::PFCandidateCollection reducedSeeds;
  std::copy(seeds.begin(), seeds.end(), std::back_inserter(reducedSeeds));
  return reducedSeeds;
}


/*
// ------------ method called once each stream before processing any runs, lumis or events  ------------
void Phase2L1TJetSeedReductionProducer::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void Phase2L1TJetSeedReductionProducer::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
void Phase2L1TJetSeedReductionProducer::beginRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a run  ------------
void Phase2L1TJetSeedReductionProducer::endRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when starting to processes a luminosity block  ------------
void Phase2L1TJetSeedReductionProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------

void Phase2L1TJetSeedReductionProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Phase2L1TJetSeedReductionProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("seeds", edm::InputTag("l1tPhase1JetSeedProducer", "seeds"));
  desc.add<std::string>("outputCollectionName", "reducedSeeds");
  
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Phase2L1TJetSeedReductionProducer);
