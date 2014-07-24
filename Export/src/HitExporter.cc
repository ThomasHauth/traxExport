// system include files
#include <memory>
#include <vector>
#include <map>
#include <deque>
#include <queue>
#include <iostream>
#include <set>
#include <limits>
#include <array>
#include <algorithm>

#include <math.h>

// cmssw
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"

#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "DataFormats/TrackingRecHit/interface/InvalidTrackingRecHit.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "traxExport/Export/interface/DetIdClassification.h"
#include "traxExport/Export/interface/tHit.h"
#include "traxExport/Export/interface/PEventStore.h"

//root
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"

//
// class declaration
//

using namespace edm;

// can be configured to generate triplets in various parts of
// the detector

class HitExporter: public edm::EDAnalyzer {
public:

	explicit HitExporter(const edm::ParameterSet&);
	virtual ~HitExporter();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

	virtual void beginJob();
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob();

	virtual void beginRun(edm::Run const&, edm::EventSetup const&);
	virtual void endRun(edm::Run const&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock const&,
			edm::EventSetup const&);
	virtual void endLuminosityBlock(edm::LuminosityBlock const&,
			edm::EventSetup const&);

	// ----------member data ---------------------------
	edm::InputTag m_inputTagRecoTracks; //used to select what tracks to read from configuration file
	edm::InputTag m_inputMatchedRecHits;
	double m_minPt;
	std::string hitFile;
	EventStoreOutput hitOutput;
//	edm::InputTag simTracks_;
//	std::vector<std::string> simHits_; //simHits to load
//	double m_minPt; //min Pt cut
//	double maxPt_; //max Pt cut
//	double simHitCut_; //ignore SimHits with low momentum
//	bool usePrediciton_; //use predicted search window
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//edm::Handle < edm::SimTrackContainer >
HitExporter::HitExporter(const edm::ParameterSet& iConfig) :

// TODO: use Simtrack so get MC truth about tracks
		m_inputTagRecoTracks(
				iConfig.getUntrackedParameter < edm::InputTag > ("recoTracks")),
		m_inputMatchedRecHits(
				iConfig.getUntrackedParameter < edm::InputTag > ("matchedRecHits")),
		m_minPt(
				iConfig.getUntrackedParameter<double>("minPt")),
		hitFile(
				iConfig.getUntrackedParameter<std::string>("hitFile")),
		hitOutput(hitFile, true){

	edm::Service < TFileService > fs;
	/*new TH1D ( name.c_str(), histInfo.title.c_str() ,
	 (histInfo.Xhigh - histInfo.Xlow) / histInfo.Xres + 1 , histInfo.Xlow, histInfo.Xhigh,
	 (histInfo.Yhigh - histInfo.Ylow) / histInfo.Yres + 1, histInfo.Ylow, histInfo.Yhigh);*/

	//TTree * myTree = fs->make < TTree > ("Triplets", "");

	//m_storeTriplet.reset(new StoreTriplets(myTree));
	/*int someInt = 23;
	 myTree->Branch("SomeInt", &someInt, "SomeInt/I");
	 myTree->Fill();*/
//myTree->Wi
	//createTrackletHistograms(m_histIdealRecoTracklets, recoIdealDir);
}

HitExporter::~HitExporter() {

// do anything here that needs to be done at desctruction time
// (e.g. close files, deallocate resources etc.)
	//myTree->

	hitOutput.Close();

}

enum class HitType {
		Unknown, Pixel, Strip1d, Strip2d, StripMatched, Invalid
	};

HitType getHitType(TrackingRecHit const* firstHit) {
	HitType hitType = HitType::Unknown;

	SiPixelRecHit const* pixelHit =
			dynamic_cast<SiPixelRecHit const*>(firstHit);
	SiStripRecHit1D const* strip1dHit =
			dynamic_cast<SiStripRecHit1D const*>(firstHit);
	SiStripRecHit2D const* strip2dHit =
			dynamic_cast<SiStripRecHit2D const*>(firstHit);
	SiStripMatchedRecHit2D const* stripMatched2dHit =
			dynamic_cast<SiStripMatchedRecHit2D const*>(firstHit);
	InvalidTrackingRecHit const* invalidHit =
			dynamic_cast<InvalidTrackingRecHit const*>(firstHit);
	if (pixelHit != nullptr) {
		hitType = HitType::Pixel;
	} else if (strip1dHit != nullptr) {
		hitType = HitType::Strip1d;
	} else if (invalidHit != nullptr) {
		hitType =HitType::Invalid;
	} else if (strip2dHit != nullptr) {
		hitType = HitType::Strip2d;
	} else if (stripMatched2dHit != nullptr) {
		hitType = HitType::StripMatched;
	} else {
		std::cout << typeid(*firstHit).name() << std::endl;
		std::cout << "Fatal: Unknow strip type" << std::endl;
		assert(false);
	}

	return hitType;
}

void HitExporter::analyze(const edm::Event& evt,
		const edm::EventSetup& iSetup) {

	edm::ESHandle < TrackerGeometry > tracker;
	iSetup.get<TrackerDigiGeometryRecord>().get(tracker);

	edm::Handle < edm::View<reco::Track> > recoTracks;
	evt.getByLabel(m_inputTagRecoTracks, recoTracks);

	//edm::Handle < SiStripMatchedRecHit2DCollection > matchedRecHits;
	//evt.getByLabel(m_inputMatchedRecHits, matchedRecHits);

	edm::ESHandle<TransientTrackingRecHitBuilder> recHitBuilder;
	iSetup.get<TransientRecHitRecord>().get("WithAngleAndTemplate",recHitBuilder);

	TrackerHitAssociator associater(evt);

	//get all RecHits
	std::vector<const TrackingRecHit*> recoHits;

	//pixel hits
	edm::Handle<SiPixelRecHitCollection> pixelRecHits;
	evt.getByLabel("siPixelRecHits", pixelRecHits);
	for(SiPixelRecHitCollection::const_iterator iDet = pixelRecHits->begin();
			iDet != pixelRecHits->end(); ++iDet){

		SiPixelRecHitCollection::DetSet det = *iDet;
		for(SiPixelRecHitCollection::DetSet::const_iterator iHit = det.begin();
				iHit != det.end(); ++iHit){

			recoHits.push_back(iHit);
		}

	}

	//matched silicon hits
	edm::Handle<SiStripMatchedRecHit2DCollection> matchedRecHits;
	evt.getByLabel("siStripMatchedRecHits","matchedRecHit", matchedRecHits);
	for(SiStripMatchedRecHit2DCollection::const_iterator iDet = matchedRecHits->begin();
			iDet != matchedRecHits->end(); ++ iDet){

		SiStripMatchedRecHit2DCollection::DetSet det = *iDet;
		for(SiStripMatchedRecHit2DCollection::DetSet::const_iterator iHit = det.begin();
				iHit != det.end(); ++iHit){

			recoHits.push_back(iHit);
		}

	}

	edm::Handle<std::vector<SimTrack> > simTracks;
	evt.getByLabel("g4SimHits",simTracks);

	PB_Event::PEvent pEvent;

	pEvent.set_lumisection(evt.id().luminosityBlock());
	pEvent.set_runnumber(evt.id().run());
	pEvent.set_eventnumber(evt.id().event());

	uint nHit = 0;
	for(auto hit : recoHits){

		DetId detID = hit->geographicalId()();

		if(!detID)
			continue;

		DetIDClassification hitClass(detID);

		if(!hitClass.isBarrel())
			continue;

		GlobalPoint globalPos = recHitBuilder->build(hit)->globalPosition();
		tHit mHit(globalPos.x(), globalPos.y(), globalPos.z(), hit->geographicalId().rawId(), *hit, &associater);

		PB_Event ::PHit * pHit = pEvent.add_hits();
		pHit->set_detectorid(mHit.detID);
		pHit->set_layer(hitClass.getLayer());
		pHit->set_hitid(nHit++);
		pHit->set_detectortype(PB_Event::DetectorType::BARREL);

		pHit->mutable_position()->set_x(mHit.globalX);
		pHit->mutable_position()->set_y(mHit.globalY);
		pHit->mutable_position()->set_z(mHit.globalZ);

		pHit->set_simtrackid(mHit.simTrackId);

		auto simTrack = std::find_if(simTracks->begin(), simTracks->end(),
				[mHit] ( SimTrack const& a ) {
			return a.trackId() == mHit.simTrackId;
		});

		if(simTrack != simTracks->end()){

			GlobalVector trackPt = GlobalVector(simTrack->momentum().x(),
					simTrack->momentum().y(),
					simTrack->momentum().z());
			pHit->set_simtrackpt(trackPt.perp());
		}


	}

	std::cout << "Storing event " << pEvent.eventnumber() << " with " << pEvent.hits_size() << " hits" << std::endl;

	hitOutput.storeElement(pEvent);

}

// ------------ method called once each job just before starting event loop  ------------
void HitExporter::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void HitExporter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void HitExporter::beginRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a run  ------------
void HitExporter::endRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when starting to processes a luminosity block  ------------
void HitExporter::beginLuminosityBlock(edm::LuminosityBlock const&,
		edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------
void HitExporter::endLuminosityBlock(edm::LuminosityBlock const&,
		edm::EventSetup const&) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HitExporter::fillDescriptions(
		edm::ConfigurationDescriptions& descriptions) {
//The following says we do not know what parameters are allowed so do no validation
// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);

//Specify that only 'tracks' is allowed
//To use, remove the default given above and uncomment below
//edm::ParameterSetDescription desc;
//desc.addUntracked("tracks","generalTracks");
//desc.addUntracked("useSimulated",false);
//desc.addOptionalUntracked("simulatedTracks", "g4SimHits");
//descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE (HitExporter);
