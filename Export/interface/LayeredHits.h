#pragma once

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "FWCore/Framework/interface/Event.h"

#include "starTrack/MCTrackAnalyzer/interface/DetIdClassification.h"
#include "starTrack/MCTrackAnalyzer/interface/tHit.h"
#include "starTrack/MCTrackAnalyzer/interface/tHits.h"

#include "starTrack/MCTrackAnalyzer/interface/Event.pb.h"

#include <map>
#include <vector>
#include <fstream>

class LayeredHits {

public:

	typedef std::map<DetIDClassification, tHits> mLayeredHits;
	typedef mLayeredHits::value_type pLayeredHits;

	LayeredHits( const TrackerGeometry& tracker, std::vector<PSimHit> simHits){
		for(auto hit : simHits){
			DetId detID = hit.detUnitId();

			if(!detID)
				continue;

			DetIDClassification detClass(detID);
			const GeomDet* det = tracker.idToDet(detID);

			if(!det)
				continue;

			GlobalPoint globalPos = det->toGlobal(hit.localPosition());
			tHit mHit(globalPos.x(), globalPos.y(), globalPos.z(), detID.rawId(), hit);

			mLayeredHits::iterator it = mHits.find(detClass);
			if( it != mHits.end())
				it->second.push_back(mHit);
			else{
				tHits layerHits;
				layerHits.push_back(mHit);
				mHits.insert(pLayeredHits(detClass, layerHits));
			}
		}
	}

	LayeredHits(const TransientTrackingRecHitBuilder& recHitBuilder, std::vector<const TrackingRecHit*> recoHits,
			TrackerHitAssociator * associater = NULL){

		for(auto hit : recoHits){
			DetId detID = hit->geographicalId()();

			if(!detID)
				continue;

			DetIDClassification detClass(detID);

			GlobalPoint globalPos = recHitBuilder.build(hit)->globalPosition();
			tHit mHit(globalPos.x(), globalPos.y(), globalPos.z(), detID.rawId(), *hit, associater);

			mLayeredHits::iterator it = mHits.find(detClass);
			if( it != mHits.end())
				it->second.push_back(mHit);
			else{
				tHits layerHits;
				layerHits.push_back(mHit);
				mHits.insert(pLayeredHits(detClass, layerHits));
			}
		}
	}

	void sortLayerHitsZ(){
		for(mLayeredHits::iterator itLayer = mHits.begin(); itLayer != mHits.end(); ++ itLayer){
			itLayer->second.sortHitsZ();
		}
	}

	void sortLayerHitsPhi(){
		for(mLayeredHits::iterator itLayer = mHits.begin(); itLayer != mHits.end(); ++ itLayer){
			itLayer->second.sortHitsPhi();
		}
	}

	const tHits getLayerHits(const DetIDClassification& layer) const {
		mLayeredHits::const_iterator it = mHits.find(layer);
		if(it != mHits.end())
			return it->second;
		else
			return tHits();
	}

	const tHits::ttHits getLayerHits(const DetIDClassification::tLayers& layers) const {
		tHits::ttHits result;

		for(DetIDClassification layer : layers){
			result.push_back(getLayerHits(layer));
		}

		return result;
	}

	const tHits getLayerHitsPacked(const DetIDClassification::tLayers& layers) const {
		tHits result;

		for(DetIDClassification layer : layers){
			tHits layerHits = getLayerHits(layer);
			result.insert(result.end(), layerHits.begin(), layerHits.end());
		}

		return result;
	}

	const tHits::ttHits getBarrelHits() const {
		return getLayerHits(DetIDClassification::getBarrelLayers());
	}

	const tHits::ttHits getForwardHits() const {
		return getLayerHits(DetIDClassification::getForwardLayers());
	}

	const tHits::ttHits getBackwardHits() const {
		return getLayerHits(DetIDClassification::getBackwardLayers());
	}

	void dumpLayeredHits(PB_Event::PEventContainer & evtContainer, const edm::Event & evt, const std::vector<SimTrack>& simTracks) const {

		PB_Event::PEvent * pEvent =  evtContainer.add_events();

		pEvent->set_eventnumber( evt.id().event());
		pEvent->set_lumisection(evt.id().luminosityBlock());
		pEvent->set_runnumber(evt.id().run());

		uint hitID = 0;

		//barrel hits
		tHits::ttHits hits = getBarrelHits();
		for(uint i = 0; i < hits.size(); ++i){
			for(tHit hit : hits[i]){

				PB_Event ::PHit * pHit = pEvent->add_hits();
				pHit->set_detectorid(hit.detID);
				pHit->set_layer(i+1);
				pHit->set_hitid(hitID++);
				pHit->set_detectortype(PB_Event::DetectorType::BARREL);

				pHit->mutable_position()->set_x(hit.globalX);
				pHit->mutable_position()->set_y(hit.globalY);
				pHit->mutable_position()->set_z(hit.globalZ);

				pHit->set_simtrackid(hit.simTrackId);

				auto simTrack = std::find_if(simTracks.begin(), simTracks.end(),
						[hit] ( SimTrack const& a ) {
					return a.trackId() == hit.simTrackId;
				});

				if(simTrack != simTracks.end()){

					GlobalVector trackPt = GlobalVector(simTrack->momentum().x(),
							simTrack->momentum().y(),
							simTrack->momentum().z());
					pHit->set_simtrackpt(trackPt.perp());
				}

				/*std::cout << "x: " << pHit->position().x()<< " y: " <<  pHit->position().y()<< " z: " <<  pHit->position().z()
						  << " layer: " << pHit->layer()<< " detID: " << pHit->detectorid()<< " hitID: " << pHit->hitid();
				std::cout << " simTrackID: " << pHit->simtrackid() << " pT: " << pHit->simtrackpt();
				std::cout << std::endl;*/

			}
		}

		//forward hits
		hits = getForwardHits();
		for(uint i = 0; i < hits.size(); ++i){
			for(tHit hit : hits[i]){

				PB_Event ::PHit * pHit = pEvent->add_hits();
				pHit->set_detectorid(hit.detID);
				pHit->set_layer(i+1);
				pHit->set_hitid(hitID++);
				pHit->set_detectortype(PB_Event::DetectorType::FORWARD);

				pHit->mutable_position()->set_x(hit.globalX);
				pHit->mutable_position()->set_y(hit.globalY);
				pHit->mutable_position()->set_z(hit.globalZ);

			}
		}

		//backward hits
		hits = getBackwardHits();
		for(uint i = 0; i < hits.size(); ++i){
			for(tHit hit : hits[i]){

				PB_Event ::PHit * pHit = pEvent->add_hits();
				pHit->set_detectorid(hit.detID);
				pHit->set_layer(i+1);
				pHit->set_hitid(hitID++);
				pHit->set_detectortype(PB_Event::DetectorType::BACKWARD);

				pHit->mutable_position()->set_x(hit.globalX);
				pHit->mutable_position()->set_y(hit.globalY);
				pHit->mutable_position()->set_z(hit.globalZ);

			}
		}
	}

private:
	mLayeredHits mHits;

};

std::ostream& operator<<(std:: ostream& s, const LayeredHits& lHits);
