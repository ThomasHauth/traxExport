#pragma once

#include <iostream>

#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

class tHit {

public:

	tHit(double x, double y, double z, unsigned int det)
	: globalX(x), globalY(y), globalZ(z), detID(det), simTrackId(0), isSim(false), simHit(PSimHit()) {}

	tHit(double x, double y, double z, unsigned int det, const PSimHit& iSimHit)
	: globalX(x), globalY(y), globalZ(z), detID(det), isSim(true), simHit(iSimHit) {

		simTrackId = iSimHit.trackId();

	}

	tHit(double x, double y, double z, unsigned int det,
			const TrackingRecHit & iRecHit, TrackerHitAssociator * associater = NULL)
	: globalX(x), globalY(y), globalZ(z), detID(det), isSim(true), simHit(PSimHit()) {

		simTrackId = 0;
		if(associater != NULL){
			std::vector<SimHitIdpr>  trackIds = associater->associateHitId(iRecHit);
			if(trackIds.size() > 0)
				simTrackId = trackIds[0].first; //if more than one track associated -> ignore it (happens seldom)
		}

	}

	tHit(const tHit& hit)
	: globalX(hit.globalX), globalY(hit.globalY), globalZ(hit.globalZ), detID(hit.detID), simTrackId(hit.simTrackId), isSim(hit.isSimulated()),
	  simHit(hit.isSimulated() ? hit.getSimHit() : PSimHit()) {}

	tHit& operator=( const tHit& rhs ) {
		globalX = rhs.globalX;
		globalY = rhs.globalY;
		globalZ = rhs.globalZ;
		detID = rhs.detID;
		isSim = rhs.isSimulated();
		simHit = isSim ? rhs.getSimHit() : PSimHit();
		return *this;
	}

	double phi() const {
		return std::atan2(globalY, globalX);
	}

	bool isSimulated() const {
		return isSim;
	}

	const PSimHit& getSimHit() const {
		assert(isSimulated());

		return simHit;
	}

	double globalX;
	double globalY;
	double globalZ;
	unsigned int detID;
	unsigned int simTrackId;

	bool isSim;
	PSimHit simHit;

};

std::ostream& operator<<(std:: ostream& s, const tHit& hit);
