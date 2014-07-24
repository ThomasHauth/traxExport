#pragma once

class HitUtilities {
public:
	static GlobalPoint trackingRecHitToGlobal(TrackingRecHitRef const& hit,
			edm::ESHandle<TrackerGeometry> const& tracker) {
		const GeomDet* hitDet = tracker->idToDet(hit->geographicalId()());
		assert(hitDet != nullptr);
		return hitDet->toGlobal(hit->localPosition());
	}

	static GlobalPoint trackingRecHitToGlobal(TrackingRecHit const* hit,
			edm::ESHandle<TrackerGeometry> const& tracker) {
		const GeomDet* hitDet = tracker->idToDet(hit->geographicalId()());
		assert(hitDet != nullptr);
		return hitDet->toGlobal(hit->localPosition());
	}
};
