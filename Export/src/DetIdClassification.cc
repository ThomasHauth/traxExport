#include "traxExport/Export/interface/DetIdClassification.h"

//stream operator
std::ostream& operator<<(std:: ostream& s, const DetIDClassification& det) {

	s << (det.isBarrel() ? "Barrel layer: " : "Endcap layer: ") << det.getLayer();
	if ( det.isBarrel() )
	{
		s << std::endl << "Ladder : " << det.getLadder();
		s << " Side : " << det.getSide();
	}
	return s;
}
