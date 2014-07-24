#include "traxExport/Export/interface/tHit.h"


std::ostream& operator<<(std:: ostream& s, const tHit& hit) {
	return s << "X: " << hit.globalX << " Y: " << hit.globalY << " Z: " << hit.globalZ;
}
