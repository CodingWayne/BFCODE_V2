#include "uclib.h"
void bodyforce(Real*, int, CoordReal*, CoordReal*, CoordReal*, CoordReal*, CoordReal*);
void USERFUNCTION_EXPORT uclib()
{
/* Register user functions here */
ucfunc(bodyforce, "RegionProfile", "bodyforce");
ucarg(bodyforce, "Cell", "Centroid", sizeof(CoordReal[3]));
ucarg(bodyforce, "Cell", "$$RelativeVelocity", sizeof(CoordReal[3]));
ucarg(bodyforce, "Cell", "$ProstarCellIndex", sizeof(CoordReal[1]));
ucarg(bodyforce, "Cell", "$Iteration", sizeof(CoordReal[1]));
ucarg(bodyforce, "Cell", "$RigidBodyAngle1Report", sizeof(CoordReal[1]));
}