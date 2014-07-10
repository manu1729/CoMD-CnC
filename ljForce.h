/// \file
/// Computes forces for the 12-6 Lennard Jones (LJ) potential.

#ifndef _LJTYPES_H_
#define _LJTYPES_H_

#ifndef _CONTEXT
#include "Context.h"
#endif

struct BasePotentialSt;
struct BasePotentialSt* initLjPot(Context *context);

#endif

