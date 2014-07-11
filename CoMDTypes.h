/// \file
/// CoMD data structures.

#ifndef __COMDTYPES_H_
#define __COMDTYPES_H_



struct SimFlatSt;
struct box;
struct info;

#include <stdio.h>
#include "mytype.h"
#include "decomposition.h"
#include "linkCells.h"
#include "initAtoms.h"
#include "haloExchange.h"


////////////////// declaring all the structures for the stub - start //////////////////
struct box {
    int i;  // box number
    double ePot;
    double eKin;
    int n;  // total number of data elements n < 10
};

struct info {
    int n;  // number of atoms
//    int data[10];
    int from;  // from box
    int to; // to box
};

struct myReduction {
    int i;
    double ePot;
    double eKin;
};

struct atomInfo {
    int id;
};

struct info allUI[1728*27];

////////////////// declaring all the structures for the stub - end //////////////////


/// The base struct from which all potentials derive.  Think of this as an
/// abstract base class.
///
/// CoMD uses the following units:
///  - distance is in Angstroms
///  - energy is in eV
///  - time is in fs
///  - force in in eV/Angstrom
///
///  The choice of distance, energy, and time units means that the unit
///  of mass is eV*fs^2/Angstrom^2.  Hence, we must convert masses that
///  are input in AMU (atomic mass units) into internal mass units.
typedef struct BasePotentialSt
{
   real_t cutoff;          //!< potential cutoff distance in Angstroms
   real_t mass;            //!< mass of atoms in intenal units
   real_t lat;             //!< lattice spacing (angs) of unit cell
   char latticeType[8];    //!< lattice type, e.g. FCC, BCC, etc.
   char  name[3];      //!< element name
   int   atomicNo;     //!< atomic number
   int  (*force)(struct SimFlatSt* s); //!< function pointer to force routine
} BasePotential;




/// species data: chosen to match the data found in the setfl/funcfl files
typedef struct SpeciesDataSt
{
   char  name[3];   //!< element name
   int   atomicNo;  //!< atomic number
   real_t mass;     //!< mass in internal units
} SpeciesData;

/// Simple struct to store the initial energy and number of atoms.
/// Used to check energy conservation and atom conservation.
typedef struct ValidateSt
{
   double eTot0; //<! Initial total energy
   int nAtoms0;  //<! Initial global number of atoms
} Validate;

///
/// The fundamental simulation data structure with MAXATOMS in every
/// link cell.
///
typedef struct SimFlatSt
{
   int nSteps;            //<! number of time steps to run
   int printRate;         //<! number of steps between output
   double dt;             //<! time step

   Domain* domain;        //<! domain decomposition data

   LinkCell* boxes;       //<! link-cell data

   Atoms* atoms;          //<! atom data (positions, momenta, ...)

   SpeciesData* species;  //<! species data (per species, not per atom)

   real_t ePotential;     //!< the total potential energy of the system
   real_t eKinetic;       //!< the total kinetic energy of the system
   HaloExchange* atomExchange;

   BasePotential *pot;    //!< the potential
} SimFlat;


typedef struct LjPotentialSt
{
   real_t cutoff;          //!< potential cutoff distance in Angstroms
   real_t mass;            //!< mass of atoms in intenal units
   real_t lat;             //!< lattice spacing (angs) of unit cell
   char latticeType[8];    //!< lattice type, e.g. FCC, BCC, etc.
   char  name[3];      //!< element name
   int   atomicNo;     //!< atomic number
   int  (*force)(SimFlat* s); //!< function pointer to force routine
   void (*print)(FILE* file, BasePotential* pot);
   void (*destroy)(BasePotential** pot); //!< destruction of the potential
   real_t sigma;
   real_t epsilon;
} LjPotential;








#endif
