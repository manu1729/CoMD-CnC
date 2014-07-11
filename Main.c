/// \file
/// Main program
///
/// \mainpage CoMD: A Classical Molecular Dynamics Mini-app
///
/// CoMD is a reference implementation of typical classical molecular
/// dynamics algorithms and workloads.  It is created and maintained by
/// The Exascale Co-Design Center for Materials in Extreme Environments
/// (ExMatEx).  http://codesign.lanl.gov/projects/exmatex.  The
/// code is intended to serve as a vehicle for co-design by allowing
/// others to extend and/or reimplement it as needed to test performance of 
/// new architectures, programming models, etc.
///
/// The current version of CoMD is available from:
/// http://exmatex.github.io/CoMD
///
/// To contact the developers of CoMD send email to: exmatex-comd@llnl.gov.
///
/// Table of Contents
/// =================
///
/// Click on the links below to browse the CoMD documentation.
///
/// \subpage pg_openmp_specifics
///
/// \subpage pg_md_basics
///
/// \subpage pg_building_comd
///
/// \subpage pg_running_comd
///
/// \subpage pg_measuring_performance
///
/// \subpage pg_problem_selection_and_scaling
///
/// \subpage pg_verifying_correctness
///
/// \subpage pg_comd_architecture
///
/// \subpage pg_optimization_targets
///
/// \subpage pg_whats_new


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <assert.h>


#include "CoMDTypes.h"
#include "decomposition.h"
#include "linkCells.h"
#include "initAtoms.h"
#include "memUtils.h"
#include "yamlOutput.h"
#include "mycommand.h"
#include "parallel.h"
#include "ljForce.h"
#include "timestep.h"
#include "haloExchange.h"


#include "Context.h"


struct DomainSt* initDecomposition(int xproc, int yproc, int zproc,
                                   real3 globalExtent, Context *c);

LinkCell* initLinkCells(const struct DomainSt* domain, real_t cutoff, Context *c);
Atoms* initAtoms(struct LinkCellSt* boxes, Context *c);

static void printSimulationDataYaml(FILE* file, SimFlat* s);

void initSubsystems(void) {
#if REDIRECT_OUTPUT
   freopen("testOut.txt","w",screenOut);
#endif

   yamlBegin();
}

/// decide whether to get LJ or EAM potentials
BasePotential* initPotential(
   int doeam, const char* potDir, const char* potName, const char* potType, Context *context)
{
   BasePotential* pot = NULL;

   pot = initLjPot(context);
   assert(pot);
   return pot;
}

SpeciesData* initSpecies(BasePotential* pot, Context *context)
{
 //  SpeciesData* species = comdMalloc(sizeof(SpeciesData));
   SpeciesData* species;
   cncHandle_t sp_handle = cncCreateItem_SPECIES(&species,sizeof(SpeciesData));

   cncPut_SPECIES(sp_handle, 1, context);

   strcpy(species->name, pot->name);
   species->atomicNo = pot->atomicNo;
   species->mass = pot->mass;

   return species;
}

void sumAtoms(SimFlat* s)
{
   // sum atoms across all processers
   s->atoms->nLocal = 0;
   for (int i = 0; i < s->boxes->nLocalBoxes; i++)
   {
      s->atoms->nLocal += s->boxes->nAtoms[i];
   }

//   startTimer(commReduceTimer);
   addIntParallel(&s->atoms->nLocal, &s->atoms->nGlobal, 1);
//   stopTimer(commReduceTimer);
}


Validate* initValidate(SimFlat* sim)
{
   sumAtoms(sim);    /////////////////////////////////ToDo Manu: need to change
   Validate* val = comdMalloc(sizeof(Validate));
   val->eTot0 = (sim->ePotential + sim->eKinetic) / sim->atoms->nGlobal;
   val->nAtoms0 = sim->atoms->nGlobal;

   if (printRank())
   {
      fprintf(screenOut, "\n");
      printSeparator(screenOut);
      fprintf(screenOut, "Initial energy : %14.12f, atom count : %d \n",
            val->eTot0, val->nAtoms0);
      fprintf(screenOut, "\n");
   }
   return val;
}

void validateResult(const Validate* val, SimFlat* sim)
{
   if (printRank())
   {
      real_t eFinal = (sim->ePotential + sim->eKinetic) / sim->atoms->nGlobal;

      int nAtomsDelta = (sim->atoms->nGlobal - val->nAtoms0);

      fprintf(screenOut, "\n");
      fprintf(screenOut, "\n");
      fprintf(screenOut, "Simulation Validation:\n");

      fprintf(screenOut, "  Initial energy  : %14.12f\n", val->eTot0);
      fprintf(screenOut, "  Final energy    : %14.12f\n", eFinal);
      fprintf(screenOut, "  eFinal/eInitial : %f\n", eFinal/val->eTot0);
      if ( nAtomsDelta == 0)
      {
         fprintf(screenOut, "  Final atom count : %d, no atoms lost\n",
               sim->atoms->nGlobal);
      }
      else
      {
         fprintf(screenOut, "#############################\n");
         fprintf(screenOut, "# WARNING: %6d atoms lost #\n", nAtomsDelta);
         fprintf(screenOut, "#############################\n");
      }
   }
}


/// Check that the user input meets certain criteria.
void sanityChecks(Command cmd, double cutoff, double latticeConst, char latticeType[8]) {
   int failCode = 0;

   // Check that domain grid matches number of ranks. (fail code 1)
   int nProcs = cmd.xproc * cmd.yproc * cmd.zproc;
   if (nProcs != getNRanks())
   {
      failCode |= 1;
      if (printRank() )
         fprintf(screenOut,
                 "\nNumber of MPI ranks must match xproc * yproc * zproc\n");
   }

   // Check whether simuation is too small (fail code 2)
   double minx = 2*cutoff*cmd.xproc;
   double miny = 2*cutoff*cmd.yproc;
   double minz = 2*cutoff*cmd.zproc;
   double sizex = cmd.nx*latticeConst;
   double sizey = cmd.ny*latticeConst;
   double sizez = cmd.nz*latticeConst;

   if ( sizex < minx || sizey < miny || sizez < minz)
   {
      failCode |= 2;
      if (printRank())
         fprintf(screenOut,"\nSimulation too small.\n"
                 "  Increase the number of unit cells to make the simulation\n"
                 "  at least (%3.2f, %3.2f. %3.2f) Ansgstroms in size\n",
                 minx, miny, minz);
   }

   // Check for supported lattice structure (fail code 4)
   if (strcasecmp(latticeType, "FCC") != 0)
   {
      failCode |= 4;
      if ( printRank() )
         fprintf(screenOut,
                 "\nOnly FCC Lattice type supported, not %s. Fatal Error.\n",
                 latticeType);
   }
   int checkCode = failCode;
//   bcastParallel(&checkCode, sizeof(int), 0);  // ToDo: Manu
   // This assertion can only fail if different tasks failed different
   // sanity checks.  That should not be possible.
   assert(checkCode == failCode);

   if (failCode != 0)
      exit(failCode);
}

/// Print information about the simulation in a format that is (mostly)
/// YAML compliant.
void printSimulationDataYaml(FILE* file, SimFlat* s)
{
   // All ranks get maxOccupancy
   int maxOcc = maxOccupancy(s->boxes);

   // Only rank 0 prints
   if (! printRank())
      return;

   fprintf(file,"Simulation data: \n");
   fprintf(file,"  Total atoms        : %d\n",
           s->atoms->nGlobal);
   fprintf(file,"  Min global bounds  : [ %14.10f, %14.10f, %14.10f ]\n",
           s->domain->globalMin[0], s->domain->globalMin[1], s->domain->globalMin[2]);
   fprintf(file,"  Max global bounds  : [ %14.10f, %14.10f, %14.10f ]\n",
           s->domain->globalMax[0], s->domain->globalMax[1], s->domain->globalMax[2]);
   printSeparator(file);
   fprintf(file,"Decomposition data: \n");
   fprintf(file,"  Processors         : %6d,%6d,%6d\n",
           s->domain->procGrid[0], s->domain->procGrid[1], s->domain->procGrid[2]);
   fprintf(file,"  Local boxes        : %6d,%6d,%6d = %8d\n",
           s->boxes->gridSize[0], s->boxes->gridSize[1], s->boxes->gridSize[2],
           s->boxes->gridSize[0]*s->boxes->gridSize[1]*s->boxes->gridSize[2]);
   fprintf(file,"  Box size           : [ %14.10f, %14.10f, %14.10f ]\n",
           s->boxes->boxSize[0], s->boxes->boxSize[1], s->boxes->boxSize[2]);
   fprintf(file,"  Box factor         : [ %14.10f, %14.10f, %14.10f ] \n",
           s->boxes->boxSize[0]/s->pot->cutoff,
           s->boxes->boxSize[1]/s->pot->cutoff,
           s->boxes->boxSize[2]/s->pot->cutoff);
   fprintf(file, "  Max Link Cell Occupancy: %d of %d\n",
           maxOcc, MAXATOMS);
   printSeparator(file);
   fprintf(file,"Potential data: \n");
//   s->pot->print(file, s->pot);

   fflush(file);
}

SimFlat* initSimulationNew(Command cmd, Context *context) {
 //  SimFlat* sim = comdMalloc(sizeof(SimFlat));


   SimFlat *sim;
   cncHandle_t sf_handle = cncCreateItem_SF(&sim,sizeof(SimFlat));

   cncPut_SF(sf_handle, 1, context);

   sim->nSteps = cmd.nSteps;
   sim->printRate = cmd.printRate;
   sim->dt = cmd.dt;
   sim->domain = NULL;
   sim->boxes = NULL;
   sim->atoms = NULL;
   sim->ePotential = 0.0;
   sim->eKinetic = 0.0;
   sim->atomExchange = NULL;


   sim->pot = initPotential(cmd.doeam, cmd.potDir, cmd.potName, cmd.potType, context);
   real_t latticeConstant = cmd.lat;
   if (cmd.lat < 0.0)
      latticeConstant = sim->pot->lat;

   // ensure input parameters make sense.
   sanityChecks(cmd, sim->pot->cutoff, latticeConstant, sim->pot->latticeType);   /////////////////////////////////ToDo Manu: need to uncomment this

   sim->species = initSpecies(sim->pot, context);

   real3 globalExtent;
   globalExtent[0] = cmd.nx * latticeConstant;
   globalExtent[1] = cmd.ny * latticeConstant;
   globalExtent[2] = cmd.nz * latticeConstant;

   sim->domain = initDecomposition(
      cmd.xproc, cmd.yproc, cmd.zproc, globalExtent, context);

   sim->boxes = initLinkCells(sim->domain, sim->pot->cutoff, context);
   sim->atoms = initAtoms(sim->boxes, context);

   // create lattice with desired temperature and displacement.
   createFccLattice(cmd.nx, cmd.ny, cmd.nz, latticeConstant, sim);
   setTemperature(sim, cmd.temperature);
   randomDisplacements(sim, cmd.initialDelta);

   sim->atomExchange = initAtomHaloExchange(sim->domain, sim->boxes);

   // Forces must be computed before we call the time stepper.
 //  startTimer(redistributeTimer);
   redistributeAtoms(sim);          /////////////////////////////////ToDo Manu: need to uncomment this
 //  stopTimer(redistributeTimer);

//   startTimer(computeForceTimer);
   computeForce(sim);              /////////////////////////////////ToDo Manu: need to uncomment this
//   stopTimer(computeForceTimer);

   kineticEnergy(sim);             /////////////////////////////////ToDo Manu: need to uncomment this

   printf("======== %d, %u\n", sim->nSteps, sim->atoms);

   return sim;
}



void cncEnvIn(int argc, char** argv, Context *context) {

    int totalBoxes;
    int MAXIT;
    int maxAtomsPerBox;
    int numNbrs;
    int numAtoms;

    PRINTF("Initialization...\n");
    initSubsystems();

    Command cmd = parseCommandLine(argc, argv);
    printCmdYaml(yamlFile, &cmd);
    printCmdYaml(screenOut, &cmd);

    SimFlat* sim = initSimulationNew(cmd, context);
    printSimulationDataYaml(yamlFile, sim);
    printSimulationDataYaml(screenOut, sim);

    printf("Initial validation ... \n");

    Validate* validate = initValidate(sim);

    // This is the CoMD main loop
    const int nSteps = sim->nSteps;
    const int printRate = sim->printRate;
    int iStep = 0;


    /////////////////
 //   totalBoxes = sim->boxes->nLocalBoxes;
 //   MAXIT = nSteps;
    totalBoxes = 1728;
    printf("Number of iterations == %d\n", nSteps);
    MAXIT = nSteps+1;
    /////////////////
    maxAtomsPerBox = 60;
    numNbrs = 27;
    numAtoms = 32000;

    PRINTF("Starting computation\n");



    int i,j;
    struct box *b;
    for (i=0; i<totalBoxes; i++){
        cncHandle_t db_handle = cncCreateItem_B(&b, sizeof(struct box));
        b->i = i;
        b->ePot = 0.0;
        b->eKin = 0.0;
        cncPut_B(db_handle, i, 0, 0, 1, context);
        cncPrescribe_advanceVelocityStep(i, 1, context);
        cncPrescribe_reduceStep(i, 1, context);
    }

    int *s;
    cncHandle_t s_handle = cncCreateItem_s(&s, sizeof(int));
    *s = 0;
    cncPut_s(s_handle, 0, 1, context);

    struct myReduction *rd;
    cncHandle_t rd_handle = cncCreateItem_redc(&rd, sizeof(int));
    rd->i = 0;
    rd->ePot = 0.0;
    rd->eKin = 0.0;
    cncPut_redc(rd_handle, 0, 1, context);

    // creating IT
    int *t;
    cncHandle_t t_handle = cncCreateItem_IT(&t);
    *t = MAXIT;
    cncPut_IT(t_handle, 0, context);

    // creating TBoxes
    t_handle = cncCreateItem_TBoxes(&t);
    *t = totalBoxes;
    cncPut_TBoxes(t_handle, 0, context);

   // cncPrescribe_generateForceTagsStep(1, context);


    // creating Nbs
    t_handle = cncCreateItem_Nbs(&t);
    *t = numNbrs;
    cncPut_Nbs(t_handle, 0, context);

    printf("     Total Energy     Potential Energy    Kinetic Energy \n");
    cncPrescribe_cncEnvOut(totalBoxes-1, MAXIT-1, context);
}


void cncEnvOut(int i, int iter, BItem B, redcItem r, Context *context) {
  //  PRINTF("Box %d, iteration %d ends\n", i, iter);
    real_t p,k,t;
    p = r.item->ePot/32000;
    k = r.item->eKin/32000;
    t = p+k;

    PRINTF("Total energy = %18.12f, Potential energy = %18.12f, Kinetic energy = %18.12f\n",t,p,k);
}


