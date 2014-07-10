// Author: Manu Shantharam (mshantharam@ucsd.edu)

// item collection declarations
[ struct box *B ];   // represents link cell
[ struct info *UI];  // update info for each cell
[ int count ];  // keeps count of number of cells updating a particular cell
[ int *s ];  // this item is used for serialization - kind of barrier
[ SimFlat *SF];  // Contains data related to the simulation system
[ LjPotential *POT];
[ SpeciesData *SPECIES];
[ Domain *DD ];
[ LinkCell *LC ];
[ int *NAtoms ];
[ Atoms *ATOMS];
[ struct myReduction *redc]; // used for reduction



[ int IT ] ; // number of iterations
[ int TBoxes ] ; // total number of link cells
[ int Nbs ]; // number of neighbors 


// tags declarations 
< int [2] AdvVelocity > ;
//< int [2] AdvPosition > ;
< int [3] UpdateBox > ;
< int [2] GenForceTags > ;
< int [2] ComputeForce> ;
< int [2] Reduce >;


// Step Prescriptions
< AdvVelocity > :: ( advanceVelocityStep ) ;
//< AdvPosition > :: ( advancePositionStep ) ;
< UpdateBox >  :: ( updateBoxStep ) ;
< GenForceTags > :: ( generateForceTagsStep ) ;
< ComputeForce > :: ( forceStep ) ;
< Reduce > :: ( reduceStep );


// Input output relationships
[ B : i/2, 0, 0, iter ], [ SF : 1 ] -> ( advanceVelocityStep : i, iter ) -> [ B : i, 1, 0, iter], < UpdateBox : i, 0, iter > ;  

[ B : i, 1, 0, iter ], [ SF : 1 ], [ s : i, iter ], [ IT : 0 ], [ TBoxes : 0 ] -> ( updateBoxStep : i, k, iter ) -> [ B : i, 3, 0, iter ],  [ s : i+1, iter ]; 
( updateBoxStep : i, k, iter ) -> < GenForceTags : TBoxes[0], iter+1>;

[ s : N, iter ], [ TBoxes : 0 ], [ Nbs : 0 ] -> ( generateForceTagsStep : N, iter ) -> < ComputeForce: { 0 .. TBoxes[0]-1 }, iter > ;  

[ B : i, 3, 0, iter ], [ SF : 1 ], [ IT : 0 ] -> ( forceStep : i, iter ) -> [ B : i, 4, 0, iter];
//( forceStep : i, iter ) -> [ B : i, 4, 0, iter ];  


[ B : i, 4, 0, iter ], [redc : i, iter ], [ IT : 0 ], [ TBoxes : 0 ] -> ( reduceStep : i, iter ) -> [ redc : i+1, iter ], [ B : i, 0, 0, iter+1 ] ;
( reduceStep : i, iter ) -> [ B : i, 5, 0, iter ] ; // is executed only when i == 1727 and iter == max_iteration
( reduceStep : i, iter ) -> [ redc : 0, iter+1 ] ;  // when i = TBoxes[0] and iter < IT[0]

// Write graph inputs and start steps
env -> [ B : i, 0, 0, 1 ], [ s : 0, 1 ], [ IT : 0 ], [ TBoxes : 0 ], [ Nbs : 0 ]; //[ TAtoms : 0 ], [ AperB : 0 ], 
env -> < AdvVelocity : i, 1 >, < GenForceTags : TBoxes[0], 1 >, < Reduce: i, IT[0] >;
env -> [ SF : 1 ], [ POT : 1 ], [ SPECIES : 1 ], [ DD : 1 ], [ LC : 1], [ NAtoms : 1 ], [ ATOMS : 1], [ redc : 0, 1];

// Return outputs to the caller
[ B : i, 5, 0, iter  ], [ redc : i+1, iter] -> ( env: i, iter );
