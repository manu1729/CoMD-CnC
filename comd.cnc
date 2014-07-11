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
[ struct atomInfo *AtomInfo];



[ int IT ] ; // number of iterations
[ int TBoxes ] ; // total number of link cells
[ int Nbs ]; // number of neighbors 


// tags declarations 
< int [2] AdvVelocity > ;
< int [3] UpdateBox > ;
< int [3] UpdateNbrs > ;
< int [2] GenForceDataTag > ;
< int [4] ForcefromNbrs >;
< int [1] GenForceTags > ;
< int [2] ComputeForce> ;
< int [2] Reduce >;


/// Step Prescriptions
< AdvVelocity > :: ( advanceVelocityStep ) ;
< UpdateBox >  :: ( updateBoxStep ) ;
< UpdateNbrs > :: ( updateNeighborsStep );
< GenForceDataTag > :: ( generateDataforForceStep);
< ForcefromNbrs > :: ( computeForcefromNeighborsStep ); 
< GenForceTags > :: ( generateForceTagsStep ) ;
< ComputeForce > :: ( forceStep ) ;
< Reduce > :: ( reduceStep );


// Input output relationships
[ B : i, 0, 0, iter ] -> ( advanceVelocityStep : i, iter ) -> [ B : i, 1, 0, iter];
( advanceVelocityStep : i, iter ) -> < UpdateBox : 0, 0, iter > ; // generated when i == 0, remaning UpdateBox tags are generated by updateBoxStep

[ B : i, 1, 0, iter ]-> ( updateBoxStep : i, k, iter ) -> < UpdateNbrs : i, j, iter > ; // only tag for the first neighbor step is generated here 
( updateBoxStep : i, k, iter ) -> [ AtomInfo : i, j, iter ];
( updateBoxStep : i, k, iter ) -> < GenForceDataTag : 0, iter>;  // generated by the last box

[ AtomInfo : i, j, iter ], [ B : j, 1, 0, iter ] -> ( updateNeighborsStep : i, j, iter ) -> < UpdateNbrs : i, k , iter >; // generate tag for the next neighbor step
( updateNeighborsStep : i, j, iter ) -> [ AtomInfo : i, k , iter ]; // "B[j...]" is updated but not inserted as an item, okay here as the updates are sequential
( updateNeighborsStep : i, j, iter ) -> < UpdateBox : i+1, 0, iter > ; // executed by the "last" neighbor step 

[ B : i, 1, 0, iter ]-> (generateDataforForceStep: i, iter) -> [ B : i, 3, 0, iter ], < GenForceDataTag : i+1, iter>;
(generateDataforForceStep: i, iter) -> < GenForceTags: iter >; // when i == last box

//[ B : i, 2, 0, iter ] -> ( genForceData : i, iter ) -> [ B : i, 3, 0, iter ];
//( genForceData : i, iter ) -> [ B : i-1, 2, 0, iter ], < GenForceDataTag : i-1, iter> ; // this is executed when i >= 0 


[ TBoxes : 0 ], [ Nbs : 0 ] -> ( generateForceTagsStep : iter ) -> < ComputeForce: { 0 .. TBoxes[0]-1 }, iter > ;  


[ B : i, 3, 0, iter ] -> ( forceStep : i, iter ) -> < ForcefromNbrs : i, j, k, iter >; // generate tag for force computation due to the first neighbor 
[ B : i, 3, k, iter ], [ B : j, 3, 0, iter ] -> ( computeForcefromNeighborsStep: i, j, k, iter ) -> [ B : i, 3, k+1, iter ], < ForcefromNbrs : i, jnext, k+1, iter >;
( computeForcefromNeighborsStep: i, j, k, iter ) -> [ B : i, 4, 0, iter ]; // when k == 27-1


[ B : i, 4, 0, iter ], [redc : i, iter ], [ IT : 0 ], [ TBoxes : 0 ] -> ( reduceStep : i, iter ) -> [ redc : i+1, iter ], [ B : i, 0, 0, iter+1 ] ;
( reduceStep : i, iter ) -> [ B : i, 5, 0, iter ] ; // is executed only when i == 1727 and iter == max_iteration
( reduceStep : i, iter ) -> [ redc : 0, iter+1 ] ;  // when i = TBoxes[0] and iter < IT[0]

//////////////////////////////////////////////////////// Done till here //////////////////////////////////////////////////////////////////////




// Write graph inputs and start steps
env -> [ B : i, 0, 0, 1 ], [ s : 0, 1 ], [ IT : 0 ], [ TBoxes : 0 ], [ Nbs : 0 ]; //[ TAtoms : 0 ], [ AperB : 0 ], 
env -> < AdvVelocity : i, 1 >, < GenForceTags : 1 >, < Reduce: i, IT[0] >;
env -> [ SF : 1 ], [ POT : 1 ], [ SPECIES : 1 ], [ DD : 1 ], [ LC : 1], [ NAtoms : 1 ], [ ATOMS : 1], [ redc : 0, 1];

// Return outputs to the caller
[ B : i, 5, 0, iter  ], [ redc : i+1, iter] -> ( env: i, iter );
