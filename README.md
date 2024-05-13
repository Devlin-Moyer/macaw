this is the roughest possible draft of a README; pls do not judge

Python package dependencies: cobra, networkx, pebble

pathway of 0 means reaction wasn't in any pathway (and probably not flagged by any tests except maybe the diphosphate test)

## Dead-End Test
identify and all dead-end reactions and reactions that are nominally reversible but only actually capable of sustaining flux in a single direction

## Dilution Test
identify reactions that become incapable of sustaining non-zero fluxes when dilution constraints are added to the model

## Diphosphate Test
identify all reversible reactions that involve diphosphate that aren't transporting it between compartments cuz they should probably be irreversible (only if lists of IDs of metabolites representing diphosphate and phosphate were provided)

## Duplicate Test
identify sets of reactions that are potentially duplicates of each other


## Loop Test
identify reactions that are capable of sustaining non-zero fluxes when all exchange reactions are blocked

## Pathways
combine the edge lists and add edges between reactions flagged by different tests. also add a column to the test results indicating which "pathway" (connected component of the network) each reaction winds up in
