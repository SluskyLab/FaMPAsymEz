Full atom binned implementation of asymetric EZ potential from Schramm et al 2012 at 10.1016/j.str.2012.03.016 for Cbeta and Cgamma energy terms. 

AsymEZPotGenerator.py: generates bin values for z range -30 to 30 and the derivatives. These output files are used as e-tables for EX potential energies in Rosetta. 

2012SamishData/: contains publication data from Table S1

EnergyGraphs/: outputs graphs for each residue's energy term values

In addition to the paper's method, we linearly interpolate between several points, especially those with a left-hand gaussian


