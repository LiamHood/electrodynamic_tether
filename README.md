# electrodynamic_tether
Only sim1.m, sim2.m, and encke_sim.m are designed to be run on their own. sim2 and encke_sim are near feature parity with encke_sim most up to date. These have the ability to set default values for all the initial conditions and satellite system parameters or run the file and modify individual values. sim1.m and sim2.m are using Gaussian Variation of Parameters to propagate the force of the perturbation forces of the tether. encke_sim.m uses Encke's method to propagate the orbit. The value of dr/rp in Vallado is suggested to be 1% but I have that .1% gives better values when looking at the ~30 hour range.

# Models
There are two magnetic models, a non-tilted magnetic model or the IGRF magnetic model. 
There tether is modeled as a stiff bar with mass with two point masses at either end
IRI is planned to be implemented for electron density in order to have orbit motion limited current
