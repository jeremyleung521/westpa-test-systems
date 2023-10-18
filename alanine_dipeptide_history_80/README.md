# WESTPA with Alanine Dipeptide in implicit solvent

This example demonstrates running WESTPA using alanine dipeptide in implicit solvent. The system is parameterized in ff14SB-onlysc in implicit solvent model igb=1. The history 10 scheme is enabled in this version such that only segments with shared history (same parent trajectory 10 iterations ago) within a bin are allowed to be merged together. The "80" in the name indicate that the gamma_ln collision frequency is set to 80 ps^{-1} to match the viscosity of water. 

The WE Tau chosen is 100ps. The starting structure is a conformation found in C7eq well of the Phi/Psi Ramachandran plot.

