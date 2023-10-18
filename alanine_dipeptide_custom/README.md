# WESTPA with Alanine Dipeptide in implicit solvent

This example demonstrates running WESTPA using alanine dipeptide in implicit solvent. The system is parameterized in ff14SB-onlysc in implicit solvent model igb=1. The custom distance scheme is enabled in this version such that segments are organized pair-wise based on similarity and merged based on that hierarchy, as indicated in sort.py. You will have to use the `custom_order_resampler` branch from https://github.com/jeremyleung/westpa to use sort.py. The lack of "80" in the name indicate that the gamma_ln collision frequency is set to 5 ps^{-1}, a value used when there is explicit water modeled. 

The WE Tau chosen is 100ps. The starting structure is a conformation found in C7eq well of the Phi/Psi Ramachandran plot.

