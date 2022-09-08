# GPCR-Landscapes

The 'mat' files can be directly loaded in MATLAB.

The document FilesAndVariables.docx has a brief description of the different file types.

WSME model outputs can be plotted for different GPCRs using the code Plot_Imp_Variables.m (input files can be modified in lines 6 and 7 of this code)

The bWSME model code used for the simulations are titled cmapCalcElecBlock.m (for reading PDB files and generating contact-map and electrostatic interaction energies), DSCcalc_Block.m (generates DSC curves at specific temperatures for a fixed set of parameters) and FesCalc_Block_full.m (calculates 1D/2D free energy profiles/landscapes, residue folding probabilities and coupling free energy matrices). Also, check https://github.com/AthiNaganathan/WSMEmodel.
