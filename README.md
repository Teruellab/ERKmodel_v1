# ERKmodel_v1
Implementation of an ODE ERK model from Sturm et al, Science Signaling 2010 in MATLAB.

To use:
Download and unpack .zip file. Navigate to this folder in your open session of MATLAB, or add it to your MATLAB path.
A simple script you can look at to get started is "sim_simple.m", which simulates and plots a basic dose response in
control conditions.


Structure of files:
- "sim_" files are scripts. You can open and run these to simulate the model in different conditions (input dose, 
    ERK and MEK level, etc.)
- "graph_" files are scripts that allow for exploration of simulation data
_ "erkInitialize.m", "erkODE.m", and "erkSimulate.m" are the core model functions - you should not need to modify these
  directly, but they set parameters and initial values, define ODE mass-action equations, perform basal equilibration,
  and simulation.
- remaining files are related to spreadsheet update functionality: if you want to change default model equations or
  parameters, you can make the changes in "ERK reactions.xlsx", then run (in the MATLAB command line)
  
  updateModel('Erk reactions.xlsx'); % this will update the core model files.
  Note: there are known issues with this functionality for versions of MATLAB prior to R2015b.
