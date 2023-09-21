The following repository has been created to reproduce the numerical results in Valega-Mackenzie et. al, 2023. 
## Optimal Control Applied to Zika Virus Epidemics in Colombia and Puerto Rico (2023) 

The folder [Optimal_Control](ZIKV_OC-CO-PR/Optimal_Control) contains the MATLAB files to reproduce the numerical simulations for Colombia and Puerto Rico in the numerical simulation section of the paper. The [ODE45zikaOC](Optimal_Control/ode45zikaOC.m) file contains the Forward Backward Sweep method which needs the state system, [zikastates](Optimal_Control/zikastates.m) and adjoint system,  [zikaadjonts](Optimal_Control/zikaadjoints.m) and [zikaadjonts_2](Optimal_Control/zikaadjoints.m) for $J_1$ and $J_2$ respectively to run. The [Figure](Optimal_Control/Figures.m) file has the parameters for each country to run all the cases illustrated in the paper. 

The folder [Par_estimation](ZIKV_OC-CO-PR/Par_estimation) contains the MATLAB files to implement the parameter estimation using fmincon and multistart functions. The [zikaoc_multistart](Par_estimation/zikaoc_multistart.m) file has the fmincon and multistart to estimate the incidence error 
