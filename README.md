# Optimal Control Applied to Zika Virus Epidemics in Colombia and Puerto Rico (2023). Wencel Valega Mackenzie, Karen Rios-Soto, Suzanne Lenhart

## This repository has been created to reproduce the numerical results in Valega-Mackenzie et. al, 2023. 

The folder [Optimal_Control](Optimal_Control) contains the MATLAB files to reproduce the numerical simulations for Colombia and Puerto Rico in the numerical simulation section of the paper. The [ODE45zikaOC](Optimal_Control/ode45zikaOC.m) file contains the Forward Backward Sweep method which needs the state system, [zikastates](Optimal_Control/zikastates.m) and adjoint system,  [zikaadjonts](Optimal_Control/zikaadjoints.m) (for $J_1$) or [zikaadjonts_2](Optimal_Control/zikaadjoints.m) (for $J_2$)  respectively to run. The [Figure](Optimal_Control/Figures.m) file has the parameters for each country to run all the cases and generate the figures illustrated in the paper. 

The folder [Par_estimation](Par_estimation) contains the MATLAB files to implement the parameter estimation using fmincon and multistart functions. The [zikaoc_multistart](Par_estimation/zikaoc_multistart.m) file has the fmincon and multistart routine to estimate the model parameter values for each country using their corresponding error functions in [errorfun_CO](Par_estimation/errorfun_CO.m) and [errorfun_PR](Par_estimation/errorfun_PR.m). The error functions need to approximate the solutions of the model in absence of control using the [ode45zika](Par_estimation/ode45zika.m) file. Once the parameters of the model have been estimated, run the files [Zika_Data_CO](Par_estimation/Zika_Data_CO.m) and [Zika_Data_PR](Par_estimation/Zika_Data_PR.m) to display the figure overlapping the incidence data and predicted incidence using the estimated parameters.

All MATLAB files run for MATLAB version R2021 or above.
