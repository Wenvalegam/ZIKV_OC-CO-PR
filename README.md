# Optimal Control Applied to Zika Virus Epidemics in Colombia and Puerto Rico (2023)
Authors: Wencel Valega Mackenzie, Karen Rios-Soto, Suzanne Lenhart

## Overview
This repository contains code to reproduce the numerical results presented in our paper (citation goes here). It provides MATLAB files and instructions for simulating and estimating parameters related to Zika virus epidemics in Colombia and Puerto Rico.

## Folder Structure
- **[Optimal_Control](Optimal_Control)**: Contains MATLAB files for numerical simulations.
  - [ode45zikaOC.m](Optimal_Control/ode45zikaOC.m): Implements the Forward Backward Sweep method.
  - [zikastates.m](Optimal_Control/zikastates.m): Defines the state system.
  - [zikaadjoints.m](Optimal_Control/zikaadjoints.m): Defines the adjoint system for J1.
  - [zikaadjoints_2.m](Optimal_Control/zikaadjoints_2.m): Defines the adjoint system for J2.
  - [Figures.m](Optimal_Control/Figures.m): Contains parameters for generating figures with controls.

- **[Par_estimation](Par_estimation)**: Contains MATLAB files for parameter estimation.
  - [zikaoc_multistart.m](Par_estimation/zikaoc_multistart.m): Implements parameter estimation using fmincon and multistart.
  - [errorfun_CO.m](Par_estimation/errorfun_CO.m): Error function for parameter estimation in Colombia.
  - [errorfun_PR.m](Par_estimation/errorfun_PR.m): Error function for parameter estimation in Puerto Rico.
  - [ode45zika.m](Par_estimation/ode45zika.m): Approximates solutions of the model without controls.
  - [Zika_Data_CO.m](Par_estimation/Zika_Data_CO.m): Displays figures for incidence data and predictions in Colombia.
  - [Zika_Data_PR.m](Par_estimation/Zika_Data_PR.m): Displays figures for incidence data and predictions in Puerto Rico.

## Requirements
- MATLAB R2021 or above.
- MATLAB Optimization Toolbox

## Citation
If you use this code in your research, please cite our paper [citation goes here].
