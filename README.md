# Satellite ADCS LEOP Simulation

Group project designing and assessing ADCS architectures for Earth-observations microsatellite LEOP applications.

## Overview

This repository contains a simulation study developed as part of a university satellite attitude dynamics project.

The work is organized in 5 parts:
1. Main MATLAB script 'satellite_ADCS_LEOP_simulation.m' collecting spacecraft, environmental model properties and ensuring simulation phases continuity
2. Satellite ADCS simulation of detumbling phase after uncontrolled detachment from launcher platform
3. Satellite ADCS simulation of slew phase, transitioning from detumbling to pointing mode
4. Satellite ADCS simulation of pointing phase for Earth-observations nadir pointing
5. Satellite ADCS simulation of pointing phase for Earth-observations nadir pointing, using alternative reaction-wheels control

## Results and Validation

Key results:
- Attitude determination through analytical solution of Wabha's Problem provides sufficiently adequate estimate for detumbling and slew phases. PD closed-loop control successfully reduces angular velocity to ~ 1e-4 rad/s in 2500 seconds.
- LQR + State Estimator + Thruster Hysteresis logic capable of maintaining nadir pointing precision of ~ 1 deg, effectively counteracting environmental disturbances.

Representative outputs:
- Control parameters histories across different phases
- Control torque histories across different phases
- Nadir pointing performance during pointing phase

Representative figures are available in 'results/' (PNG format).
See 'results/results.txt' for figure-by-figure explanations.
The full methodology and results are documented in 'docs/report.pdf'.

## Repository structure

- 'src/' - MATLAB, Simulink implementations of each study
- 'docs/' - Project report
- 'results/' - Key result figures (PNG)

## Reproducibility and external dependencies

Run the code section by section in the given order. When prompted, run appropriate Simulink simulation e.g.: when 'execute detumble', run 'project_detumble.slx'
Once the simulation is completed, advance to next section to visualize related plots.

Some scripts call proprietary helper functions provided by the course staff, not included in this repository.

Missing dependencies:
- 'astroConstants.m' (used in 'satellite_ADCS_LEOP_simulation.m' returns a row-vector of astrodynamics-related physical constants, used for Earth gravitational parameter, mean Earth radius and mean Sun-Earth distance)

Some scripts call external MATLAB function obtained from a public online source, not included in this repository.

Missing dependencies:
- 'icosphere.m' (used in 'satellite_ADCS_LEOP_simulation.m' to create a unit geodesic sphere created by subdividing a regular icosahedron with normalised vertices. Source: https://it.mathworks.com/matlabcentral/fileexchange/50105-icosphere/files/icosphere.m)
