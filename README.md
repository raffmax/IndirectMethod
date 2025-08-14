# IndirectMethod

This repository contains MATLAB code for simulating the compass-gait walker and robot RABBIT, based on the paper:

> **"The Indirect Method for Generating Optimal Periodic Trajectories and Its Application to Economical Bipedal Walking."**

The code demonstrates how to generate optimal gaits from passive gaits using both the **indirect method** and the **direct method** described in the paper.

## Contents

Both the **indirect** and **direct** approaches are implemented for:
- A 2-DOF compass-gait walker  
- A 5-DOF biped robot (*RABBIT*)

### Main Scripts
- **`indirectPassive2levelGround.m`** – Generates a 1D family of passive gaits using the indirect method.  
- **`directPassive2levelGround.m`** – Generates a 1D family of passive gaits using the direct method.

## Dependencies
- **MATLAB** (tested on version 2023a)  
- **CasADi** – Required for RABBIT’s dynamics computation.  
  Download here: [https://github.com/casadi/casadi](https://github.com/casadi/casadi)