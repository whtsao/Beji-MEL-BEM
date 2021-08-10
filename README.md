# potential-flow-wave-problem

Developed by Wen-Huai Tsao at LSU-CE
Aug. 2021

Solve 2D wave problem by Eulerian - Lagrangian method

Eulerian step solves mixed type boundary value problem by boundary element method
- Linear element used for boundary integral
- Kinematic and dynamic boundary conditions for free surface
- Impermeable boundary condition for wave maker and bottom
- Radiation condition for absorption side

Lagrangian step tracks the free surface particle by Taylor series expansion
- 2nd order expansion
- Update velocity potential and location of a free surface particle

Output variables:
- Pressure on boundary node and domain node
- Velocity of boundary node and domain node

Next version to be continued:
- Regularized boundary integral method
- Output variables force and energy
- Floating body coupling
