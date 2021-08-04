# potential-flow-wave-problem

Developed by Wen-Huai Tsao at LSU
Aug. 2021

Solve wave problem by Eulerian - Lagrangian method

Eulerian step solves mixed type boundary value problem by boundary element method
- Linear element used for boundary integral
- Kinematic and dynamic boundary condition for free surface
- Impermeable boundary condition for wave maker and bottom
- Radiation condition for absorption side

Lagrangian step tracks the free surface particle by Taylor series expansion
- 2nd order expansion
- Update velocity potential and location of a free surface particle

Output variables:
- Pressure
- Velocity

Next version to be continued:
- 3D case
- Regularized boundary integral method
- Output variables force and energy
- Floating body coupling
