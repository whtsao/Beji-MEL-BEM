************************************************************************************************
Solve 3D wave problem by Eulerian - Lagrangian method
Developed by Wen-Huai Tsao at LSU-CE Aug. 2021

************************************************************************************************
Description:

1. Domain dimension = 50m x 5m x 8m. Time interval = 0.001s.
2. Six boundaries: Free surface (z+), Wave maker (x-), Side absorption (x+, y+, y-), Flat bottom (z-).
3. The wave maker is a rotating baffle (period T=3s). sita = A*sin(omega*t)

************************************************************************************************
Method:

A. Eulerian step solves mixed type boundary value problem by boundary element method
  1. Linear 4-node element is used for boundary integral, 7 Gaussian quadrature points are used in one element.
  2. Kinematic and dynamic boundary conditions for free surface
  3. Impermeable boundary condition for wave maker and bottom
  4. Radiation condition for absorption sides (No reflection)
  5. An iterative solver for the radiation BC on the absorption sides is embedded

B. Lagrangian step tracks the free surface particle by 2nd-order Taylor series expansion
  1. Update velocity potential and location of a free surface particle

C. Output variables:
  1. Pressure on boundary node and domain node
  2. Velocity of boundary node and domain node

************************************************************************************************
Usage:
1. 1.ipt = input file
2. io.dat = external excitation displacement, velocity, and acceleration
3. s.dat = boundary node
4. w.dat = wave elevation on the wetted vertical boundaries
5. domain.dat = location, velocity, and pressure of all nodes in the entire domain
6. err.dat = converge error of the iterative solver

************************************************************************************************
Reference:
S. T. Griili, J. Skourup and I. A. Svendsen. An efficient boundary element method for nonlinear water waves. Engineering Analysis with Boundary Elements, 1989, Vol. 6, No. 2, 97-107.

************************************************************************************************
NO WARRANTY:
There is no warranty, expressed or implied, for the use of this package.
The authors are not responsible for any possible damages in using the software.



