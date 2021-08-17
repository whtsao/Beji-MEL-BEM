************************************************************************************************
Solve 2D wave problem by Eulerian - Lagrangian method
Developed by Wen-Huai Tsao at LSU-CE Aug. 2021

************************************************************************************************
Description:

1. Domain dimension = 50m x 8m. Time interval = 0.001s.
2. Four boundaries: Free surface (y+), Left side wave maker (x-), Right side absorption (x+), Bottom (y-).
3. The wave maker is a vertical impermeable baffle.
4. Generate periodic waves by x = A*sin(omega*t), A=amplitude(m), omega-forcing frequency(rad/s)
5. Generate solitary wave by half-sine displacement-time history with one stroke(m) in one duration(s).
6. The bottom can have any kind of topograghy by giving nodes between each section
************************************************************************************************
Method:

A. Eulerian step solves mixed type boundary value problem by boundary element method
  1. Linear 2-node element is used for boundary integral, 8 Gaussian quadrature points are used in one element.
  2. Kinematic and dynamic boundary conditions for free surface
  3. Impermeable boundary condition for wave maker and bottom
  4. Radiation condition for absorption side (No reflection)
  5. An iterative solver for the radiation BC on the absorption side is embedded

B. Lagrangian step tracks the free surface particle by 2nd-order Taylor series expansion
  1. Update velocity potential and location of a free surface particle

C. Output variables:
  1. Pressure on boundary node and domain node
  2. Velocity of boundary node and domain node

************************************************************************************************
Usage:
1. In 1.ipt
  Enter coordinate of corner nodes from left-bottom node, then in counterclockwise
  Enter number of element on each planes from free surface, then in clockwise
  Enter boundary condition type of each planes from free surface, then in clockwise
    1=phi is given (usually used for free surface)
    2=pphi is given (usually used for absoption and wall)
2.ipt = input file
  Enter number of plane, e.g. 4 for simple flat-bottom channel
  Enter wave generation type
    1=periodic waves
    2=solitary wave
    Note: the variable in code amp=amplitude=stroke, omega=omega=duration, ignore psi if using solitary wave
3. io.dat = external excitation displacement, velocity, and acceleration
4. s.dat = boundary node
5. w.dat = wave elevation on the wetted vertical boundaries
6. domain.dat = location, velocity, and pressure of all nodes in the entire domain
7. err.dat = converge error of the iterative solver

************************************************************************************************
Reference:
S. T. Griili, J. Skourup and I. A. Svendsen. An efficient boundary element method for nonlinear water waves. Engineering Analysis with Boundary Elements, 1989, Vol. 6, No. 2, 97-107.

************************************************************************************************
NO WARRANTY:
There is no warranty, expressed or implied, for the use of this package.
The authors are not responsible for any possible damages in using the software.