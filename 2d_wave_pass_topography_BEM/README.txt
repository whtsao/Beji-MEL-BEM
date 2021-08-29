************************************************************************************************
Solve 2D wave problem by Eulerian - Lagrangian method
Developed by Wen-Huai Tsao at LSU-CE Aug. 2021

************************************************************************************************
Description:

1. Domain dimension = [2]. Time interval = 0.001s.
2. Four kinds of boundaries: Free surface (y+), Left side wave maker (x-), Right side absorption (x+), Bottom (y-).
3. The wave maker is a vertical impermeable baffle. Amplitude of the piston wavemaker = [3]
4. Generate periodic waves by x = A*sin(omega*t), A=amplitude(m), omega-forcing frequency(rad/s)
5. Generate solitary wave by half-sine displacement-time history with one stroke(m) in one duration(s).
6. The bottom can have any kind of topograghy by giving nodes between each section
************************************************************************************************
Method [1]:

A. Eulerian step solves mixed type boundary value problem by boundary element method
  1. Linear 2-node element is used for boundary integral, 8 Gaussian quadrature points are used in one element.
  2. Kinematic and dynamic boundary conditions for free surface
  3. Impermeable boundary condition for wave maker and bottom
  4. Radiation condition for absorption side (No reflection). Wave speed C is calculated via dispersion relation before program start
  5. An iterative solver (due to radiation BC) for the velocity potentials on the absorption side is embedded

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
    amp=amplitude (if 1)
    amp=stroke (if 2)
    omega=omega (if 1)
    omega=duration (if 2)
    ignore psi if using solitary wave
3. io.dat = external excitations
    (time, displacement, velocity, acceleration)
4. s.dat = boundary node
    [x coordinate; y coordinate] at ith step
    [x coordinate; y coordinate] at i+1th step
    ...
6. domain.dat = scatter data of all nodes in the entire domain
    [x coordinate; y coordinate; u, v, p] at ith step
    [x coordinate; y coordinate; u, v, p] at i+1th step
    ...
7. wg.dat = wave elevations of the wave gauges set along x-direction.
   (time, elevation 1, elevation 2,... elevation n)
8. err.dat = converged error of velocity potentials of the nodes on the radiation outlet
    (time, iteration count, ABS error, RMS error)

************************************************************************************************
Reference:
[1] S. T. Griili, J. Skourup and I. A. Svendsen. An efficient boundary element method for nonlinear water waves. Engineering Analysis with Boundary Elements, 1989, 6(2), 97-107.
[2] S. Beji, J.A. Battjes. Numerical simulation of nonlinear wave propagation over a bar. Coastal Engineering, 1994, 23, 1-16.
[3] F. Ursell, R. G. Dean, Y. S. Y u. Forced small-amplitude water waves: a comparison of theory and experiment. Journal of Fluid Mechanics, 1960, 7(1) , 33-52.

************************************************************************************************
NO WARRANTY:
There is no warranty, expressed or implied, for the use of this package.
The authors are not responsible for any possible damages in using the software.
