This project is for PHYS 251 (Computational Physics). It approximates the solutions to ordinary differential equations (ODEs) governing the gravitational interactions between bodies. 
The Sun is assumed to be the reference frame due to its significantly larger mass compared to other bodies, and it is treated as stationary in this frame.

The gravitational force between all bodies is calculated, and the system's dynamics are solved using three different numerical methods:

1. Euler-Cromer Method (ec): A simple, explicit method for time integration, updating both position and velocity.
2. Velocity Verlet Method (vv): An improvement over Euler's method, providing more accurate results by updating positions and velocities in a symplectic manner.
3. SciPy `solve_ivp` Method (solve_ivp): Uses high-order solvers to numerically solve the ODEs with built-in methods from the SciPy library.

These methods are used to approximate the motion of bodies in the system, and the accuracy and performance of each method are compared.
