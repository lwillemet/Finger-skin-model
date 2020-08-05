# Finger-skin-model
2D parsimonious finite-difference model including frictional effects

# finger-skin-model-static.m
This script solves stresses and strains on the finger surface for an initial pressure P and a friction coefficient mu.
At the end of the execution, you will obtain two matrices U and F corresponding respectively to the displacement and the force.
The length of the matrices is 2*(Ndx+1) where Ndx is the number of particles (default value : 501). The number of columns of matrices is equal to the length of the time vector.

# RK4_f_U.m
This function is called in the script finger-skin-model-static to compute the displacement of each particles knowing the previous displacement of these elements and the forces exerted on each.
To solve this, we use Runge-Kutta solver at the 4th order.
