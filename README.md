# Finger-skin-model
This code provides a 2D parsimonious finite-difference model of human skin capturing the deformation of the skin, viscoelastic effects through the damper, and local friction behavior at the interface.
The model is constructed with a bottom-up approach, using as few parameters to explain a wide array of observed phenomena. It is composed of a chain of massless elements maintained together by elastic springs. This chain can be assimilated to the external layer of the skin (the epidermis). Its shape is maintained using other elastic springs that connect the massless elements to a virtual bone, analogous to the mechanical behavior of the subcutaneous tissues. The two elements on the outside of the membrane are also attached to the bone and model the effect of the rigid nail. Overall, the model resemble the discrete version of a curved elastic membrane on a spring foundation. Viscosity of the skin is modeled by dampers, connecting each particles to the mass of the system. 

![model_connectivity](https://user-images.githubusercontent.com/58175648/115708505-61f45800-a370-11eb-9294-5a7a5baaf0d0.png)

External forces are applied on the bone element and the frictional forces are computed with the Dahl's friction model. The displacements were solved using Runge-Kutta at order 4. For the equations, please refer to the submission.

# finger-skin-model-static.m
This script solves stresses and strains on the finger surface for an initial pressure P and a friction coefficient mu.
At the end of the execution, you will obtain two matrices U and F corresponding respectively to the displacements and the forces acting on each element. The first element corresponds to the bone and the other from the left one to right one.
The length of the matrices is 2*(Ndx+1) where Ndx is the number of particles (default value : 501). The first half of the vector is the normal displacement (resp. force) and the second half, the tangential displacement (resp. force). The number of columns of matrices is equal to the length of the time vector.

# curvspace.m
This function is called at the beginning of the script finger-skin-model-static and creates an evenly spaced points distribution along the external layer of the finger in 2D. Then, at the beginning of the stimulation each spring has the same length at rest.

# RK4_f_U.m
This function is called in the script finger-skin-model-static to compute the displacement of each particles knowing the previous displacement of these elements and the forces exerted on each. As Runge-Kutta is used at the 4th order, the function is called 4 times.
The matrices of springs and dampers dependency are built in this function as following:

![matrices_dependencies](https://user-images.githubusercontent.com/58175648/115710399-99fc9a80-a372-11eb-9e20-84fb4c629270.png)


