%%%%%README%%%%%
%%%%%%%%%%%%%%%%

This homework covers number 8-11.
All the executable files in the folder are scripts.
The code should be pretty well commented within the file.

cannonBall_num8

	The parameters on lines 10-24 affect parts b

	The outputs of the script are 4 figures and two outputs in the command window. 
		figure 1: numerical solution using the parameters in lines 10-24 (given by problem statement)
		figure 2: plot of the trajectory using varying tolerances and using the same parameter values as part b
		figure 3: plot of the trajectory with varying starting velocities
		figure 4: plot of the trajectory with v0 = infinity (not actually but some large v0)
	
		errors is a matrix of the errors at time t=2 for relative and absolute errors of the tolerances on line 54.
			tolerances = [1e-10 1e-8 1e-6 1e-4 1e-2 1];
		thetaStar is the theta for max range at the initial velocities given on line 138.	

	
	
massOnSpring3D_num9
	parameters on lines 11-17 can be toggled
	add initial conditions by adding to matrix in lines 21-26
	each figure i outputted represents the initial conditions given by row i of ic 
	Note that the images are stretched! Look at the actual scale on the plots
	
	
	
centralForce_num10
	parameters on lines 10-24 can be toggled. 
	to change the central force law, change it within rhs equation
	to see some of my attempts, check within the rhs function
	
	i tried a whole bunch of central force forces and got a bunch of pretty plots, but couldn't find anything periodic	
	
	figure 1: nice looking non periodic plot
	
	
	
cannonBall_num11
	The parameters on lines 10-23 can be changed. They're on 
	The number of steps for part b can also be toggled. It is on line 48. 
	
	figure 1 is the plot of the cannon ball using the parameters in 10-23. 
	figure 2 is the difference between KE and drag work as the step sizes increase
	