%%%%%README%%%%%
%%%%%%%%%%%%%%%%

This homework covers number 10b,12-14.
All the executable files in the folder are scripts.
The code should be pretty well commented within the file.


centralForce_num10
	I DID IT! At least good enough I think. 
	I tried using the fmincon function, setting constraints on it.   	

	parameters on lines 10-19 can be toggled. 
	to change the central force law, change it within rhs equation
	The commented out code in lines 21-41 are showing what it looks like before fmincon is applied. This is figure 1.
	The plot after fmincon is applied is figure 

	Process: 
	After setting up fmincon, I played around with initial conditions. Initially, I had fmincon minimizing the difference between r(1)-r(end); theta(1)- theta(end); rDot(1)- rDot(end); thetaDot(1)-  thetaDot(end). Playing around with the initial conditions was giving me solutions that was linear-like, but didn't meet. 
	From this, I realized that the constraints on the initial conditions should be that at least rDot > 0 and thetaDot >0. If thetaDot>0 initially, then it can't be a straight line. If rDot >0 along with thetaDot>0, then it can't be a circle. These conditions were set in line 51 (lb).
	Afterward, I switched fmincon to minimize x(1)- x(end); y(1)- y(end); xDot(1)-xDot(end); yDot(1)-yDot(end). This got me a much more awesome graph! After changing some initial conditions, I got something pretty good (although not perfect I think but good enough). 

	 
	Figure 1 is the initial one before fmincon.
	Figure 2 is the one after fmincon! 
	Very satisfying.
	

	
mechanicsParticles_num14b
	The parameters on lines 17-20 can be changed.
	The number of steps for part b can also be toggled. It is on line 27. 
	The equations of motion were done in polar.
	The plot is set to go for 1 period (derived theoretically), if it makes 2 circles, then thats the actual period.
	
	Figure 1 is the complete trajectory of the two masses.
	Figure 2 is the animation of the two masses as it moves.

mechanicsParticles_num14b
	The parameters on lines 15-19 can be changed, though G, m1,m2,m3,d are all given by the problem.
	Line 41 can also be modified: you can make it go the other direction.
	Lines 85-91 contain commented out code- this does the animation. 
	The output is a vector containing the error. Because each mass only moves for a third of the period, 
		the error terms are taken as the final position of m1 - starting position of m2, etc. 
	Figure 1 shows the plot of the three masses.
		
