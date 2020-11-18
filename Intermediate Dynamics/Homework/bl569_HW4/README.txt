%%%%%README%%%%%
%%%%%%%%%%%%%%%%

This homework covers number 7,12-14. Also 10b.
All the executable files in the folder are scripts.
The code should be pretty well commented within the file.

polarDynamics_num7
    Lines 16-30 can be changed. Should be pretty self explanatory as to what they represent.
    Each figure produced represent the plots of lines using the initial conditions in lines 16-24.
    The difference between each figure is the tolerance of ode45, with figure i corresponding to tolerances(i).



montgomery8_num15
    makes some pretty graphs! Makes the number of plots as length(durArray).
    Figure 1 = part a- run for 2.1 time units
    Figure 2 = part b- run for 10 time units
    Figure 3 = part c- run for 10 time units and IC changed a little
    Figure 4 = part d- run for 100 time units and IC changed a little
    
    For figure 3 and 4, i changed x1 and y1 by -0.1
    
    Can play with IC (lines 10-22), parameters (line 27), duration (line 30) how much IC changes (lines 35,36)



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
	
	

		
