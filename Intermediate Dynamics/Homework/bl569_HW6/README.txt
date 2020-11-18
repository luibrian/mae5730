%%%%%README%%%%%
%%%%%%%%%%%%%%%%

This homework covers number 29,23,24,27.

All the executable files in the folder are scripts.
The code should be pretty well commented within the file.


dumbbell_num23
	Lines 14-20 can be changed. As can lines 25-27. They refer to the parameters, timescale, and initial conditions of the problem. 
	Figure 1 is the animation of the two masses moving under the EOM. 
	Figure 2 is the error in the constraint equation over time. We see that it increases as numerical errors accumulate. 
	Figure 3 is the error of the conservation of energy over time. We see that it increases as numerical errors accumulate. 
	Figure 4 is the error in the angular momentum conservation. We see that it increases as numerical errors accumulate. 
	Figure 5 is the error in the linear momentum conservation. We see that it increases as numerical errors accumulate. 

	The overall error of these are pretty small- seems to suggest that the EOM is good. 
	Increasing the tolerance on ode45, we see that the magnitude of the errors produced decrease. You can change the tolerance on ode45 on line 36. This also suggests that the errors plotted are numerical errors
	The magnitude of the error matters; the sign of the magnitude doesn't matter. This is because the values of these checks should equal 0 and the fact that it isn't is the error. In other words, doing the energy at time t - energy at time 0 is equal in terms of the error as energy at time 0 - energy at time t.  



shakingPendulum_num27
	Lines 12-28 can be changed. These are the parameters, timescale, and IC. 
	Figure 1 is the plot of the end of the pendulum in the x-y plane, as well as how theta changes with time, done with the minimal coordinates approach.
	Figure 2 (you need to uncomment it out in lines 50-63) is an animation of the end of the pendulum in the x-y plane. 
	Figure 3 is the plot of the end of the pendulum in the x-y plane, as well as how theta changes with time, done with the DAE approach. 
	Figure 4 is the difference in theta between the two approaches. 

The pendulum oscillates back and forth forever because there is no damping on it. If there was, it would tend to the peak value at (0,10), as seen by the demonstration in class. It isn't unstable because the motion of the theta is oscillates between the +/- the initial theta supplied. 