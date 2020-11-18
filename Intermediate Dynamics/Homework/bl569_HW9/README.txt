%%%%%README%%%%%
%%%%%%%%%%%%%%%%

This homework covers number 30, 32, 33, 26.

num30_carBreaking
FRONT BRAKES CASE:
	Figure 1 is a plot of position in space for the case when the rear brakes are applied, the front is modeled as a skate. The star represents the initial and final position of the car's front wheel, CM, and rear wheel. Can see that it's unstable because it does a loopdeloop as it slows down.
	Figure 2 is a plot of velocity vs t, and omega vs t for the case when the rear brakes are applied, the front is modeled as a skate. This is of v of the front wheels. 
	Figure 3 is an animation of the front brakes case.
	Figure 4 is a plot of distance between the front and rear wheels as a function of time. This was used to just check the animation.

REAR BRAKES CASE:
	Figure 5 is a plot of position in space for the case when the front brakes are applied, the rear is modeled as a skate. The star represents the initial and final position of the car's front wheel, CM, and rear wheel. Stable because it basically travels in a straight line as it slows down.
	Figure 6 is a plot of velocity vs t, and omega vs t for the case when the front brakes are applied, the front is modeled as a skate. This is of v of the back wheels. 	
	Figure 7 is an animation of the front brakes case.

Mu = 0:
	Figure 8 is a plot of the rear brakes case without friction (muu =0). It goes further than the friction case. Can see that it's unstable because it does a loopdeloop as it slows down.
	Figure 9 is a plot of the front brakes case without friction (muu = 0). It goes further than the friction case.
Stable because it basically travels in a straight line as it slows down.	

	Parameters on 18-29 can be changed.
	To see the animation of the position of the rear brakes case, uncomment 102-118.
	To see the animation of the position of the front brakes case, uncomment lines 165-175. 
	To see special case of mu=0 ;set line 22=0.




num32_massSlotTurntable
	Figure 1 is a plot of theta vs time for multiple initial parameters/conditions. Can see that it doesn't approach


	

num26_montgomery8_pt2
	This uses fmincon to work some magic to making it periodic. I know the solution should be around what question 15 of the hw should get, ie. Montgomery's eight part 1. It doesn't quite get there right now because fmincon does not terminate when the cost function equals 0. 

	Figure 1 shows the system prior to fmincon doing things. 
	Figure 2 shows the system after fmincon, 1 period of it. 
	Figure 3 shows the system after fmincon, 5 periods of it. 
