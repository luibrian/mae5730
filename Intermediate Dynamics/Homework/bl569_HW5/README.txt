%%%%%README%%%%%
%%%%%%%%%%%%%%%%

This homework covers number 6,18-22.

All the executable files in the folder are scripts.
The code should be pretty well commented within the file.

twoMasses_num18
    Lines 16-30 can be changed. Should be pretty self explanatory as to what they represent.
    Line 66 can also be changed. This is the ks used in part d.

    Figure 1 is the position vs time graph of x1, x2, and xG (position of CM) for k=1
    Figure 2 is the animation of x1 and x2
    Figure 3-7 is the error vs time (x2-xG vs t) and position vs time (x1,x2,xG vs t)
        for ks in line 66. Each succeeding figure is of an increasing k value. 
        The k values I have shown are kArray = [1 10 100 1000 10000];

    The magnitude of the error decrease as k increases- you can see this in
        the scale of the y axis in figures 3-7. 
    The errors due to numerical inaccuracies increase over time for higher k;
        this is more visible for higher k because of the smaller magnitude of 
        the difference between x2-xG. 
    Looking at the second plot within figures 3-7, you can see that the lines
        get increasingly more straight. 


twoMassesConstrained_num19
    Lines 12- 26 are changeable- pretty self explanatory what they refer to. 
    Line 64-65 also are changeable- this is the initial speed difference between mass 1 and mass 2. 

    Figure 1: x1, x2, xG vs time
    Figure 2: x1-xG vs time, x2-x1 vs time
    Figure 3: x1 and x2 for different initial velocities.

    This is much more accurate than the previous question. The distance from a mass to the CM is basically constant, and the distance between the masses between each other is basically satisfied for all time (some miniscule variation in these values but we talking 10^-13 ish). This can be seen in figures 1 and 2. 
    Looking at figure 3, with unequal IC, it doesn't satisfy the constraint (x2Dot =/= x1Dot). As a result, mass 2 moves further and further away from mass 1. 


pendulum_num21
    Parameters on lines 19-28 are changeable, as are 55-57. 
    
    Figure 1 is the motion of the pendulum solved through DAE.
    FIgure 2 is the motion of the pendulum solved from the solution of the simple pendulum equations (thetaDoubleDot = -g/el sin theta.
    Figure 3 is the difference in the xposition of the pendulum.
    Figure 4 is the animation of figure 1- have to uncomment 44-50 first.

    As t increases, the drift away from the satisfying the kinematic constraint increases, as can be seen in figure 3. 
    Solving linear equations (DAE method) leads to numerical inaccuracies that add up over time. 


pendulumAwkward_num22
    Parameters on lines 19-28 are changeable, as are 53-55. 

    Figure 1 is the motion of the pendulum solved through DAE.
    FIgure 2 is the motion of the pendulum solved from the solution of the simple pendulum equations (thetaDoubleDot = -g/el sin theta.
    Figure 3 is the difference in the yposition of the pendulum.
    Figure 4 is the animation of figure 1- have to uncomment 44-50 first.
    
    As t increases, the drift away from the satisfying the kinematic constraint increases, as can be seen in figure 3. 
    Solving linear equations (DAE method) leads to numerical inaccuracies that add up over time. 
    The scale of the error in this part is less than the previous part. 

