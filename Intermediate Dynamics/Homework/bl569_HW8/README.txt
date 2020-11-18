%%%%%README%%%%%
%%%%%%%%%%%%%%%%

This homework covers number 28 and 29.
Number 28 is the double pendulum and can be found in the brianlui_doublePendulum folder. 
Number 29 is the Bead on parabolic wire and can be found in the brianlui_beadOnParabolicWIre folder. 

The code should be pretty well commented within the file. (Sorry my file isn't called root)

Process: derive the EOM using a derive file, then get the EOM into a file to output the animation. 


Double Pendulum
num28_doublePendulum
	This is THE file to execute. It produces the simulation. 
	There are three derivations within it: AMB, Lagrange, and DAE. 
	
	FIGURES PRODUCED:
	Figure 1. The simulation of the double pendulum derived using AMB.
	Figure 2. Just the path of the end of the double pendulum for all 30 secs, AMB. 
	Figure 3. The simulation of the double pendulum derived using DAE.  
	Figure 4. Just the path of the end of the double pendulum for all 30 secs, DAE. 
	Figure 5. The simulation of the double pendulum derived using Lagrange.
	Figure 6. Just the path of the end of the double pendulum for all 30 secs, Lagrange. 
 

To derive the stuff to put in the rhs equation, ie. the EOM, 
dPendDeriveEOMamb, dPendDeriveEOMLagrange, and dPendDeriveEODAE are used. They are for each of the three methods. 

These dPendDeriveEOM files respectively produce a rhsStuffmassMatrixY and rhsStuffbVectorY file, in which Y is AMB, Lagrange, or DAE. These two files get called in the rhs function within the num28_doublePendulum file.

Interesting functions called within the dPendDeriveEOM files: 
	A = jacobian(f, thdots).
		If f = [f1 f2 f3], thdots = [th1, th2, th3]. 		Produces a matrix A containing the partial derivatives of f with respect to thdots; ie. the first row of A: A(1,1) = df1/dth1, A(1,2) = df1/dth2, A(1,3) = df1/dth3. 
	matlabFunction(m, 'file','nameOfFile')
		Produces a file named 'nameOfFile that outputs the value of m.
	simplify(m)- simplifies expression m. Real useful.


Bead on Parabolic Wire
num29_beadOnParabolicWire
	This is the file to execute. It produces the motions and simulation of the situation. 

	FIGURES PRODUCED: 
	Figure 1. Path of bead, NE method
	Figure 2. Path of bead, DAE method (currently doesn't work :( )
	Figure 3. Path of bead, Lagrange
	Figure 4. Errors between the methods
	Figure 5. Nothing. Simulation keeps rewritting my earlier figures so this is meant to be blank.
	FIgure 6. Animation 

Other files in this folder: 
	beadDeriveEOM_x, where x=DAE,NE, or Lagrange. These files derive the EOM. 
	rhs_x, where x=DAE, NE, or Lagrange. THese are the derived EOM; output of beadDeriveEOM_x. 

So at the moment, I think my NE and Lagrange are working. 
Bead on wire by DAE doesn't seem to be. I think it's on how I'm defining my et/en vectors. (tangent and normal vectors)
	Potential error1: en always has to be pointing inward, so it's not always the cross of e3 and et. I put a sign(xdot) in there that I think should fix this. 
	Potential error 2: et shouldn't be defined with xdot and ydot terms. But shouldn't et always be pointing in the direction of the velocity vector?