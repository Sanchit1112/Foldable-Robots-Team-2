# System Dynamics
## Team 2
## Kevin Julius, Romney Kellogg, Sanchit Singhal, Siddhaarthan Akila Dhakshinamoorthy


## 1. Robot Figures:

### 1.1 Dimensioned Figure:
Below is a dimensioned figure of the solidworks model of the robot.
![solidworkdraw1](attachment:solidworkdraw1)

![solidworkdraw2](attachment:solidworkdraw2)

### 1.2 Dynamics Figure:
Below are 2 figures. The first of which is a detailed diagram of the dynamics of the system, the second of which is a simplified diagram showing what aspects of the dynamics were modeled. 

![DynamicsIIFig1.PNG](attachment:DynamicsIIFig1.PNG)


![DynamicsIIFig2Simple.PNG](attachment:DynamicsIIFig2Simple.PNG)


```python
%matplotlib inline
```


```python
#Importing Libraries
import pynamics
from pynamics.frame import Frame
from pynamics.variable_types import Differentiable,Constant,Variable
from pynamics.system import System
from pynamics.body import Body
from pynamics.dyadic import Dyadic
from pynamics.output import Output,PointsOutput
from pynamics.particle import Particle
import pynamics.integration
import pynamics.tanh

import sympy
import numpy
import matplotlib.pyplot as plt
plt.ion()
from math import pi

from math import degrees, radians
from pynamics.constraint import Constraint
import scipy.optimize
```

## 2. Constants

### 2.1 Link Lengths:

All link lengths in the saurus linkage are the same length of 60mm or 0.06m. lDi1 (.02m) and lDi2 (.04m) are lengths that makeup link D, that has been given a virtual joint that will be discussed below.


```python
# Initializing Pynamics
system = System()
pynamics.set_system(__name__,system)
# Defining Link Constants
lAi=.060 #all in m
lBi=.060
lCi=.060
lDi1=.015
lDi2=.045
lEi=.060
lFi=.060

lA = Constant(lAi,'lA',system)
lB = Constant(lBi,'lB',system)
lC = Constant(lCi,'lC',system)
lD1 = Constant(lDi1,'lD1',system)
lD2 = Constant(lDi2,'lD2',system)
lE = Constant(lEi,'lE',system)
lF = Constant(lFi,'lF',system)

```

### 2.2 Rigid Body Lengths:

The below lengths makeup the two endcaps on the left and right side of the robot housing mechanical and electrical components. Left Endcap = G>H>I, Right Endcap =  J>K>L


```python
lGi=.060 #all in m
lHi=.060
lIi=.060
lJi=.060
lKi=.060
lLi=.060

lG = Constant(lGi,'lG',system)
lH = Constant(lHi,'lH',system)
lI = Constant(lIi,'lI',system)
lJ = Constant(lJi,'lJ',system)
lK = Constant(lKi,'lK',system)
lL = Constant(lLi,'lL',system)
```

### 2.3 Link Masses:

Below are the dimensions of each link of .06mx.06mx.004m. These measurements were used to calculate volume and multiplied by the density of the cardstock to determine the mass of each link. Once again link D has been split into two masses that will be discussed further below.


```python
#mass calculated using density of cardstock and volume
a=.06 #m
b2=.06 #m
c=.004 #m
rho=689 #kg/m^3
m=a*b2*c*rho


#links
mA = Constant(m,'mA',system) #in kg
mB = Constant(m,'mB',system)
mC = Constant(m,'mC',system)
mD1 = Constant(m/4,'mD1',system)
mD2 = Constant(3*m/4,'mD2',system)
mE = Constant(m,'mE',system)
mF = Constant(m,'mF',system)
```

### 2.4 Component and Endcap Lumped Masses:

Below are the masses of the endcaps and all components within each side summed into a lump mass on each side. This simplification can be performed because the components on each rigid end of the robot should not experience rotation, therefor their moment of inertia will be irrelevant to the motion of the robot.

Table of Masses:

| Item | Mass(kg) |
| --- | --- | 
| Overall System | .1101 |
| Endcap x2 | .00121 |
| Saurus Links x6 | .0003 |
| Power Supply + Microcontroller | .045 | 
| Motor 1 | .009 |
| Motor 2 | .009 |
| Stopper 1 | .0033 |
| Stopper 2 | .0041 |


```python
#other parts
endcap=.00121
psmc=.045
motor=.009
stop1=.0033
stop2=.0041
mlefts=endcap+stop1
mrights=endcap+motor+stop2+motor+psmc
mleft = Constant(mlefts,'mleft',system)
mright = Constant(mrights,'mright',system)
```

### 2.5 Joint Stiffness:

These values were calculated using motion tracking of a prototype joint that will be used in the construction fo the system. The method that these values were calculated can be seen further in the damping.pdf document

It is important to include joint stiffness because each joint has an inherent elasticity to it that needs to be modeled.

Unfortantely data was collected, but was the joint stiffness value was unable to be calculated from the data after attempting different methods.


```python
k = Constant(.01,'k',system)
```

### 2.6 Joint Damping:

These values were calculated using motion tracking of a prototype joint that will be used in the construction fo the system. The method that these values were calculated can be seen further in the damping.pdf document

It is important to include this because each joint has an inherent damping to it that needs to be modeled to calculate energy lost through the joints in the system.

Unfortantely data was collected, but was the damping value was unable to be calculated from the data after attempting different methods.


```python
b = Constant(0.001,'b',system)
```

### 2.7 Adding Compliance:

Compliance was added to link D because it is to be actuated by a rigid link connected to the motor. This will cause link D to likely experience more force than the other links and a good candidate to study the effect of compliance on the system.

One can see further calculations within the ComplianceRK.ipynb document.

The damping constant was found to be 0/negligable
The stiffness/spring constant was found to be 6.05656 N/m

The location that the link was to be placed at was found to be .25 along the length of the beam which was represented in the constants defined above. lD=.06 lD1=lD*.25=.015 lD2=lD-lD2=.045


```python
bD = Constant(0,'bD',system)
kD = Constant(6.05656,'kD',system)
```

### 2.8 Preloading:

The system is to be fold to be relaxed at the current initial condition so all preload values were set to zero.

If an alternative initial condition was used for the rotational states, these preloads would be altered to represent the change in those states from the current initial position states.


```python
preload1 = Constant(0*pi/180,'preload1',system)
preload2 = Constant(0*pi/180,'preload2',system)
preload3 = Constant(0*pi/180,'preload3',system)
preload41 = Constant(0*pi/180,'preload41',system)
preload42 = Constant(0*pi/180,'preload42',system)
preload5 = Constant(0*pi/180,'preload5',system)
preload6 = Constant(0*pi/180,'preload6',system)
```

## 3. Defining Frames and State/Differential Variables

### 3.1 Defining Time Frame, Time Step, and Integration Tolerance:

The time frame is set to end once the system has maintained a Steady State. Time step is set to be 1 step per frame at 30 frames a second (1/30).
Integration Tolerance is set to 1e-4 to keep some integration accuracy, while allowing a quicker integration time. 


```python
tinitial = 0
tfinal = 5
fps = 30
tstep = 1/fps
t = numpy.r_[tinitial:tfinal:tstep]
tol = 1e-4
```

### 3.2 Establishing State Variables and Their Derivatives:

Below a state variable and their derivatives are created for each link/frame. Link D has an added virtual joint for compliance requiring two differential variables.

Additional variables x1 and y1 (and their derivatives) are created to break free from the Newtonian Frame and create a moving system. These measure the distance of the joint of links A and F (pNA) from the origin in the Nx and Ny direction.


```python
# Defining State Variables and their derivatives
x1,x1_d,x1_dd = Differentiable('x1',system)
y1,y1_d,y1_dd = Differentiable('y1',system)
qA,qA_d,qA_dd = Differentiable('qA',system)
qB,qB_d,qB_dd = Differentiable('qB',system)
qC,qC_d,qC_dd = Differentiable('qC',system)
qD1,qD1_d,qD1_dd = Differentiable('qD1',system)
qD2,qD2_d,qD2_dd = Differentiable('qD2',system)
qE,qE_d,qE_dd = Differentiable('qE',system)
qF,qF_d,qF_dd = Differentiable('qF',system)
```

### 3.3 Creating Frames and Their Relationships:
Below the Newtonian frame is declared along with one frame for each link within the saurus linkage A-F. Liink D once again has two frames due to the virtual joint within link D

Afterwards the Newtonian frame is set.

Sequentially after the Newtonian frame the frames are defined by a z rotation equal to each state variable.
A is rotated qA from N

B is rotated qB from A

C is rotated qC from B
"
D1 is rotated qD1 from C

D2 is rotated qD2 from D1

E is rotated qE from D2

F is rotated qF from E


```python
# Declaring Frames
N = Frame('N')
A = Frame('A')
B = Frame('B')
C = Frame('C')
D1 = Frame('D1')
D2 = Frame('D2')
E = Frame('E')
F = Frame('F')

```


```python
# Placing Newtonian Frame
system.set_newtonian(N)
```


```python
# Establishing Frame Rotation Relationships
A.rotate_fixed_axis_directed(N,[0,0,1],qA,system)
B.rotate_fixed_axis_directed(A,[0,0,1],qB,system)
C.rotate_fixed_axis_directed(B,[0,0,1],qC,system)
D1.rotate_fixed_axis_directed(C,[0,0,1],qD1,system)
D2.rotate_fixed_axis_directed(D1,[0,0,1],qD2,system)
E.rotate_fixed_axis_directed(D2,[0,0,1],qE,system)
F.rotate_fixed_axis_directed(E,[0,0,1],qF,system)
```

## 4. Defining Robot Kinematic Relationships and Masses

### 4.1 Defining Location of Saurus Linkage Points using Kinematic Relationships
Here the point/joint locations are defined.

The point NA is defined by a displacement x1 and y1 from the origin.

The points that follow are all defined by the previous point and the lenght and orientation of each link in the order:

pNA> (Link A) >pAB> (Link B) >pBC> (Link C) >pCD1> (Link D1) >pD1D2> (Link D2) >pD2E> (Link E) >pEF> (Link F) > pFtip = pNA


```python
# Defining Point Locations based on kinematics of the system
pN = 0*N.x+0*N.y+0*N.z
pNA = x1*N.x+y1*N.y
pAB = pNA + lA*A.x
pBC = pAB + lB*B.x
pCD1 = pBC + lC*C.x
pD1D2 = pCD1 + lD1*D1.x
pD2E = pD1D2 + lD2*D2.x
pEF = pD2E+ lE*E.x
pFtip= pEF + lF*F.x
points = [pNA,pAB,pBC,pCD1,pD1D2,pD2E,pEF,pFtip]
```

### 4.2 Defining Location of Links Forming Rigid Endcap:

Below are the points defining the corner of the two rigid end caps. These are defined by: 

Left: pEF> (Link G) >pGH> (Link H) >pHI> (Link I) >pNA

Right: pBC> (Link J) >pJK> (Link K) >pKL> (Link L) >pCDD1

These links do not have frames as they are rigidly defined by the F and C frames for the left and right endcap respectively.


```python
# Points used for graphical representation end caps holding components
pGH = pEF-lG*F.y
pHI = pGH+lH*F.x
pJK = pBC-lJ*C.y
pKL = pJK+lK*C.x
```

### 4.3 Defining Center of Mass of Saurus Linkage Links:

Here the center of masses are defined for each link to be at the center of each link.


```python
#Center of Masses of Links
pAcm=pNA+lA/2*A.x
pBcm=pAB+lB/2*B.x
pCcm=pBC+lC/2*C.x
pD1cm=pCD1+lD1/2*D1.x
pD2cm=pD1D2+lD2/2*D2.x
pEcm=pD2E+lE/2*E.x
pFcm=pEF+lF/2*F.x
```

### 4.4 Defining Center of Mass of Lumped Masses:

Below the center of masses of the components within and links composing the endcaps are placed at the center of the left and right endcap boxes. This is an approximation, but as long as the center of mass is within the boxes and these masses are't experiencing large amount of rotation this approximation should hold true.


```python
#Center of Masses Each Endcap Holding Components
pmright=(pNA+pEF)*.5-.03*F.y
pmleft= (pBC+pCD1)*.5-.03*F.y
```

### 4.5 Creating Particles Using Center of Masses and Masses:

Particles masses were used to model all bodies in the system to limit complexity and computation time. This idealization should have little effect on results due to the incredibly small moment of inertia of the links, and the lack of rotation of the large masses in each end cap (saurus linkage leads to linear actuation of both lumped masses).

Below how the moment of inertia would be calculate for each link is listed, but these rotational inertias were not implemented in the final model due to the use of particles because of their small size shown below. (Izz was shown because it would be the most impactful moment, since it is along the links only axis of rotation).


```python
#6 equal size links have same Ixx,Iyy,Izz
Ixx=(1/12)*m*(b2**2+c**2)
Iyy=(1/12)*m*(a**2+c**2)
Izz=(1/12)*m*(a**2+b2**2)
```


```python
Izz
```




    5.952959999999999e-06




```python
PA = Particle(pAcm,mA,'PA',system)
PB = Particle(pBcm,mB,'PB',system)
PC = Particle(pCcm,mC,'PC',system)
PD1 = Particle(pD1cm,mD1,'PD1',system)
PD2 = Particle(pD2cm,mD2,'PD2',system)
PE = Particle(pEcm,mE,'PE',system)
PF = Particle(pFcm,mF,'PF',system)

Pright= Particle(pmright,mright,'Pright',system)
Pleft=  Particle(pmleft,mleft,'Pleft',system)
```

### 4.6 Defining Valid Initial Condition of System:

Below is a defined valid initial condition for the system, this was retrieved from angle measurements within the solidworks model of the system and confirmed using the kinematics code within the zip folder this was submitted alongside.

The system starts at rest at the origin x1,x1_d=0, y1,y1_d=0

An additional initial condition was added that not in the solidworks model of qD2/qD2_d = 0 because the virtual joint will begin alined like the original link D.



```python
statevariables = system.get_state_variables()
# Initial "Guess" for state values
initialvalues = {}
initialvalues[x1]=0
initialvalues[x1_d]=0
initialvalues[y1]=0
initialvalues[y1_d]=0
initialvalues[qA]=-30*pi/180
initialvalues[qA_d]=0*pi/180
initialvalues[qB]=60*pi/180
initialvalues[qB_d]=0*pi/180
initialvalues[qC]=60*pi/180
initialvalues[qC_d]=0*pi/180
initialvalues[qD1]=60*pi/180   
initialvalues[qD1_d]=0*pi/180
initialvalues[qD2]=0*pi/180   
initialvalues[qD2_d]=0*pi/180
initialvalues[qE]=60*pi/180
initialvalues[qE_d]=0*pi/180
initialvalues[qF]=60*pi/180
initialvalues[qF_d]=0*pi/180
ini = [initialvalues[item] for item in statevariables]
```

### 4.7 Defining Angular Velocities based on Frame Relationships:

Here the angular velocities are defined similarly to how each frame is related to eachother.

A is rotated at angular velocity wNA from N

B is rotated at angular velocity wAB  from A

C is rotated at angular velocity wBC from B
"
D1 is rotated at angular velocity wCD1 from C

D2 is rotated at angular velocity wD1D2 from D1

E is rotated at angular velocity wD2E from D2

F is rotated at angular velocity wEF from E


```python
#Angular Velocities
wNA = N.getw_(A)
wAB = A.getw_(B)
wBC = B.getw_(C)
wCD1 = C.getw_(D1)
wD1D2 = D1.getw_(D2)
wD2E = D2.getw_(E)
wEF = E.getw_(F)

```

## 5. Adding Forces:

### 5.1 Creating Joint Stiffness/Spring Forces:

Spring forces are added to each joint to model the inherent elasticity and energy storage of each joint, this is calculated using the displacement of the state q"" minus the preload applied initially.

A different spring force is applied between the virtual joint between D1 D2 to mimic the inherent stiffness and compliance of the material.


```python
#Adding Spring/Joint Stiffness Forces
system.add_spring_force1(k,(qA-preload1)*N.z,wNA) 
system.add_spring_force1(k,(qB-preload2)*A.z,wAB)
system.add_spring_force1(k,(qC-preload3)*B.z,wBC)
system.add_spring_force1(k,(qD1-preload41)*C.z,wCD1)
system.add_spring_force1(k,(qE-preload5)*D2.z,wD2E)
system.add_spring_force1(k,(qF-preload6)*E.z,wEF)
#Virtual Joint Stiffness
system.add_spring_force1(kD,(qD2-preload42)*D1.z,wD1D2)
```




    (<pynamics.force.Force at 0x212f35ecca0>,
     <pynamics.spring.Spring at 0x212f35c9460>)



### 5.2 Creating Joint Damping Forces:

Below the damping forces are added for each joint proportional to the angular velocity of each joint.

There is a separate damping value for the virtual joint wD1D2 used when applying the damping force.


```python
#Adding Dampers
system.addforce(-b*wNA,wNA)
system.addforce(-b*wAB,wAB)
system.addforce(-b*wBC,wBC)
system.addforce(-b*wCD1,wCD1)
system.addforce(-bD*wD1D2,wD1D2)
system.addforce(-b*wD2E,wD2E)
system.addforce(-b*wEF,wEF)
#Adding Damper for Virtual Joint
system.addforce(-bD*wD1D2,wD1D2)
```




    <pynamics.force.Force at 0x212f3830370>



### 5.3 Adding Time Dependent Static Friction Force on Stoppers: (Actuator 1)

Here the stoppers raising and lowering causing static static friction are simulated by a large force damping the system at the stoppers locations. The static friction of the rubber on the stoppers will be high enough to prevent movement so this model should true.

For the code there is a force applied at the midpoint of pNA and PEF restricting movement until cutoff time tswitch. Then this force is released and a force is added at the then (fully extended) midpoint of pBC and pCD1 restricting movement.

This is used in combination with the torques discussed further below to create a "snake-like" rectilinear movement by:

Left Grabbing >Extending > Right Grabbing > Contracting

This is to be controlled by a motor actuating a string pulling up and down the stopped, but this would be difficult to model in this 2-D simulation, and wouldn't impact the robots motion.


```python
#Add static Friction force
bfrs= Constant(1000,'bfrs',system)

v1=((pNA+pEF)*.5).time_derivative(N,system)


tswitch = system.t-1
on = (tswitch+abs(tswitch))/(2*tswitch+1e-3)
off= abs(1-on)

v2 = ((pBC+pCD1)*.5).time_derivative(N,system)

system.addforce(-bfrs*v1*off,v1)
system.addforce(-bfrs*v2*on,v2)
```




    <pynamics.force.Force at 0x212f38297c0>



### 5.4 Adding Time Dependent Motor Torque: (Actuator 2)

Below are the motor torques applied to the system at point wCD1. Torque begins at t=0, expanding the linkage. Then the torque stops at time t=1 (tswitch) and reverses at time t=2 (tswitch2), contracting the linkage. 

A motor will interface with with link D using a rigid link attached to a slider on link D (Depicted in the Dynamic Figure Portion at the begining of this Document). This would allow actuation of the link via the motor rotation. Link D was made a compliant linkage due to the force this interface will apply on the linkage requiring compliance in the model. 

The interfacing of this actuated rigid link with link D was not modeled due to the difficulty of placing a slider along the linkage. A torque applied at the joint pCD1 with added compliance should be a suitible model for this actuation.

A torque constant was selected that fully extended the linkage and fell within the specs of the selected motor (Max torque .24 N*m)


```python
#Add motor force
Torque= Constant(.03,'Torque',system)

tswitch2 = system.t-2
on2 = (tswitch2+abs(tswitch2))/(2*tswitch2+1e-3)
off2= abs(1-on2)

system.addforce(Torque*N.z*off,wCD1)
system.addforce(-Torque*N.z*on2,wCD1)
```




    <pynamics.force.Force at 0x212f3749910>



### 5.5 Adding Directionally Dependent Kinetic Friction:

Below is the force of kinetic friction applied for each end cap dragging against the ground as the system is extending and retracting.

uk was found to be 0.3057 using methods found in the motor_powerSupply_Friction.pdf document.

Friction force for both end caps is calculated and shown below, but unfortunately was unable to be implemented due to running into computional problems, but the small force would likely have limited effect on this system.


```python
#Add kinetic friction force
uk=Constant(.3057,'uk',system)
Ffrright=9.81*uk*mrights
Ffrleft=9.81*uk*mlefts
#system.addforce(-uk*.25*pBC*off*N.x,v1)
#system.addforce(-uk*.25*pBC*on*-N.x,v1)
```


```python
Ffrright
```




    0.6701211⋅uk




```python
Ffrleft
```




    0.0442431⋅uk



### 5.6 Gravity's Role:

Because gravity is in the z plane of this model and has no effective impact on the motion of the model it was removed from the system. Gravity however was used to calculate the normal force external to the model and consequentially kinetic friction force.

#Gravity in -z direction
system.addforcegravity(-g*N.z)

## 6. Adding Constraints

### 6.1 Closing the Loop Using pFtip=pNA Constraints (Closing the Loop)

This constraint requires that the distance between pFtip and pNA is zero in the N.x and N.y directions. This closes the loop of the system and makes pFtip and pNA be coinciding.


```python
eq = []
eq.append((pFtip-pNA).dot(N.y))
eq.append((pFtip-pNA).dot(N.x))
```

### 6.2 C and F Parrallel Constraints

These constraints force C.x and F.x to be parrallel to the N.y axis. This is an idealized condition, but the model continiously ran into singularities and calculalation issues when making C.x and F.x parralel directly. 

With further testing, this constraint would be changed to C.x.dot(F) with another constraint to limit the degrees of freedom to 3 (saurrus linkage extension and (x,y) position).


```python
eq.append(C.x.dot(N.x))
eq.append(F.x.dot(N.x))
```

### 6.3 pNA and pBC Equal Y Component Constraints

This final constraint prevents the right half of the saurrus linkage from shifting up and down with the applied torque by restricting the movement of pBC up and down in the C.x direction from pNA.

#### Both this constraint and the constraints in 6.2 would not be placed in a 3-D simulation/real world model as they would be represented by a 3rd set of 2 links above the linkage applying these 3 constraints inherently.


```python
eq.append((pNA-pBC).dot(C.x))
```

### 6.4 Creating Constraint Derivatives

Derivatives of constraints to be used when solving for acceleration.


```python
eq_d=[(system.derivative(item)) for item in eq]
eq_dd=[(system.derivative(item)) for item in eq_d]
```

## 7. Solution:

### 7.1 Getting F and ma from System Dynamics

F,ma are retrieved from the system.


```python
#get system dyamics
f,ma = system.getdynamics()
```

    2021-03-22 12:31:12,821 - pynamics.system - INFO - getting dynamic equations
    

### 7.2 Getting Acceleration Equation by Dividing F by m

State acceleration function and lamdas are retrieved by dividing F by an inverted m.


```python
#solve for acceleration
func1,lambda1 = system.state_space_post_invert(f,ma,eq_dd,return_lambda = True)
```

    2021-03-22 12:31:18,425 - pynamics.system - INFO - solving a = f/m and creating function
    2021-03-22 12:31:18,517 - pynamics.system - INFO - substituting constrained in Ma-f.
    2021-03-22 12:32:06,487 - pynamics.system - INFO - done solving a = f/m and creating function
    2021-03-22 12:32:06,487 - pynamics.system - INFO - calculating function for lambdas
    

### 7.3 Integrating Accleration to Find System States

State acceleration function is integrated to determine system states.


```python
#integrate
states=pynamics.integration.integrate(func1,ini,t,rtol=tol,atol=tol, args=({'constants':system.constant_values},))


```

    2021-03-22 12:32:06,564 - pynamics.integration - INFO - beginning integration
    2021-03-22 12:32:06,566 - pynamics.system - INFO - integration at time 0000.00
    2021-03-22 12:33:16,364 - pynamics.system - INFO - integration at time 0001.00
    2021-03-22 12:34:24,858 - pynamics.system - INFO - integration at time 0001.15
    2021-03-22 12:35:36,290 - pynamics.system - INFO - integration at time 0002.00
    2021-03-22 12:36:38,241 - pynamics.integration - INFO - finished integration
    

### 7.4 Plotting State Values 

System State values are plotted x1 and y2 in m and qA,qB,qC,qD1,qD2,qE,qF in radians


```python
#Plot
plt.figure()
artists = plt.plot(t,states[:,:9])
plt.legend(artists,['x1','y1','qA','qB','qC','qD1','qD2','qE','qF'])
```




    <matplotlib.legend.Legend at 0x212f90810a0>




    
![png](SystemDynamicsII_files/SystemDynamicsII_81_1.png)
    


### 7.5 Plotting Total System Energy over Time

Total system energy is calculated by retrieving kinetic energy and potential energy from the system then summing them.

This is then plotted in Joules against each frame of the system.

System Energy is seen increasing due to work input by Acuator 2 actuating the linkage and decreasing due to the damping in each joint.


```python
#Energy
KE = system.get_KE()
PE = system.getPEGravity(pNA) - system.getPESprings()
energy_output = Output([KE-PE],system)
energy_output.calc(states)
energy_output.plot_time()
```

    2021-03-22 12:36:39,334 - pynamics.output - INFO - calculating outputs
    2021-03-22 12:36:39,667 - pynamics.output - INFO - done calculating outputs
    


    
![png](SystemDynamicsII_files/SystemDynamicsII_83_1.png)
    


### 7.6 Plotting Motion via Frames

Here the frames of the system are plotted over the time the system is simulated. The plot is made by lines between the points previously defined in the order pEF,pGH,pHI,pNA,pAB,pBC,pCD1,pKL,pJK,pBC,pCD1,pD1D2,pD2E,pEF,pFtip.

One can see the extending and grabbing and contracting in this diagram, but is easier to visualize in the animation within 7.7.


```python
#Motion

points = [pEF,pGH,pHI,pNA,pAB,pBC,pCD1,pKL,pJK,pBC,pCD1,pD1D2,pD2E,pEF,pFtip]

points_output = PointsOutput(points,system)
y = points_output.calc(states)
points_output.plot_time(20)
```

    2021-03-22 12:36:40,166 - pynamics.output - INFO - calculating outputs
    2021-03-22 12:36:40,369 - pynamics.output - INFO - done calculating outputs
    


    
![png](SystemDynamicsII_files/SystemDynamicsII_85_1.png)
    


### 7.7 Plotting Motion via Animation

Below the animation and final position of the robot is shown. Here one can see the robot extend switch stoppers then contract similar to a snake's rectilinear locomotion.

The compliant link D is noticable and shows how a link deforming may effect performance, especially if more virtual joints were applied to the system. This may indicate that we need to add another layer of material to the linkage, and will be tested in prototyping. 

It is also noticable how the system contracts twice, the first contraction is due to the removal of the applied torque allowing the system to "spring" back into its relaxed followed by the applied torque in reverse, contracting the system further.


```python
from matplotlib import animation, rc
from IPython.display import HTML
points_output.animate(fps = fps,movie_name = 'render.mp4',lw=2,marker='o',color=(1,0,0,1),linestyle='-')
HTML(points_output.anim.to_html5_video())
```




<video width="432" height="288" controls autoplay loop>
  <source type="video/mp4" src="data:video/mp4;base64,AAAAIGZ0eXBNNFYgAAACAE00ViBpc29taXNvMmF2YzEAAAAIZnJlZQAAcoxtZGF0AAACrgYF//+q
3EXpvebZSLeWLNgg2SPu73gyNjQgLSBjb3JlIDE2MSByMzA0OCBiODZhZTNjIC0gSC4yNjQvTVBF
Ry00IEFWQyBjb2RlYyAtIENvcHlsZWZ0IDIwMDMtMjAyMSAtIGh0dHA6Ly93d3cudmlkZW9sYW4u
b3JnL3gyNjQuaHRtbCAtIG9wdGlvbnM6IGNhYmFjPTEgcmVmPTMgZGVibG9jaz0xOjA6MCBhbmFs
eXNlPTB4MzoweDExMyBtZT1oZXggc3VibWU9NyBwc3k9MSBwc3lfcmQ9MS4wMDowLjAwIG1peGVk
X3JlZj0xIG1lX3JhbmdlPTE2IGNocm9tYV9tZT0xIHRyZWxsaXM9MSA4eDhkY3Q9MSBjcW09MCBk
ZWFkem9uZT0yMSwxMSBmYXN0X3Bza2lwPTEgY2hyb21hX3FwX29mZnNldD0tMiB0aHJlYWRzPTkg
bG9va2FoZWFkX3RocmVhZHM9MSBzbGljZWRfdGhyZWFkcz0wIG5yPTAgZGVjaW1hdGU9MSBpbnRl
cmxhY2VkPTAgYmx1cmF5X2NvbXBhdD0wIGNvbnN0cmFpbmVkX2ludHJhPTAgYmZyYW1lcz0zIGJf
cHlyYW1pZD0yIGJfYWRhcHQ9MSBiX2JpYXM9MCBkaXJlY3Q9MSB3ZWlnaHRiPTEgb3Blbl9nb3A9
MCB3ZWlnaHRwPTIga2V5aW50PTI1MCBrZXlpbnRfbWluPTI1IHNjZW5lY3V0PTQwIGludHJhX3Jl
ZnJlc2g9MCByY19sb29rYWhlYWQ9NDAgcmM9Y3JmIG1idHJlZT0xIGNyZj0yMy4wIHFjb21wPTAu
NjAgcXBtaW49MCBxcG1heD02OSBxcHN0ZXA9NCBpcF9yYXRpbz0xLjQwIGFxPTE6MS4wMACAAAAR
62WIhAAz//727L4FNf2f0JcRLMXaSnA+KqSAgHc0wAAAAwAAeB0oXug4SqvQwAAbcAB3KdixJz/8
AExk9nC766hr8cyqauS98M6la4FDcFcmzXw10zHrLzON5s8ZUOP0ubLvkNoPOLOUufuzVS6t5I0Z
PJeD89n4Rp8Wjfe3Dk6k+PgwqbcylnmvDqnDiJfcCRj1CN9xyC9oPlaAAu1tF2lL1d62Mzah8FHN
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
AAAAAQAABAAAAAABAAAAAAAAAAEAAAIAAAAAAQAACgAAAAABAAAEAAAAAAEAAAAAAAAAAQAAAgAA
AAABAAAKAAAAAAEAAAQAAAAAAQAAAAAAAAABAAACAAAAAAEAAAoAAAAAAQAABAAAAAABAAAAAAAA
AAEAAAIAAAAAAQAACgAAAAABAAAEAAAAAAEAAAAAAAAAAQAAAgAAAAAEAAAEAAAAAAEAAAoAAAAA
AQAABAAAAAABAAAAAAAAAAEAAAIAAAAAAQAACAAAAAACAAACAAAAAAEAAAoAAAAAAQAABAAAAAAB
AAAAAAAAAAEAAAIAAAAAAgAABAAAAAABAAAGAAAAAAEAAAIAAAAAAQAACgAAAAABAAAEAAAAAAEA
AAAAAAAAAQAAAgAAAAABAAAKAAAAAAEAAAQAAAAAAQAAAAAAAAABAAACAAAAAAEAAAgAAAAAAgAA
AgAAAAABAAAEAAAAAAEAAAYAAAAAAQAAAgAAAAABAAAIAAAAAAIAAAIAAAAAAQAACgAAAAABAAAE
AAAAAAEAAAAAAAAAAQAAAgAAAAADAAAEAAAAAAEAAAoAAAAAAQAABAAAAAABAAAAAAAAAAEAAAIA
AAAAAQAACgAAAAABAAAEAAAAAAEAAAAAAAAAAQAAAgAAAAABAAAKAAAAAAEAAAQAAAAAAQAAAAAA
AAABAAACAAAAAAEAAAoAAAAAAQAABAAAAAABAAAAAAAAAAEAAAIAAAAAAQAACgAAAAABAAAEAAAA
AAEAAAAAAAAAAQAAAgAAAAABAAAKAAAAAAEAAAQAAAAAAQAAAAAAAAABAAACAAAAAAEAAAoAAAAA
AQAABAAAAAABAAAAAAAAAAEAAAIAAAAAAQAACgAAAAABAAAEAAAAAAEAAAAAAAAAAQAAAgAAAAAB
AAAKAAAAAAEAAAQAAAAAAQAAAAAAAAABAAACAAAAAAEAAAoAAAAAAQAABAAAAAABAAAAAAAAAAEA
AAIAAAAAAQAACgAAAAABAAAEAAAAAAEAAAAAAAAAAQAAAgAAAAABAAAKAAAAAAEAAAQAAAAAAQAA
AAAAAAABAAACAAAAAAEAAAoAAAAAAQAABAAAAAABAAAAAAAAAAEAAAIAAAAAAQAACgAAAAABAAAE
AAAAAAEAAAAAAAAAAQAAAgAAAAABAAAKAAAAAAEAAAQAAAAAAQAAAAAAAAABAAACAAAAAAEAAAoA
AAAAAQAABAAAAAABAAAAAAAAAAEAAAIAAAAAAQAACgAAAAABAAAEAAAAAAEAAAAAAAAAAQAAAgAA
AAABAAAKAAAAAAEAAAQAAAAAAQAAAAAAAAABAAACAAAAAAEAAAgAAAAAAgAAAgAAAAAcc3RzYwAA
AAAAAAABAAAAAQAAAJYAAAABAAACbHN0c3oAAAAAAAAAAAAAAJYAABShAAAEtwAAA9cAAAM1AAAE
LgAAAQcAAAL1AAAAtAAAApMAAACLAAAAYgAAAFMAAAFmAAAAWgAAADEAAAAlAAAAxQAAADYAAAAk
AAAAIwAAAFIAAAAtAAAAIQAAAB4AAABAAAAAKgAAAB0AAAAcAAACJwAAAD8AAAAdAAAALAAAAtsA
AAMOAAAC8AAAAkkAAAQCAAABCAAAAOcAAAC+AAACJQAAAJoAAABiAAAB4AAAAFQAAAA6AAAAIAAA
AJgAAACpAAAAfQAAACUAAABpAAAAJgAAAB0AAAAVAAAASgAAACMAAAAYAAAAFgAAAssAAABHAAAA
LwAAAyEAAAU1AAAA9gAABQAAAAETAAAAqAAAAzYAAADdAAAAjgAAADkAAAFRAAAA5AAAAMIAAACt
AAAAIgAAABwAAAAXAAAAHwAAABwAAAAXAAAAFwAAACYAAAAcAAAAFwAAABcAAAAzAAAAHAAAABcA
AAAXAAAAGAAAABwAAAAXAAAAFwAAABgAAAAcAAAAFwAAABcAAAAYAAAAHAAAABcAAAAXAAAAGAAA
ABwAAAAXAAAAFwAAAC8AAAAeAAAAGQAAABcAAAAvAAAAHAAAABcAAAAXAAAAGAAAABwAAAAXAAAA
FwAAABgAAAAcAAAAFwAAABcAAAAYAAAAHAAAABcAAAAXAAAAGAAAABwAAAAXAAAAFwAAABgAAAAc
AAAAFwAAABcAAAAYAAAAHAAAABcAAAAXAAAAFwAAABwAAAAXAAAAFwAAABcAAAAcAAAAFwAAABcA
AAAYAAAAHAAAABcAAAAUc3RjbwAAAAAAAAABAAAAMAAAAGJ1ZHRhAAAAWm1ldGEAAAAAAAAAIWhk
bHIAAAAAAAAAAG1kaXJhcHBsAAAAAAAAAAAAAAAALWlsc3QAAAAlqXRvbwAAAB1kYXRhAAAAAQAA
AABMYXZmNTguNjguMTAw
">
  Your browser does not support the video tag.
</video>




    
![png](SystemDynamicsII_files/SystemDynamicsII_87_1.png)
    


### 7.8 Additional Plots of Energy Decreasing when Actuators are off and Initial Velocity is Given

Initial Angular Velocity wNA=5rad/s

Here the natural damping of the system is shown via the energy decrease in the total energy plot when an initial angular velocity is applied with actuators off.

States Plot with Actuators Off(y: m for x1,y1 rad for qA-F, x: time s)

![noacuationstates.png](attachment:noacuationstates.png)

Total Energy Plot with Actuators Off (y: Joules x: Frames)

![noacuationenergy.png](attachment:noacuationenergy.png)

Animation with Initial Velocity and No Actuators On

![noactuationanimation.gif](attachment:noactuationanimation.gif)

### 7.9 Equation of States and Derivatives (Commented Out)

Here are the equations describing each state and their derivatives commented out to reduce computation time.


```python
#eq
```


```python
#eq_d
```


```python
#eq_dd
```

## 8. Discussion

### 8.1 Model Assumptions, Idealizations, and Limitations

There are quite a few assumptions made within this model, but are well-grounded and should leave this model to accurately depict movement in real-world application. The primary limitation of this model is that it is mostly represented in 2 dimensions, this requires constraints to be used instead of links for the 3rd set of 2 linkages that form a 1 degree of freedom saurus linkage. The 2-D simulation also makes it difficult to simulate the snake robot traversing uphill, which is one of the primary uses of rectilinear locomotion in snakes.

Some other assumptions that should yield very little consequence when compared to real-world performance. These include the lumped masses and particle modeling. The lumped masses should not differ from a simulation using precisely defined masses, as their locations would not effect the dynamics as long as they are located on the correct side of the saurrus linkage. This also applies to their inertias because the ends of the robot should not be turning for the rotational inertias to impact the dynamics. The links however do rotate, but their rotational inertias are negligible due to how light the links are and how slow they are rotating. Other assumptions include assuming static friction would hold, no compliance in 5 of the links, and multiple assumptions made during the collection of damping and stiffness data collection.



### 8.2 Discussion of Difficulties
Two issues were unable to be resolved after being able to overcome other difficulties that appeared throughout the creation of this simulation. The damping and stiffness values were unable to be calculated due to difficulties translating the collected data into the damping and stiffness values.

Another issue was running into computation difficulties when implementing kinetic friction. The kinetic friction was calculated, but adding the time dependant component ended up taking the simulation a very long time to run (>20mins) and could be resolved in the future by optimizing the model with 2 3 bar linkages as opposed to the 1 6 bar linkage currently used. 

### 8.3 Discussion of Results

One of the main results that can be seen in this experiment is one of the main trials to face when using foldable robotics, the effect compliance in materials used. It can be seen that the angles of each side of the saurus linkage are different and the point in the middle (pD2E) is closer to the center than the opposite side (pAB). This could cause the robot to tilt slightly and start moving in a path that is not straight. This means that special attention will need to be paid attention to link stiffness when constructing the prototype.

Another main result that can be seen is the inherent elasticity and damping of each joint in the system. These joints lose some energy due to damping and cause quite a bit of movement due to elasticity (springiness) that will need to be factored in when planning the movement or programing the actuator control for the system. 

Lastly, the animation displays that we were able to successfully simulate the rectilinear motion we were attempting to model using foldable robotics principles and dynamics calculations. This simulation will be useful to test what-if scenarios in case we wanted to see how different materials or adding different masses to the system might affect it without constructing or altering a prototype each time.

### 8.4 Future Plans for Simulations

In the future, a fully functional prototype will be constructed to compare to the results of this simulation. To assist with this simulation will also be optimized to run faster using multiple chains of frames as opposed to a single chain from A to F. With a more optimized simulation more compliance in linkages will be added to determine the effect of compliance on every link in the system will have. Also hopefully the interface of the actuators with their actuated components robot will be improved. 

## 9. Bibliography

[1] D.Aukes, “Approximating compliant beams with the pseudo-rigid-body model,”egr557. Accessed: March. 19, 2021 [Online]. Available: https://egr557.github.io/modules/compliance/generated/prbm.html

[2] D.Aukes, “Cantilever Beams,”egr557. Accessed: March. 19, 2021 [Online]. Available: https://egr557.github.io/modules/compliance/cantilever-beams.html 

[3] D.Aukes, “Euler-Bernoulli Beams,”egr557. Accessed: March. 19, 2021 [Online]. Available: https://egr557.github.io/modules/compliance/generated/euler-bernoulli-beams/euler-bernoulli-beams.html 

[4] “Servo Motor SG-90 Basics, Pinout, Wire Description, Datasheet,” Components101. [Online]. Available: https://components101.com/motors/servo-motor-basics-pinout-datasheet. [Accessed: 20-Mar-2021]. 
