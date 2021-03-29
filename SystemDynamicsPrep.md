<p style="color: red; font-weight: bold">>>>>>  gd2md-html alert:  ERRORs: 0; WARNINGs: 0; ALERTS: 3.</p>
<ul style="color: red; font-weight: bold"><li>See top comment block for details on ERRORs and WARNINGs. <li>In the converted Markdown or HTML, search for inline alerts that start with >>>>>  gd2md-html alert:  for specific instances that need correction.</ul>

<p style="color: red; font-weight: bold">Links to alert messages:</p><a href="#gdcalert1">alert1</a>
<a href="#gdcalert2">alert2</a>
<a href="#gdcalert3">alert3</a>

<p style="color: red; font-weight: bold">>>>>> PLEASE check and correct alert issues and delete this message and the inline alerts.<hr></p>


1b) Motor Selection

 

SG90 Specification


    ●	Rotation : 0°-180°


    ●	Weight of motor : 9gm


    ●	Operating Voltage : +5 V


    ●	Stall Current : 350 mA


    ●	Stall Torque : 2.5 kg . cm = 0.24 N.m


    ●	Speed : 0.12s / 60° = 82.8 rpm


    ●	Power : Torque * Speed / 9.5488 = 2.08 kW

This motor was selected due to the requirement of having a very lightweight actuator which did not lack the torque to move the entire load of the robot.

 

1c) Power Supply Selection:

 

Given the stall current of the motors being used, a power supply capable of supplying far more than 350 mA is required. An Alkaline 9v battery is chosen as the power supply for the system which satisfies this requirement and has enough battery capacity to sustain the robot to operate for a reasonable span of time.

 

 

 

 

 

 

 

 

 

 

 

Static Friction Calculation:

 

The Angle of Repose experiment was conducted to determine the static friction coefficient of the cardstock on a smooth surface. The platform having the cardstock is lifted slowly with increasing angle to find that critical angle at which the paper starts sliding down.

 

Trial 1:



<p id="gdcalert1" ><span style="color: red; font-weight: bold">>>>>>  gd2md-html alert: inline image link here (to images/image1.png). Store image on your image server and adjust path/filename/extension if necessary. </span><br>(<a href="#">Back to top</a>)(<a href="#gdcalert2">Next alert</a>)<br><span style="color: red; font-weight: bold">>>>>> </span></p>


![alt_text](images/image1.png "image_tooltip")


Coefficient of Static Friction : 0.3141

 

 

 

 

 

 

 

 

 

 

Trial 2 :



<p id="gdcalert2" ><span style="color: red; font-weight: bold">>>>>>  gd2md-html alert: inline image link here (to images/image2.png). Store image on your image server and adjust path/filename/extension if necessary. </span><br>(<a href="#">Back to top</a>)(<a href="#gdcalert3">Next alert</a>)<br><span style="color: red; font-weight: bold">>>>>> </span></p>


![alt_text](images/image2.png "image_tooltip")


Coefficient of Static Friction : 0.4575

 

Trial 3 :

 

<p id="gdcalert3" ><span style="color: red; font-weight: bold">>>>>>  gd2md-html alert: inline image link here (to images/image3.png). Store image on your image server and adjust path/filename/extension if necessary. </span><br>(<a href="#">Back to top</a>)(<a href="#gdcalert4">Next alert</a>)<br><span style="color: red; font-weight: bold">>>>>> </span></p>


![alt_text](images/image3.png "image_tooltip")


Coefficient of Static Friction : 0.3803

Average Coefficient of Static Friction : 0.38399

 

---------------------------------------------------------------------------------------------

 

Average Coefficient of Static Friction : 0.3057

 

Approximation : The coefficient of kinetic friction could be determined during the same Angle of Repose experiment where the angle corresponding to when the cardstock slips down the inclined surface at constant velocity is noted. An angle of ~17 degrees let the cardstock slip down in an approximately constant velocity. Thus coefficient of kinetic friction = tan (17) = 0.3057

