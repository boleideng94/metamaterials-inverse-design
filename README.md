# metamaterials-inverse-design

***discrete_model*** contains MATLAB code for discrete model simulation. It takes 16 design parameters of the geometry, generate the 10by8 metamaterial and simulate its response under uniaxial compression.

***NAES*** contains jupyter notebook for neural accelerated evolution strategy (NAES). In the notebook, a neural network (NN) is first trained based on training data and then a evoltuion strategy (ES) is used to conduct inverse design. 

***dataset*** contains input geometrical data and output stress data from 30,000 randomly generated caeses.

***autocad_script*** contains MATLAB code that generates AUTOCAD script. User can input any unit cell design and the MATLAB code will generate AUTOCAD script of the mold, which we used to physically fabricate our samples. 
