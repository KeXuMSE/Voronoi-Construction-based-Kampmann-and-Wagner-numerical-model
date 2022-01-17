# Voronoi Construction based-Kampmann and Wagner numerical (KWN) model
This software package was developed for predicting the precipitation kinetics in multi-component alloy systems, based on Kampmann and Wagner numerical (KWN) model and CALPHAD. Nucleation is implemented utilizing the classical nucleation theory (CNT). Growth and coarsening are modeled by a single growth kinetics equation, which is constructed based on the interfacial diffusion flux balance and the capillarity effect. A new feature of the model is the incorporation of a more realistic spatial site distribution via Voronoi construction in one characteristic cell. Two approaches are presented here, one is the classical Lagrange like PSD implementation, and the other one is the modified version based on Voronoi construction.
## Highlights of the model
* Accurate predictions of dimensions and composition of the precipitates
* Ability to handle sophisticated alloy chemistries by a simple and quantitative growth kinetics equation
* Visualization of the spatial distribution of precipitates via Voronoi construction 
## How to run the code
#### Install [Thermo-Calc](https://thermocalc.com/products/thermo-calc/) and [TQ-interface](https://thermocalc.com/products/software-development-kits/)
#### Install [Multi-Parametric Toolbox](https://www.mpt3.org/) in Matlab
#### Build the main framework and extract thermodynamics data from Thermo-Calc in Fortran with Visual Studio
#### Hybrid-programming of Fortran and Matlab
* set the Environment Variables in Windows system
* set the Property of Fortran in Visual Studio
## Representative results
* The evolution of spatial distribution of precipitates and Voronoi cells within the characteristic cell in alloy Ni-7.5Al-8.5Cr at. % during isothermal ageing
![image](https://github.com/KeXuMSE/Voronoi-Construction-based-Kampmann-and-Wagner-numerical-model/blob/main/Fig1.png)
## License
The code is released under the MIT license.
## Reference
The paper is still under review. Code will be posted online after final accept. Please cite our work if getting inspired! Thank you!
