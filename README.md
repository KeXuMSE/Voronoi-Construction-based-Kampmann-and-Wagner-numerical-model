# Voronoi Construction based-Kampmann and Wagner numerical (KWN) model
This software package was developed for predicting the precipitation kinetics in multi-component alloy systems, based on Kampmann and Wagner numerical (KWN) model and CALPHAD. Nucleation is implemented utilizing the classical nucleation theory (CNT). Growth and coarsening are modeled by a single growth kinetics equation, which is constructed based on the interfacial diffusion flux balance and the capillarity effect. A new feature of the model is the incorporation of a more realistic spatial site distribution via Voronoi construction in one characteristic cell. Two approaches are presented here, one is the classical Lagrange like PSD implementation, and the other one is the modified version based on Voronoi construction.
## Highlights of the model
* Accurate predictions of dimensions and composition of the precipitates
* Ability to handle sophisticated alloy chemistries by a simple and quantitative growth kinetics equation
* Visualization of the spatial distribution of precipitates via Voronoi construction 
## How to run the code
### Install [Thermo-Calc](https://thermocalc.com/products/thermo-calc/) and [TQ-interface](https://thermocalc.com/products/software-development-kits/)
### Install [Multi-Parametric Toolbox](https://www.mpt3.org/) in Matlab
### Build the main framework and extract thermodynamics data from Thermo-Calc in Fortran with Visual Studio
### Hybrid-programming of Fortran and Matlab
#### Set the Environment Variables in Windows system  
* Right click on My Computer >> Properties >> Advanced system settings >> Advanced >> Environment Variables >> find PATH in System Variables, double click and add new PATH, “C:\Program Files\Matlab\R2010b\bin\win64”
#### Set the Properties of Fortran in Visual Studio  
* Project >> Properties >> Fortran >> General >> Additional Include Directories >> add “C:\ProgramFiles\Matlab\R2010b\extern\include”  
* Project >> Properties >> Fortran >> Preprocessor >> Preprocess Source File >> choose “Yes”  
* Project >> Properties >> Linker >> General>> Additional Library Directories >> add “C:\ProgramFiles\Matlab\R2010b\extern\lib\win64\microsoft”   
* Project >> Properties >> Linker >> Input >> Additional Dependencies >> add “libmx.lib”, “libmat.lib” and “libeng.lib”  
#### Open Matlab and set path in Fortran framework
* Open Matlab
'''
      ep=engOpen('')
      if (ep == 0) then
          write(*,*) 'can''t start matlab engine'
          stop
      end if
'''
* set path
## Representative results
* The evolution of spatial distribution of precipitates and Voronoi cells within the characteristic cell in alloy Ni-7.5Al-8.5Cr at. % during isothermal ageing
![image](https://github.com/KeXuMSE/Voronoi-Construction-based-Kampmann-and-Wagner-numerical-model/blob/main/Fig1.png)
## License
The code is released under the MIT license.
## Reference
The paper is still under review. Code will be posted online after final accept. Please cite our work if getting inspired! Thank you!
