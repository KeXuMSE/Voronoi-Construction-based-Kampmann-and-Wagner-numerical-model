# Voronoi Construction based-Kampmann and Wagner numerical (KWN) model
This software package was developed for predicting the precipitation kinetics in multi-component alloy systems, based on Kampmann and Wagner numerical (KWN) model and CALPHAD. Nucleation is implemented utilizing the classical nucleation theory (CNT). Growth and coarsening are modeled by a single growth kinetics equation, which is constructed based on the interfacial diffusion flux balance and the capillarity effect. A new feature of the model is the incorporation of a more realistic spatial site distribution via Voronoi construction in one characteristic cell. Two approaches are presented here, one is the classical Lagrange like PSD implementation, and the other one is the modified version based on Voronoi construction.
## Highlights of the model
* Accurate predictions of dimensions and composition of the precipitates
* Ability to handle sophisticated alloy chemistries by a simple and quantitative growth kinetics equation
* Visualization of the spatial distribution of precipitates via Voronoi construction 
## How to run the code
### 1. Install [Thermo-Calc](https://thermocalc.com/products/thermo-calc/) and [TQ-interface](https://thermocalc.com/products/software-development-kits/)
### 2. Install [Multi-Parametric Toolbox](https://www.mpt3.org/) in Matlab
### 3. Build the main framework and extract thermodynamics data from Thermo-Calc in Fortran with Visual Studio
#### 3.1. Configuration of Fortran for TQ-interface in Visual Studio
* Set up application builds to 64-bit platforms.
* Project >> Properties >> MIDL >> General >> Target Environment >> Microsoft Windows 64-bit on x64 (/env x64)
* Project >> Properties >> Fortran >> General >> Additional Include Directories >> C:\Users\Public\Documents\Thermo-Calc\2017b\SDK\TQ
* Project >> Properties >> Fortran >> Data >> Default Integer KIND >> 8 (/integer_size:64)
* Project >> Properties >> Fortran >> Data >> Default Real KIND >> 8 (/real_size:64)
* Project >> Properties >> Linker >> General >> Additional Library Directories >> C:\Users\Public\Documents\Thermo-Calc\2017b\SDK\TQ
* Project >> Properties >> Linker >> Input >> Additional Dependencies >> tq-win-x64-X.X.XXXXX.lib (license of TQ-interface)
* Project >> Properties >> Linker >> System >> Sub System >> Console (/SUBSYSTEM:CONSOLE)

### 4. Hybrid-programming of Fortran and Matlab
#### 4.1. Set the Environment Variables in Windows system  
* Right click on My Computer >> Properties >> Advanced system settings >> Advanced >> Environment Variables >> find PATH in System Variables, double click and add new PATH, “C:\Program Files\Matlab\R2017b\bin\win64”
#### 4.2. Configuration of Fortran for Matlab programming in Visual Studio  
* Project >> Properties >> Fortran >> General >> Additional Include Directories >> add “C:\ProgramFiles\Matlab\R2017b\extern\include”  
* Project >> Properties >> Fortran >> Preprocessor >> Preprocess Source File >> choose “Yes”  
* Project >> Properties >> Linker >> General>> Additional Library Directories >> add “C:\ProgramFiles\Matlab\R2017b\extern\lib\win64\microsoft”   
* Project >> Properties >> Linker >> Input >> Additional Dependencies >> add “libmx.lib”, “libmat.lib” and “libeng.lib”  
#### 4.3. Open Matlab and set path in Fortran framework
* Open Matlab
```
      ep=engOpen('')
      if (ep == 0) then
          write(*,*) 'can''t start matlab engine'
          stop
      end if
```
* Set path
```
      if (engEvalString(ep," addpath('C:\Program Files\MATLAB\R2017b\toolbox\mpt\mpt')") /= 0) then
          write(*,*) 'engEvalString failed'
          stop
      end if
```
## Representative results
* The evolution of spatial distribution of precipitates and Voronoi cells within the characteristic cell in alloy Ni-7.5Al-8.5Cr at. % during isothermal ageing
![image](https://github.com/KeXuMSE/Voronoi-Construction-based-Kampmann-and-Wagner-numerical-model/blob/main/Fig1.png)
## License
The code is released under the MIT license.
## Reference
The paper is still under review. Code will be posted online after final accept. Please cite our work if getting inspired! Thank you!
