# Copyright
```
//================================================================================================//
//------------------------------------------------------------------------------------------------//
//    MPH-ElasticPlate : Moving Particle Hydrodynamics for Elastic Plate                          //
//------------------------------------------------------------------------------------------------//
//    Developed by    : Masahiro Kondo                                                            //
//    Distributed in  : 2023                                                                      //
//    Lisence         : GPLv3                                                                     //
//    For instruction : see README                                                                //
//    For theory      : see the following references                                              //
//     [1] JSCES Paper No.20100016,   https://doi.org/10.11421/jsces.2010.20100016                //
//    Copyright (c) 2010  Masahiro Kondo                                                          //
//    Copyright (c) 2023  National Institute of Advanced Industrial Science and Technology (AIST) //
//================================================================================================//

This program is for conducting particle based elastic plate simulation
using MPH-ElasticPlate (Moving Particle Hydrodynamics for Elastic Plate) method. 
```

# Directories
The directories in the repository is as follows:  
```
MphExplicit ---- generator (not ready)
            |--- results
            |--- source
```          

# How to execute samples

1. "source" contains the main solver for MPH-Elastic calculations. 
To complie it, run
```
> make 
```
in the source directory. 

2. Then, in the case directory (results/sample/), launch
```
> ./execute.sh
```

The solver starts with reading parameter file (*.data) 
and particle file (*.grid). 
The series of calculation results (*.vtk) are output,  
and they can be visualized with a viewer like ParaView. 
(https://www.paraview.org/)


# Using openMP and openACC (not ready)


# Generating particles (not ready)


# Changing physical properties
The physical properties and calculation conditions are given in 
the parameter file (*.data). By editting the file, the physical 
properties can be changed. 





