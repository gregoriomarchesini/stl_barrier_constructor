# STL Barrier Functions Constructor 

This repo contains the code to create a set of time-varying barrier functions for a multi-agent system tasked with the achievement of a time-varying formation dictated by a set of Signal temporal Logic tasks.

## Installation

### Mosek License
In order to run the tool a [MOSEK](https://www.mosek.com/) license is necessary. This is freely available for students. We recomment following the instructions [here](https://www.mosek.com/products/academic-licenses/) to get your license 

### Download the mosek binaries 
Download the mosek binaris for your platform [here](https://www.mosek.com/downloads/#:~:text=Default%20Installers%20(Preferred)) and place the extracted mosek folder in the lib directory such that you have 

```
lib/mosek/11/..
```

### Compile
Compile the code by running 

```
mkdir build
cd build
cmake ..
make
```




