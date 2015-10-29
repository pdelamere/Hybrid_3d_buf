# Hybrid space plasma simulator
It only simulates Pluto, but that may change in the future.

## Compiling
First build and install OpenMPI. We're using version 1.4. You might be able to ge more recent versions to work, but I've had problems with 1.8.

To build the hybrid code just run `make`
This will produce a binary called hybrid that can be executed using `mpirun -n 24`
Right now it works well with the Intel compiler, but we're working on improving portability.

## Usage
The most commonly used settings are in inputs.dat, and dimensions.f.
Yes you will have to recompile if you change the dimensions of the simulation.

Right now that's about all there is to it. It simulates Pluto, if you want something else you 
have to dig into the codebase and edit it. There's currently no convenient way to change
most simulation parameters.

In fact it probably won't work at all unless you use 24 processes. (i.e. `mpirun -n 24 hybrid`)

If you get a segmentation fault shortly after startup try increasing the shell stack size. In bash or zsh try `ulimit -s unlimited` some systems have a hard limit on stack size. In that case try `ulimit -s 65532`

## Organization
Program starts in maind.f, but the bulk of the code is in gutsp.f.
Gutsp.f is the guts of the program for the kinetic simulation.
Gutsf.f is the guts of the program for the fluid simulation.
Comments in maind.f explain many of the variable names.

## Contact
Project lead: Peter.Delamere@gi.alaska.edu
