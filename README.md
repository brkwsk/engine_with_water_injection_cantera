# engine with water injection python cantera
This program allows simplified simulations of 0-dimensional internal combustion engine's work with water injection.

The code is a modified version of an example aviable on Cantera's page: 
https://cantera.org/examples/python/reactors/ic_engine.py.html

A large part of the comments is in Polish, which I plan to change in the future for easier reading.

The parameters for the simulation have to be set manually as values for the variables.

The program's output is a printed line of parameters for every 10th time step, printed condition of injected water (just to check that everything is in order), a file "wyniki.txt", containing each line, the maximum value of pressure, and p-V and p-CA plots.

The simulation runs for 3 cycles of engine work.
