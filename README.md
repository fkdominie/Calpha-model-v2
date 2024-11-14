# Calpha-model

This C code implements a one bead-per-amino acid coarse-grained
model for protein folding, with a structure-based potential.
Conformational sampling is carried out using Langevin dynamics.

For a detailed description of the model and sampling method,
see:

S Bazmi, B Seifi, S Wallin. Simulations of a protein fold switch
reveal crowding-induced population shifts driven by disordered
regions. Communications Chemistry 6 191 (2023).

S Wallin, HS Chan. Conformational entropic barriers in topology-
dependent protein folding: perspectives from a simple native-
centric polymer model. Journal of Physics: Condensed Matter 18
S307-S328 (2006).

A minimal sequence of commands to compile the code and run a
simulation of a single chain at a fixed temperature is:

make constants

./constants input 1

make fixtemp

./main

Contact: Stefan Wallin, swallin@mun.ca

