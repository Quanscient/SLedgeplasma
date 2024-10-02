# SLedgeplasma
---

This is a finite element based 1D edge plasma (SOL) solver. The code implementation is currently carried out through example simulations.

The repository contains the following simulation examples:
- A_isothermalmodel
- B_selfconsistenttemperature
- Case01
- Case02
- Case03
- Case04

# Running a simulation

To run the simulation, the [Sparselizard](https://github.com/halbux/sparselizard) library is required. Follow the instructions in the library to configure and build it.

After building the Sparselizard library, do the following to run any of the examples:
1. Copy the example folder, say `Case01`, to `sparselizard/simulations`.
1. Add line `add_subdirectory(Case01)` to `simulations/CMakeLists.txt`.
1. Configure and build. Executable file will be located in build/simulations/newsim folder.