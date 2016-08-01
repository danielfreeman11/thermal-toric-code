# thermal-toric-code
A collection of programs which analyze finite temperature behavior and error correction strategies for Kitaev's toric code (and some related stabilizer codes).

In approximate chronological order (as they were developed by me):

1. TCPython contains code for performing finite temperature simulation of the dynamics of Alexei Kitaev's toric code.  This code is now deprecated, but it was what I used to generate all of the data in http://arxiv.org/pdf/1405.2315.pdf

2. CyclingAlgorithm contains code for finite temperature simulation + a simple measurement-free error correction protocol for the 1D Ising model (a 1D analogue of the toric code).  This was used to generate the data in: http://arxiv.org/pdf/1603.05005.pdf

3. FastCpp is functionally the same code as TCPython, but entirely in C++.  It is currently being used to investigate entanglement dynamics (among other things) in the toric code.

4. ToricLearning is my new foray into the world of machine learning, where I'm trying to build an OpenAI-gym-style environment for error-correcting the toric code.
