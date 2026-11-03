This is an implementation agenda for future developments of PPFM.

1) Re-factoring of all of the I/O features with Json files, tags, and metadata.

2) Implementation of the Hard-sphere approximation in order to always provide a very first approximation of transport coefficients for 

3) Re-factoring of the code to use smart pointers in order to avoid memory leaks

4) Re-factoring of the template classes, avoid using variants so that the build does no longer depend on the main file and on file AcceptedSpeciesGuard.py

5) Re-factoring of the transport module: 
    - bracket integrals matrices on initialization
    - possibility to include or not the reactive and internal contribution to thermal conductivity
    - ease the customization of the order of transport coefficients
    - Implement 2T reactive contributions to thermal conductivity 

6) Implement at least the rigid rotor and harmonic oscillator approximations for molecular rotational and vibrational partition functions.