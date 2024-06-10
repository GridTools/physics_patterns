In this example, the turbulent mixing scale is computed according to the distance of the layer to the
tropopause.

- l91 to l98 : the vertical index _JK_ of the tropopause level is retrieved from temperature inversion and stored in _ITPL(IJ)_ 
- l101 to 113 : the mixing length _ZL_ is computed according to height of the tropopause _PZZ(IJ, JKP)_ 


Indices conventions are the following :
- IKL : order of the vertical levels
    1 if level are numbered from ground to space
    -1 if level are numbered from space to ground
- IKB : nearest level to the ground
- IKE : nearest level to the space
- IKTB is IKB + padding
- IKTE is IKE + padding

In this case, loop at l101 is always executed from ground to space