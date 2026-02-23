The following is a list of things that are "wrong" or suboptimal things in this code base and could/should be changed 
by the next person.

## 1. In varsub-PionPhotoProdThresh.twobody/varsub-2Bkernel.PionPhotoProdThresh.f:
  **The parameter `Jacobian` shouldn't be multiplying the diagram contributions by -1.**
  Since this is an integral over all space (twice) you actually need
  to take the absolute value of the Jacobian, so it should be just +1 (not -1 as it is coded). 
  It has not been changed because the file parsing code has been generated in such a way that compensates for this, 
  and there have been many calculations done with this method, and it is too close to the end of my thesis to redo all of this.

## 2. In the varsub-density-modules
  The densities come in a grid with momentum values $P=\{p_i\}$ (set notation). The integration grid requires
  density values at grid points $Q=\{q_i\}$, that are distinct from $\{p_i\}$. In the `varsub` version of this however, we do not
  know the values of the density that will be required in advance so we have to interpolate on the fly.
  In order to not take an interpolation of an interpolation, the "raw" density values are also loaded on top of the
  standard value of rho. 
  **As a result, we are probably creating the probably initializing the density twice, and perhaps not deleting the allocated array afterwards.**

## 3. varsub in general
  There are many files that are repeated in the varsub version of density modules, and other folders that have not been modified at all from the non-varsub versions
  **This should be fixed so that if a change is made it only has to be made in one place.**
  The easiest way to fix this is to make a symlink to each of the files that hasn't been changed instead of actually copying it over.
  The "right" way to do this is to make a new folder inside `density-modules` called `varsub` and then change the make files accordingly.

## 4. Use the symmetries
Right now the code does not have a way to take advantage of symmetries in the calculation. Many matrix elements are basically 
identical to each other, or factors of `i` times each other, and this is not being used right now. 
Properly using this will significantly improve runtime (only do this for twobody, onebody is super fast).
This (probably) happens easiest by changing the 2Bkernel.f file and pretending like the `Kernel2B` array is smaller than
it actually is, and then applying a mapping function to the actualy `Kernel2B` array
