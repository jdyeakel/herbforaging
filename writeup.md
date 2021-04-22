# Simulation description

The simulation works by drawing travel distances for the forager moving from one resource to the next... m and alpha play into the means and variances of travel distances.

1) First we enter the environmental parameters:
        1) rho = how we tune resource availability
        2) m = mean encounter rates of resources (proportional to resource density)
        3) alpha = dispersion parameter; high value means even dispersion; low value means patchy.
        4) Mean travel distance = m
        5) SD travel distance = m/alpha (check)

The resource parameters are expected to change with body size... for instance, grassland resources look very patchy to a rodent, but very even to an elephant... we capture this with the variability-scaling parameter zeta. This type of scaling is given by zeta = 1. If zeta = 2, the patchiness/variability of resource availability for smaller organisms matches those of larger organisms.

2) We enter the herbivore parameters
body size
tooth type
gut type

3) Within-day calculations: Calculate a distribution for daily return for the gut based on a) chewing rates and travel between foods during foraging period (6 hours)

4) Across-day calculations: Draw from the within-day gut return distribution and track the energetic state of the consumer's gut and fat.
gut(t) = gut(t-1) + deltagut, where deltagut = food consumed - food passed
fat(t) = fat(t-1) + deltafat, where deltafat = conversionrate*food passed - fat used for metabolic processes

All of the rates are allometrically determined, so we should be able to capture general energetic constraints across body sizes for different types of environments, where environment is characterized by energy availability, patchiness, and the scaling of patchiness with body size.


# Notes for Justin
- check within-day foraging for weird clumping artifacts