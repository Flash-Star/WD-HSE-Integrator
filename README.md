# WD-HSE-Integrator

HSE Integrator for WD Models using F. Timmes' Helmholtz EOS and user-supplied T, X profile


## Original version by Dean Townsley, 2009/3/30

This code builds a WD in hydrostatic equilibrium on an evenly spaced radial
grid.  The finite-difference form of the HSE equations are taken from Zingale
etal 2002, ApJS, 143, 539, which account in some sense for the PPM
reconstruction of the fluid quantities.

This uses Frank Timmes' standalone EOS from cococubed and is therefore clean
of DOE-paid-for code.  However it uses GSL, and therefore technically must be
GPL.  Currently it is unreleased, and therefore a bit vague on this point.

The main program is SNproj, which builds a WD profile, density and temperature
as a function of radius given the central density and temperature (as command
line arguments).  It reports the resulting mass at the end.  'make' should
build this with the existing makefile.  You will need to have GSL installed,
but this usually is available with your linux distribution.

If your target is an input file for the SNIa_ddt setup in flash, you will
need to add the carbon and neon abundance columns manually and put the number
of lines at the top.  Also only one comment line is allowed.  There are a few
example files for this.  (the cc_* files)

## Revisions by Donald Willcox, 2016:

SNproj looks for a file in the current directory called
`ref_profile.dat` at the same uniform grid resolution as the model
SNproj integrates. During the HSE integration, SNproj uses the
temperature and species mass fractions from `ref_profile.dat` as inputs
to the EOS call. This helps make sure the energetics of the model are
consistent with a uniform grid model derived from another calculation,
e.g. 1-D stellar evolution codes.

SNproj is set up to read (from `ref_profile.dat`) and use the following
species: C12, O16, Ne20, Ne22. Only C12, Ne20, Ne22 are written in the
output, where it is assumed 1 = X(C12) + X(O16) + X(Ne20) + X(Ne22).

The EOS is Frank Timmes' Helmholtz EOS, and the table is the version
provided on [Frank's website](http://cococubed.asu.edu) as of December
8, 2014. This corresponds to the EOS used for the WD models and Flash
simulations in Willcox, et al 2016.

Look at the file compute_4km_flash_profile.sh for an example of how to
use SNproj.

The GSL routines in SNproj.c have been removed in favor of manually doing a
combination of Newton iteration and brute force searching to get the
HSE density in each zone.

NOTICE RE. COPYRIGHT STATUS:

Because this code uses the GSL, we must use the GPL (included in the
LICENSE file), and this applies to the following 6 source files:

- compute_4km_flash_profile.sh
- Teos.h
- cgs.h
- SNproj.c
- composition.h
- Teos_f.c

However, the GPL does NOT apply to the following, as they are from
Frank Timmes' Helmholtz EOS distribution:

- helm_eos.f
- helm_table.dat
- vector_eos.dek