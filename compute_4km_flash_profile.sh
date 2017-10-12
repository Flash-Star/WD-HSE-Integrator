# Calling SNproj with a uniform grid profile (ref_profile.dat)
# from which to take temperature and mass fractions.
ln -s profile75_Dr-400k_qic_renormed.dat ref_profile.dat
./SNproj > 400k_flash.dat
