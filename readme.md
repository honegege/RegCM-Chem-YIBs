The RegCM-Chem-YIBs model enables coupling between climate, chemistry and ecology.

Following the steps below, the user can install the model and start using it.

code replace and add
1. Replace the source code under RegCM-Chem-YIBs/Main with the file source code mod_lm_interface.F90, mod_ncio.F90, mod_output.F90 in the code folder
2. Replace the source code under RegCM-Chem-YIBs/Main/chemlib with the file source code in the code/chemlib folder
3. Replace the source code under RegCM-Chem-YIBs/Main/mpplib 
4. Add the yibs directory (the directory contains all the YIBs program source code)in code folder to the RegCM-Chem-YIBs/Main directory 

install
1. ./configure --enable-clm45 CC=icc FC=ifort --with-netcdf=/your netcdf directory
***CC and GCC depend on your compiler
2. make
3. make install

Running
1. Create a run directory
2. ../bin/terrainCLM45 namelist_example.in
3. ../bin/sstCLM45 namelist_example.in
4. ../bin/icbcCLM45 namelist_example.in
5. ../bin/emcre_gridCLM45 namelist_example.in
6. ../bin/interp_emissions namelist_example.in
7. ../bin/chem_icbcCLM45 namelist_example.in
8. ../bin/mksurfdataCLM45 namelist_example.in
9. Add CO2 species to CHBC file
***CO2 species can be added to the CHB(file generated in step 7) using the cdo command (which you need to install first)
10. Create a CO2 source input file
***You can use the cdo command (which you need to install first) in adding CO2 emission source files,which you can make with reference to the sample data
11. Create the yibs input file
***You can refer to the example data yibs_driver.nc to create yibs input files using ncl
12. ../bin/regcmMPICLM45 namelist_example.in 