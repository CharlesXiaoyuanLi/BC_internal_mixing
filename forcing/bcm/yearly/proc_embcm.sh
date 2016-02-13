#!/bin/bash

for ncfile in {1990,1970,1950,1930,1910,1890,1860}
do

    mv ${ncfile}_imbcm.nc ${ncfile}_imbcm_0.05.nc
    mv ${ncfile}_csbcm.nc ${ncfile}_csbcm_0.05.nc
done
