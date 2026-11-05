#!/bin/bash

cd data
ncgen -o test_pressure_level.nc test_pressure_level.cdl
ncgen -o test_single_level.nc   test_single_level.cdl
ls -la *.nc
cd ..
