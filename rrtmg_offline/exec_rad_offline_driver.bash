#!/bin/bash -x
mkdir -p rad
rad_exec=./rad
data_dir=../data
exec_dir=../build
cd ${rad_exec}
cp -f ${exec_dir}/rad_offline_driver .
cp -f ${data_dir}/RRTMG_* .
cp -f ${data_dir}/era5_pressure_level.nc   era5_pressure_level.nc
cp -f ${data_dir}/era5_single_level.nc     era5_single_level.nc
cp -f ../config.nam.sample ./config.nam
cp -f ../namelist.offline.sample ./namelist.offline
export OMP_NUM_THREADS=1
./rad_offline_driver namelist.offline
