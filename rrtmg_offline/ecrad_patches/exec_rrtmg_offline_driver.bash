#!/bin/bash
data_dir=../data
exec_dir=../build
mkdir -p rrmtg
cd rrmtg
cp -f ${data_dir}/test_pressure_level.nc .
cp -f ${data_dir}/test_single_level.nc .
cp -f ${data_dir}/RRTMG_* .
cp -f ${exec_dir}/rrtmg_offline_driver .
./rrtmg_offline_driver test_pressure_level.nc test_single_level.nc 1 1 1 2 saida.csv
