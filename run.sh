#!/bin/bash  

#gfortran main_cplx.f90
gfortran main.f90
./a.out
python plots.py
rm fort.2
