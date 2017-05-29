#!/usr/bin/env bash
gfortran -fcheck=all *.f -o ../bin/ctest 2>c_err.log
