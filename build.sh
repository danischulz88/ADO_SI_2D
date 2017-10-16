#!/bin/sh

gfortran adosibidimensional.f95 -L . -llapack+blas_optim -o adosi2d
./adosi2d
