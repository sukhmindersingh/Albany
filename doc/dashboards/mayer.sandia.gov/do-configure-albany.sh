#!/bin/bash

cmake \
 -D ALBANY_TRILINOS_DIR:PATH=/home/ikalash/Trilinos/build/install \
 -D ENABLE_LCM:BOOL=ON \
 -D ENABLE_MOR:BOOL=ON \
 -D ENABLE_GOAL:BOOL=OFF \
 -D ENABLE_LANDICE:BOOL=ON \
 -D ENABLE_HYDRIDE:BOOL=ON \
 -D ENABLE_AMP:BOOL=OFF \
 -D ENABLE_ATO:BOOL=ON \
 -D ENABLE_SCOREC:BOOL=OFF \
 -D ENABLE_QCAD:BOOL=ON \
 -D ENABLE_SG:BOOL=OFF \
 -D ENABLE_ENSEMBLE:BOOL=OFF \
 -D ENABLE_ASCR:BOOL=OFF \
 -D ENABLE_AERAS:BOOL=ON \
 -D ENABLE_64BIT_INT:BOOL=OFF \
 -D ENABLE_LAME:BOOL=OFF \
 -D ENABLE_DEMO_PDES:BOOL=ON \
 -D ENABLE_KOKKOS_UNDER_DEVELOPMENT:BOOL=ON \
 -D ALBANY_CTEST_TIMEOUT=400 \
 -D ENABLE_CHECK_FPE:BOOL=OFF \
 ..\
