
clear; clc;
warning off;

display('Compling Main Mex Functions...')
cd 'mexFunctions';
compileMex;
cd ..

display('Compling ANN Mex Functions...')
cd 'funcs/ann'
make;
cd ../..