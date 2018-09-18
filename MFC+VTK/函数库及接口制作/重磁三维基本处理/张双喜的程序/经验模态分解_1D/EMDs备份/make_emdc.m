%MAKE_EMDC  Compiles the C codes for Empirical Mode Decomposition
%
% Note: The compilation can fail on some systems (e.g. MacOS) if Matlab cannot find the C compiler.
% In this case, you should either install a C compiler or check Matlab configuration.
% The configuration files for compilation are mexopts.sh (Unix / MAC OS) and mexopts.bat (Windows)
% use "mex -setup" to choose a configuration file (Unix) or select a compiler (Windows).

function make_emdc
clc;
oldpwd = pwd;
path = fileparts(which('make_emdc'));
cd(path)

mex src/emdc.c
mex src/emdc_fix.c
mex src/cemdc.c
mex src/cemdc_fix.c
mex src/cemdc2.c
mex src/cemdc2_fix.c

cd(oldpwd)