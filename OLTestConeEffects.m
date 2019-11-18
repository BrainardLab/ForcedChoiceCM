lambdaMaxes = [558.9 530.3 420.7]';     % Normal trichromat
dphotopigments = [0 0 0]'; 
T_cones = findConeFundamentals(lambdaMaxes, dphotopigments);
load OLSampleMatches11_14_DB.mat;
DBMatches = matches; 
load OLSampleMatches11_14_DE.mat; 
DEMatches = matches; 