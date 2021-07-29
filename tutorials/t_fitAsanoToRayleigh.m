% t_FitAsanoToRayleigh
%
% Description:
%    Demonstration/test tutorial to check/develop fitting methods to
%    Rayleigh match data.
%

% History:
%    07/29/21  dhb  Started.

%% Clear
clear; close all;

%% Load in the sample data
theDir = fileparts(mfilename('fullpath'));
theData = load(fullfile(theDir,'sampleData','MELC_0004_spds.mat'));

%% Structure data 
S = [380 2 201];
wls = SToWls(S);
[nWls,nReps,nRefs] = size(theData.measPrimarySpdsByWl);
refSpds = reshape(theData.measRefSpdsByWl,S(3),nReps*nRefs);
primarySpds = reshape(theData.measPrimarySpdsByWl,S(3),nReps*nRefs);

%% Call Deena fit routine
% findObserverParameters(testSpds,primarySpds,varargin)
findObserverParameters(refSpds,primarySpds);