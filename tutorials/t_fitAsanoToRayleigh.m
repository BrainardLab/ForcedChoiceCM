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

