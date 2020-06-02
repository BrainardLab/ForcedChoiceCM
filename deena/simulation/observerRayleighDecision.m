function [p1_up, t_up] = observerRayleighDecision(observer, primarySpd, testSpd) 
% Simulated observer decision making for a Rayleigh match experiment 
%
% Syntax:
%   observerRayleighDecision(T_observer, primarySpd, testSpd); 
%
% Description:
%    Takes in a simulated observer's cone fundamentals and a pair of
%    spds (primary and test). 
%
%    Computes the opponent contrast of the two spectra based on the
%    observer's cone fundamentals. Based on the contrast values of
%    individual channels, determines whether the test light should be made
%    brighter or dimmer and whether the mixing light should be shifted
%    towards the red primary or the green primary. 
%
%   
% Inputs:
%     observer          -Struct with observer parameters. Must contain the
%                        following fields: colorDiffParams, T_cones. 
%     primarySpd        -201x1 vector representation of the primary spd
%     testSpd           -201x1 vector represetation of the test spd  
%
% Outputs:
%     p1_up             -Logical indicating whether  primary ratio should 
%                        be shifted towards p1 (red). 
%     t_up              -Logical indicating whether test intensity should  
%                        be increased. 
%
% Optional key-value pairs:
%    None

% History:
%   06/02/20  dce       Wrote initial code

% Cone responses for the given spectra 
tLMS = observer.T_cones * testSpd; 
pLMS = observer.T_cones * primarySpd; 

% Opponent contrasts for the given spectra 
opponentContrast = LMSToOpponentContrast(observer.colorDiffParams,...
    tLMS, pLMS);

% Check luminance 
if opponentContrast(1) > 0 % Primary is brighter 
    t_up = true;
else 
    t_up = false; 
end 

% Check R/G ratio
if opponentContrast(2) < 0 % Too green/too close to p2
    p1_up = true; 
else 
    p1_up = false; 
end 
end 