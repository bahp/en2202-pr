% Author : Bernard Hernandez and Fernando Iglesias
% Date : 20/09/2012

% ------------------------------------------------------------------------
%                    LAB3 - FORWARD IMPLEMENTATION
% ------------------------------------------------------------------------
clear; clc; close all;

% Create a Markov Chain.
mc = MarkovChain([1; 0], [0.9 0.1 0; 0 0.9 0.1]);
f_state1 = GaussD('Mean', 0, 'StDev', 1); % N(mu=0,sigma=1)
f_state2 = GaussD('Mean', 3, 'StDev', 2); % N(mu=3,sigma=2)
obs_seq = [-0.2, 2.6, 1.3]; % Observed finite-duration sequence. 

% Create observation sequences as gaussians.
[pX logS] = prob([f_state1, f_state2],obs_seq);
% Other way:
% g = [f_state1; f_state2]; pX = g.prob(obs_seq)

% Test forward.
disp('Forward result:')
[alfaHat c] = forward(mc,pX)
% Other way:
% [alfaHat c] = mc.forward(pX)




% Create the finite-duration test HMM.
h = HMM(mc, [f_state1; f_state2]);

% Test logprob.
disp('logprob result:');
logP = logprob(h,obs_seq)