function p = getParameters()
%GETPARAMETERS Basic parameters

% Define the population
load('./data/popdata20.mat');
p.N = sum(pop,2);
p.Na = pop;
p.Nltlas = length(p.N);
 
p.sigma = [1 0.4 0.25]; % Reduction in susceptibility by vaccine
p.d = [1 0.4 0.17]; % Probability of symptoms by vaccine
p.tt = [1 0.55 0.55]; % Reduction in transmission from vaccine
p.gamma = 0.52; % Recovery rate
p.epsilon = 1/5.5; % Rate of moving through exposed compartments
p.delta = 0.7; % Proportion of day in daytime

ageprop = p.Na./sum(p.Na,2);

end
