function n_a = axelrod_solve(g, F, q)
%axelrod_solve Solves a set of coupled ODEs using ode45.
%   Solves a set of coupled ODEs using ode45. The solution corresponds to 
%   the mean-field analysis of the Axelrod culture dissemination model
%     g is lattice coordination number of number of participating agents
%     F is length of culture vector
%     q is number of traits
%
% Supplementary material for the paper:
% Stivala, A. & Keeler, P. "Another phase transition in the Axelrod model"
% 2016 (submitted to arXiv).

mmValues=(0:F); %calculate m index values


tspan=[0,10^3]; %time span for the ODE solution
%initial values  for numerical method
rho0=1/q; %probability of two cultural uniform elements coinciding
%rho0=besseli(0,2*q)*exp(-2*q); %if elements are Poisson distributed

%binomial coefficent ie F choose m
F_choose_m=gamma(F+1)./gamma(mmValues+1)./gamma(F-mmValues+1);
%intial values -- binomial variable
PmIntial=F_choose_m.*rho0.^(mmValues).*(1-rho0).^(F-mmValues);
%PmIntial=fliplr(PmIntial);
%PmIntial= ones(1,F+1)/(F+1);

%running ode 45, which takes function, time span and initial conditions
[t,PmMatrix]=ode45(@(u,v)axelrod_meanfield(g,u,v),tspan,PmIntial);
%each coloumn corresponds to a set of Pm values
n_a=sum(PmMatrix(:,2:end-1),2);
