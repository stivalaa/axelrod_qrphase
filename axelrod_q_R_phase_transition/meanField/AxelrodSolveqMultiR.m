% Solves a set of coupled ODEs using ode45. The solution corresponds to the
% mean-field analysis of the Axelrod culture dissemination model.
% This version solves for multiple values of q and plots phase diagram
% with q on the x axis and n_a(t) number of active links on y axis.
%
% Supplementary material for the paper:
% Stivala, A. & Keeler, P. "Another phase transition in the Axelrod model"
% 2016 (submitted to arXiv).

clear all; close all; clc;

F=5; %Length of Axelrod cultral vector
Rv=1:4; % vector of R = von Neumann radius
qv = 1:500; % vector of q = maximum value for entries of culture vector
n_a_matrix = zeros(length(Rv), length(qv));
tic
for R = Rv
    R
    for q  = qv
        g = 2*R*(R+1)+1; % von Neumann neigborhood radius R (+1 for focal agent)
        n_a = axelrod_solve(g, F, q);
        n_a_matrix(R, q) = n_a(end);
    end
    toc
end
%%plotting numerical solution
linestyles = {'-g', '--r', '-.b', ':m', '--c'};
legendlabels = {'R = 1', 'R = 2', 'R = 3', 'R = 4'};
figure;
hold on;
for R = Rv
    plot(qv, n_a_matrix(R,:), linestyles{R}, 'LineWidth', 2);
end
axis([0 500 0.0 1.0])
xlabel('q'); ylabel('n_a(t)');
legend(legendlabels);
legend('boxoff');
box on;
print('meanfield_q_multiradius', '-depsc');