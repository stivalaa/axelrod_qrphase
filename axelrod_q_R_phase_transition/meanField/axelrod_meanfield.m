% This function defines a set of (strongly) coupled linear differential
% equations for the mean-field analysis of the Axelrod culture model
%
% Supplementary material for the paper:
% Stivala, A. & Keeler, P. "Another phase transition in the Axelrod model"
% 2016 (submitted to arXiv).

function dPm=axelrod_meanfield(g,t,y)
% This function script defines a set of nonlinear coupled differential
% equations derived in the single-bond mean field analysis by Castellano, 
% Marsili and Vespignani (2000) of the Axelrod model for cultural
% dissemination proposed by Axlerod (1997)
% g is the lattice coordination number, or number of participating agents
yLength=length(y); 
F=yLength-1; %culture vector number
dPm=zeros(yLength,1);
kkValues=(0:F)';
%rho value -- can also be fixed to a constant (Castellano et al 2000).
rho_t=sum(y(2:end).*kkValues(2:end))/F;
%rho_t=0.01;

%create W transition probabilities for Axelrod mean-field analysis.
% Wmm is a (F+1) x (F+1) square matrix
%row number is initial state, column is final state
%eg W_{2,1} corresponds to third row, second column ie Wmm(3,2)
Wnneg=(0:F)/F; %W_{n,n-1}
%Wnzero=(1-Wnneg)*(1-rho_t); % W_{n,n}
Wnpos=(1-Wnneg)*(rho_t); %W_{n,n+1}
%Wmm=diag(Wnzero)+diag(Wnneg(2:end),-1) +  diag(Wnpos(1:end-1),1);

%calculating the first (ie linear) term
for m=2:yLength
    mm=m-1; %m is array index, mm corresponds to the m subscript in the paper
    firstTerm=(mm-1)/F*y(m-1)-mm/F*y(m);
    dPm(m)= firstTerm;
    
end

%calculating second (ie nonlinear) term
kkReducedValues=(1:F)';
coeff=(g-1)*sum(y(2:end).*kkReducedValues)/F;
for m=2:yLength-1
    secondTerm=y(m-1)*Wnpos(m-1)+y(m+1)*Wnneg(m+1)-y(m)*(Wnneg(m)+Wnpos(m));
    dPm(m)=dPm(m)+coeff*secondTerm ;
end

m=yLength;
secondTerm=y(m-1)*Wnpos(m-1)-y(m)*Wnneg(m);
dPm(m)= dPm(m)+coeff*secondTerm ;
dPm(1)=-sum(dPm(2:end));
