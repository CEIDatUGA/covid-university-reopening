%% ODE system for Trait-based MBP model
% Made for input into a MATLAB ODE solver like ODE45

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% t - "time", needed for ODE solver but actually unused
% y - output vector name, a 1x4 vector
% p - parameter vector, a 1xN vector of model parameters
%
% Outputs:
% dydt - a 1x6 vector output to ODE45
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function dydt = University_Model(t,y,p)
% Initialize the vector of state variables
% y = [S,L,A,I,R]
dydt = zeros(5,1);

%% Assign parameter names
% In order to more easily reference them in the ODE system equations

% Initially, we will write them in this order:
% b0, b1, mu_H, eta_H, gamma_H, f, K_L, rho_L, mu_V, eta_V, sigma_V,
% sigma_H, beta_H, beta_V

% Create a cell of parameter values
Params = num2cell(p);
% Assign names all at once
[beta,a,bL,bA,gammaI,gammaA,sigma,xi,S0] = Params{:};

%% Assign variable names
S = y(1);
L = y(2);
A = y(3);
I = y(4);
R = y(5);

%% Incidence terms
force_of_infection = beta*(I+bL*L+bA*A);

%% System of DEs
% Susceptible
dydt(1) = -force_of_infection*S;
% Latent
dydt(2) = (1-a)*force_of_infection*S-(sigma+xi)*L;
% Asymptomatic
dydt(3) = a*force_of_infection*S - (gammaA+xi)*A;
% Infected
dydt(4) = sigma*L-gammaI*I;
% Removed/recovered
dydt(5) = gammaI*I+gammaA*A+xi*L+xi*A;



end
