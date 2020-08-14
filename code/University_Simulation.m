% Parameter values
a = 0.45;
bL = 0.3;
bA = 0.5;
gammaI = 1;
gammaA = 1/14;
sigma = 1/4;
Ntot = 49776;
imports = 85;
initial_removed = 456;
S0 = Ntot - imports - initial_removed;
xi = 0;

% Solve for beta with R0=3
R0 = 3;
beta = R0./(S0.*((1-a).*(bL./sigma + 1./gammaI)+a.*bA./gammaA));

% Set testing rate
xi = 0.006;
RL = beta.*S0.*(1-a).*bL./(sigma+xi);
RI = beta.*S0.*(1-a).*sigma./((gammaI).*(sigma+xi));
RA = beta.*S0.*a.*bA./(gammaA+xi);
R0_xi = RL+RI+RA; % R0 with testing
beta_generalized = (1-0.245)*beta;
R0_generalized = R0_xi*beta_generalized/beta; % R0 with testing and generalized interventions

% Parameter vectors
p = [beta,a,bL,bA,gammaI,gammaA,sigma,xi,S0];
p_generalized = p; p_generalized(1) = beta_generalized;

% Compute trajectories
tspan = [0,1000];
ICs = [S0-85, 45, 40,0, 400];
[time,N]= ode45( @(t,N) University_Model(t,N,p),tspan,ICs);
[time_generalized,N_generalized]= ode45( @(t,N) University_Model(t,N,p_generalized),tspan,ICs);
%% Plotting trajectories
% No generalized interventions
figure
title(['Infected individual time series: no generalized intervention (R0=',num2str(R0_xi,3),')'])
hold on
plot(time,N(:,2:4),'LineWidth',3) % Only plot Latent, Asymptomatic, Infected
legend({'L','A','I'})
xlabel('Days')
ylim([0,max(max(N(:,2:4)))])

% with generalized interventions
figure
title(['Infected individual time series: with generalized intervention (R0=',num2str(R0_generalized,3),')'])
hold on
plot(time_generalized,N_generalized(:,2:4),'LineWidth',3) % Only plot Latent, Asymptomatic, Infected
legend({'L','A','I'})
xlabel('Days')
ylim([0,max(max(N_generalized(:,2:4)))])


