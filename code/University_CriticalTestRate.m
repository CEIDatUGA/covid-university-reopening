%% Parameter values
% assuming 15% variation from baseline
var = 0.15;
var_vec = [1,1-var,1+var];
a = [0.45 0.4 0.6];%0.45.*var_vec;
bL = 0.3.*var_vec;
bA = 0.5.*var_vec;
gammaI = 1.*var_vec;
gammaA = 1/14.*var_vec;
sigma = [1/4 1/5 1/3];%1/4.*var_vec;
Ntot = 49776;
S0 = Ntot-[907+456,1606+456,208+456];
p_ranges = [a;bL;bA;gammaI;gammaA;sigma;S0];
Parameter_names = {'\beta','a','b_L','b_A','\gamma_I','\gamma_A','\sigma','S_0'};
% xi = 0;

% Solve for beta with fixed R0=3
R0 = 3;
beta = R0./(S0.*((1-a).*(bL./sigma + 1./gammaI)+a.*bA./gammaA));
beta = beta(1).*[1 1 1]; % don't need a range of beta values, just baseline
beta_generalized = (1-0.245).*beta;
p_ranges = [beta_generalized; p_ranges];

%% Set up sensitivity analysis design
Num_samples = 200;
Num_params = size(p_ranges,1);
Samples = lhsdesign(Num_samples,Num_params);
Parameter_samples = Samples'.*range(p_ranges(:,2:3)')' + p_ranges(:,2);
for i = 1:Num_samples
    p = Parameter_samples(:,i);
    xi_c(i) =  get_crit_test_rate(p);
end

% Compute partial rank correlation coefficients
% Rank parametervectors
x=Parameter_samples';
[a, b]=size(x);
for j=1:b
    [~,i]=sort(x(:,j));
    r(i,j)=[1:a]';
end
Xranked = r;
% Rank outputs
y = xi_c;
n=length(y);
[s,i]=sort(y);                                                
Yranked(i,1)=[1:n];

% Compute PRCCs
[xi_c_prcc, xi_c_pvalues, ~]=PRCC(Parameter_samples',xi_c',[],Parameter_names,0.05);


%% Collect and save tables
% Variable_names = {'a','bL','bA','gammaI','gammaA','sigma','S0','xi_c'};
% CriticalTesting_samples = [Parameter_samples',xi_c'];
% CriticalTesting_Ranked_samples = [Xranked, Yranked];
% CriticalTesting_PRCC = [xi_c_prcc;xi_c_pvalues];
% 
% tCriticalTesting_samples = array2table(CriticalTesting_samples,'VariableNames',Variable_names);
% tCriticalTesting_Ranked_samples = array2table(CriticalTesting_Ranked_samples,'VariableNames',Variable_names);
% tCriticalTesting_PRCC = array2table(CriticalTesting_PRCC,'VariableNames',Parameter_names);
% 
% writetable(tCriticalTesting_samples,'CriticalTesting_LHSsamples.csv')
% writetable(tCriticalTesting_Ranked_samples,'CriticalTesting_Rankedsamples.csv')
% writetable(tCriticalTesting_PRCC,'CriticalTesting_PRCCvalues.csv')

%% Plots
% Scatter plot of parameter and xi_c values
for i = 1:8
    subplot(3,3,i)
    plot(Parameter_samples(i,:),xi_c,'.')
    xlabel(Parameter_names{i})
    ylabel('\xi_c')
end

% Scatter plot of RANKED parameter and xi_c values with PRCC measure
PRCC_Plot(Parameter_samples',xi_c,Parameter_names,'\xi_c')

%% Subfunctions

% Compute critical testing rate
function xi_c = get_crit_test_rate(p)
% Assign parameter names for ease of reading
Params = num2cell(p);
[beta,a,bL,bA,gammaI,gammaA,sigma,S0] = Params{:};
% Compute reproduction numbers and beta
RL0 = beta.*S0.*(1-a).*bL./sigma;
RI0 = beta.*S0.*(1-a)./gammaI;
RA0 = beta.*S0.*a.*bA./gammaA;
R0 = RL0+RI0+RA0;
% beta = R0./(S0.*((1-a).*(bL./sigma + 1./gammaI)+a.*bA./gammaA));
% RA0 = beta.*S0.*a.*bA./(gammaA);
% Determine critical testing rate (as a function of other parameters)
xi_c = -0.5.*(sigma.*(1-R0+RA0)+gammaA.*(1-RA0))+0.5.*...
    sqrt((sigma.*(1-R0+RA0)+gammaA.*(1-RA0)).^2 + 4.*sigma.*gammaA.*(R0-1));
end
