% Colored Noise Input
% Function to find the mutual information for a 3-state ring (e.g. ChR2)

function [MI_full, MI_part] = MI_3state_ring_colored(p)

%% Input Parameters

%M = 10000; %Number of input sources.
N = 1;

%p = alpha/(alpha+beta);
gamma = 1e-1; %gamma = alpha + beta
alpha = p*gamma;
beta = gamma - alpha;
sigma = sqrt((2*alpha*beta)/((alpha+beta)^2));

%% Build the transition matrix A

a0 = zeros(3,3);% rates with no input
a1 = zeros(3,3);% rates proportional to input

% rates with no input
a0(2,1)=0; % transition rate from state 1 to state 2
a0(3,2)=50; % transition rate from state 2 to state 3
a0(1,3)=17; % transition rate from state 3 to state 1
a0=a0-diag(sum(a0,1)); % cols sum to zero

% rates proportional to input 
a1(2,1)=5e3; % transition rate from state 1 to state 2
a1(3,2)=0; % transition rate from state 2 to state 3
a1(1,3)=0; % transition rate from state 3 to state 1
a1=a1-diag(sum(a1,1)); % cols sum to zero

% average transition rates (not conditioned on input)
abar=a0+p*a1;
% note cols of abar sum to zero.


%% Compute Steady-State, b, and B
        
% normfactor = abar(3,2)*abar(1,3)+abar(2,1)*abar(1,3)+abar(2,1)*abar(3,2);
% pi1 = abar(3,2)*abar(1,3)/normfactor;
% pi2 = abar(1,3)*abar(2,1)/normfactor;
% pi3 = abar(2,1)*abar(3,2)/normfactor;

[v_A,d_A]=eig(abar);
d_A=diag(d_A);

% Find the steady state from the zero eigenvector
steadyStateDist = v_A(:,find(abs(d_A)<1e-9));
steadyStateDist = steadyStateDist/(sum(steadyStateDist));
if steadyStateDist < 0,
    steadyStateDist = -1*steadyStateDist;
end
pi1 = steadyStateDist(1);
pi2 = steadyStateDist(2);
pi3 = steadyStateDist(3);


% [Ybar1;Ybar2;Ybar3] is the steady-state dist and a right eigenvector of A 
if abs(pi1 + pi2 + pi3 - 1) > 1e-6
    warning('Steady-state distribution does not sum to one')
end

% input vector b captures effects of the input; scale by sqrt(p*(1-p))
b=sqrt(N)*[-pi1*a1(2,1)+pi3*a1(1,3);
    pi1*a1(2,1)-pi2*a1(3,2);
    pi2*a1(3,2)-pi3*a1(1,3)];
    
B = [-sqrt(pi1*abar(2,1)), 0, sqrt(pi3*abar(1,3));
    sqrt(pi1*abar(2,1)), -sqrt(pi2*abar(3,2)), 0;
    0, sqrt(pi2*abar(3,2)), -sqrt(pi3*abar(1,3))];

Ar = [-abar(1,3)-abar(2,1), -abar(1,3); abar(2,1), -abar(3,2)];
Br = [1, 0, 0; 0, 1, 0] * B;
br = [1, 0, 0; 0, 1, 0] * b;

BrBrt = Br*Br';
brbrt = br*br';


%% Partially Observed System
C = [1, 0]';


omega_list = [-fliplr(logspace(-9, 8, 1000)), logspace(-9, 8, 1000)]; %creates a range of 2000 frequencies to numerically integrate over
SE_part = nan(1, length(omega_list)); %integrand in the spectral efficiency integral
SE_full = nan(1, length(omega_list)); %integrand in the spectral efficiency integral
powerSpectrum_part_uncnd = nan(1, length(omega_list));
powerSpectrum_part_cnd = nan(1, length(omega_list));
powerSpectrum_full_uncnd = nan(1, length(omega_list));
powerSpectrum_full_cnd = nan(1, length(omega_list));
for k = 1:length(omega_list)
    omega = omega_list(k);
    
    powerSpectrum_part_uncnd(k) = det(C'*(((Ar-eye(2)*sqrt(-1)*omega)\(sigma^2/(gamma^2+omega^2)*brbrt + BrBrt))/(Ar'+eye(2)*sqrt(-1)*omega))*C/(2*pi)); %power spectrum unconditioned on input
    powerSpectrum_part_cnd(k) = det(C'*(((Ar-eye(2)*sqrt(-1)*omega)\(BrBrt))/(Ar'+eye(2)*sqrt(-1)*omega))*C/(2*pi)); %power spectrum conditioned on input
    SE_part(k) = log(powerSpectrum_part_uncnd(k)) - log(powerSpectrum_part_cnd(k)); %integrand in the spectral efficiency integral

    powerSpectrum_full_uncnd(k) = det((((Ar-eye(2)*sqrt(-1)*omega)\(sigma^2/(gamma^2+omega^2)*brbrt + BrBrt))/(Ar'+eye(2)*sqrt(-1)*omega))/(2*pi)); %power spectrum unconditioned on input
    powerSpectrum_full_cnd(k) = det((((Ar-eye(2)*sqrt(-1)*omega)\(BrBrt))/(Ar'+eye(2)*sqrt(-1)*omega))/(2*pi)); %power spectrum conditioned on input
    SE_full(k) = log(powerSpectrum_full_uncnd(k)) - log(powerSpectrum_full_cnd(k)); %integrand in the spectral efficiency integral

end

MI_part = 1/2*real(trapz(omega_list, SE_part));
MI_full = 1/2*real(trapz(omega_list, SE_full));






end
