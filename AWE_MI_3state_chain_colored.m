% Function to find the mutual information in a three state chain model with
% colored noise. Function takes p as an input and fixes gamma = alpha +
% beta

function [MI_full, MI_part, fsv_dot] = AWE_MI_3state_chain_colored(p)


%% Input Parameters

%M = 1; %Number of input sources. NO LONGER NECESSARY WITH OUR SCALING
N = 1; %Number of receptors

%p = alpha/(alpha+beta);
gamma = 1e-2; %gamma = alpha + beta
alpha = p*gamma;
beta = gamma - alpha;
sigma = sqrt((2*alpha*beta)/((alpha+beta)^2));


%% Build the transition matrix A

a0 = zeros(3,3);% rates with no input
a1 = zeros(3,3);% rates proportional to input

% rates with no input
a0(2,1)=1; % transition rate from state 1 to state 2
a0(1,2)=1; % transition rate from state 2 to state 1
a0(3,2)=1; % transition rate from state 2 to state 3
a0(2,3)=1; % transition rate from state 3 to state 2
a0=a0-diag(sum(a0,1)); % cols sum to zero

% rates proportional to input 
a1(2,1)=0; % transition rate from state 1 to state 2
a1(1,2)=0; % transition rate from state 2 to state 1
a1(3,2)=1; % transition rate from state 2 to state 3
a1(2,3)=0; % transition rate from state 3 to state 2
a1=a1-diag(sum(a1,1)); % cols sum to zero

% average transition rates (not conditioned on input)
abar=a0+p*a1;
% note cols of abar sum to zero.

%% Compute Steady-State, b, and B

normfactor = abar(2,3)*abar(1,2) + abar(2,3)*abar(2,1) + abar(2,1)*abar(3,2);
pi1=abar(2,3)*abar(1,2)/normfactor;
pi2=abar(2,3)*abar(2,1)/normfactor;
pi3=abar(2,1)*abar(3,2)/normfactor;

J12 = abar(2,3)*abar(1,2)*abar(2,1)/normfactor;
J23 = abar(2,1)*abar(3,2)*abar(2,3)/normfactor;

% note [pi1;pi2;pi3] is the ss dist and a right eigenvector of A' 
if abs(pi1 + pi2 + pi3 - 1) > 1e-6
    warning('Steady-state distribution does not sum to one')
end


% input vector b captures effects of the input; scale by sqrt(p*(1-p))
% %b=sqrt(p*(1-p))*[-pi1*a1(2,1)+pi2*a1(1,2);
%     pi1*a1(2,1)-pi2*(a1(1,2)+a1(3,2))+pi3*a1(2,3);
%     pi2*a1(3,2)-pi3*a1(2,3)];
b=sqrt(N)*[-pi1*a1(2,1)+pi2*a1(1,2);
    pi1*a1(2,1)-pi2*(a1(1,2)+a1(3,2))+pi3*a1(2,3);
    pi2*a1(3,2)-pi3*a1(2,3)];


% B is an m by n matrix where m is the number of states and n is the number
% of edges.  Each column corresponds to an edge
B = [-sqrt(pi1*abar(2,1)), sqrt(pi2*abar(1,2)), 0, 0;
    sqrt(pi1*abar(2,1)), -sqrt(pi2*abar(1,2)), -sqrt(pi2*abar(3,2)), sqrt(pi3*abar(2,3));
    0, 0, sqrt(pi2*abar(3,2)), -sqrt(pi3*abar(2,3))];


%% Calculating the Mutual Information for Fully Observed (Third State obs)
C = [0;1];

Ar = [-abar(2,1) - abar(1,2), -abar(1,2); -abar(3,2), -abar(3,2) - abar(2,3)];
br = sqrt(N)*[-pi1*a1(2,1)+pi2*a1(1,2); pi2*a1(3,2)-pi3*a1(2,3)];
Br = [sqrt(2*J12), 0; 0, sqrt(2*J23)];

omega_list = [-fliplr(logspace(-9, 8, 1000)), logspace(-9, 8, 1000)]; %creates a range of 2000 frequencies to numerically integrate over

SE_part = nan(1, length(omega_list)); %integrand in the spectral efficiency integral
powerSpectrum_part_uncnd = nan(1, length(omega_list));
powerSpectrum_part_cnd = nan(1, length(omega_list));

SE_full = nan(1, length(omega_list)); %integrand in the spectral efficiency integral
powerSpectrum_full_uncnd = nan(1, length(omega_list));
powerSpectrum_full_cnd = nan(1, length(omega_list));

%fsv = nan(3,length(omega_list);
fsv_dot = nan(2,length(omega_list));

for k = 1:length(omega_list)
    omega = omega_list(k);
    
    powerSpectrum_part_uncnd(k) = det(C'*(((Ar-eye(2)*sqrt(-1)*omega)\(sigma^2/(gamma^2+omega^2)*(br*br') + Br*Br'))/(Ar'+eye(2)*sqrt(-1)*omega))*C/(2*pi)); %power spectrum unconditioned on input
    powerSpectrum_part_cnd(k) = det(C'*(((Ar-eye(2)*sqrt(-1)*omega)\(Br*Br'))/(Ar'+eye(2)*sqrt(-1)*omega))*C/(2*pi)); %power spectrum conditioned on input
        
    SE_part(k) = log(powerSpectrum_part_uncnd(k)) - log(powerSpectrum_part_cnd(k)); %integrand in the spectral efficiency integral
    
    
%     powerSpectrum_full_uncnd(k) = det((((Ar-eye(2)*sqrt(-1)*omega)\(sigma^2/(gamma^2+omega^2)*(br*br') + Br*Br'))/(Ar'+eye(2)*sqrt(-1)*omega))/(2*pi)); %power spectrum unconditioned on input
%     powerSpectrum_full_cnd(k) = det((((Ar-eye(2)*sqrt(-1)*omega)\(Br*Br'))/(Ar'+eye(2)*sqrt(-1)*omega))/(2*pi)); %power spectrum conditioned on input
    powerSpectrum_full_uncnd(k) = det((sigma^2/(gamma^2+omega^2)*(br*br') + Br*Br')); %power spectrum unconditioned on input
    powerSpectrum_full_cnd(k) = det((Br*Br')); %power spectrum conditioned on input
    
    fsv_dot(1,k) = sum((inv(Ar + eye(2)*sqrt(-1)*omega)*br) .* C);
    fsv_dot(2,k) = sqrt(sum((inv(Ar + eye(2)*sqrt(-1)*omega)*br).*(inv(Ar - eye(2)*sqrt(-1)*omega)*br)));
    fsv_dot(2,k) = fsv_dot(2,k) * sqrt(sum(C.*C));

SE_full(k) = log(powerSpectrum_full_uncnd(k)) - log(powerSpectrum_full_cnd(k)); %integrand in the spectral efficiency integral
    
end

MI_part = 1/2*real(trapz(omega_list, SE_part));
MI_full = 1/2*real(trapz(omega_list, SE_full));


%% Calculating the Mutual Information for Fully Observed

% 2 to 3 edge sensitive
%MI_full = pi * (-gamma + sqrt(gamma^2 + b(3)^2*sigma^2/(2*J23)));

% 1 to 2 edge sensitive
%MI_full = pi * (-gamma + sqrt(gamma^2 + b(1)^2*sigma^2/(2*J12)));



end
