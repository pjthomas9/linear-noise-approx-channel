% Colored Noise Input
% Function to find the mutual information for ACh

function [MI_full, MI_part] = MI_ACh_colored(p)

%% Input Parameters

%M = 10000; %Number of input sources.

%p = alpha/(alpha+beta);
gamma = 1e-1; %gamma = alpha + beta
alpha = p*gamma;
beta = gamma - alpha;
sigma = sqrt((2*alpha*beta)/((alpha+beta)^2));

%% Build the transition matrix A

a0 = zeros(5,5);% rates with no input
a1 = zeros(5,5);% rates proportional to input

% rates with no input
a0(2,1) = 0; % transition rate from state 1 to state 2
a0(1,2) = 0.66; % transition rate from state 2 to state 1
a0(4,1) = 3e3; % transition rate from state 1 to state 4
a0(1,4) = 15; % transition rate from state 4 to state 1
a0(3,2) = 5e2; % transition rate from state 2 to state 3
a0(2,3) = 1.5e4; % transition rate from state 3 to state 2
a0(4,3) = 4e3; % transition rate from state 3 to state 4
a0(3,4) = 0; % transition rate from state 4 to state 3
a0(5,4) = 2e3; % transition rate from state 4 to state 5
a0(4,5) = 0; % transition rate from state 5 to state 4
a0 = a0-diag(sum(a0,1)); % cols sum to zero

% rates proportional to input 
a1(2,1) = 5e8*1e-6; % transition rate from state 1 to state 2
a1(1,2) = 0; % transition rate from state 2 to state 1
a1(4,1) = 0; % transition rate from state 1 to state 4
a1(1,4) = 0; % transition rate from state 4 to state 1
a1(3,2) = 0; % transition rate from state 2 to state 3
a1(2,3) = 0; % transition rate from state 3 to state 2
a1(4,3) = 0; % transition rate from state 3 to state 4
a1(3,4) = 5e8*1e-6; % transition rate from state 4 to state 3
a1(5,4) = 0; % transition rate from state 4 to state 5
a1(4,5) = 1e8*1e-6; % transition rate from state 5 to state 4
a1 = a1-diag(sum(a1,1)); % cols sum to zero

% average transition rates (not conditioned on input)
abar = a0 + p*a1;
% note cols of abar sum to zero.


%% Compute Steady-State, b, and B
        
% normfactor = abar(3,2)*abar(1,3)+abar(2,1)*abar(1,3)+abar(2,1)*abar(3,2);
% pi1 = abar(3,2)*abar(1,3)/normfactor;
% pi2 = abar(1,3)*abar(2,1)/normfactor;
% pi3 = abar(2,1)*abar(3,2)/normfactor;
% pi4 = 0;
% pi5 = 0;

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
pi4 = steadyStateDist(4);
pi5 = steadyStateDist(5);


% [Ybar1;Ybar2;Ybar3] is the steady-state dist and a right eigenvector of A 
if abs(pi1 + pi2 + pi3 + pi4 + pi5 - 1) > 1e-6
    warning('Steady-state distribution does not sum to one')
end

% input vector b captures effects of the input; scale by sqrt(p*(1-p))
b=[-pi1*(a1(2,1)+a1(4,1))+pi2*a1(1,2)+pi4*a1(1, 4);
    pi1*a1(2,1)-pi2*(a1(1,2)+a1(3,2))+pi3*a1(2,3);
    pi2*a1(3,2)-pi3*(a1(2,3)+a1(4,3))+pi4*a1(3,4);
    pi3*a1(4,3)-pi4*(a1(1,4)+a1(3,4)+a1(5,4))+pi5*a1(4,5);
    pi4*a1(5,4)-pi5*a1(4,5)];

B = [-sqrt(pi1*abar(2,1)), sqrt(pi2*abar(1,2)), 0, 0, 0, 0, -sqrt(pi1*abar(4,1)), sqrt(pi4*abar(1,4)), 0, 0;
    sqrt(pi1*abar(2,1)), -sqrt(pi2*abar(1,2)), -sqrt(pi2*abar(3,2)), sqrt(pi3*abar(2,3)), 0, 0, 0, 0, 0, 0;
    0, 0, sqrt(pi2*abar(3,2)), -sqrt(pi3*abar(2,3)), -sqrt(pi3*abar(4,3)), sqrt(pi4*abar(3,4)), 0, 0, 0, 0;
    0, 0, 0, 0, sqrt(pi3*abar(4,3)), -sqrt(pi4*abar(3,4)), sqrt(pi1*abar(4,1)), -sqrt(pi4*abar(1,4)), -sqrt(pi4*abar(5,4)), sqrt(pi5*abar(4,5));
    0, 0, 0, 0, 0, 0, 0, 0, sqrt(pi4*abar(5, 4)), -sqrt(pi5*abar(4,5))];

Ar = abar(1:4, 1:4) - [0,0,0,0;0,0,0,0;0,0,0,0; abar(4,5) abar(4,5), abar(4,5), abar(4,5)];
Br = [1, 0, 0, 0, 0; 0, 1, 0, 0, 0; 0, 0, 1, 0, 0; 0, 0, 0, 1, 0] * B;
br = [1, 0, 0, 0, 0; 0, 1, 0, 0, 0; 0, 0, 1, 0, 0; 0, 0, 0, 1, 0] * b;

BrBrt = Br*Br';
brbrt = br*br';


%% Partially Observed System
C = [1, 1, 0, 0]';


omega_list = [-fliplr(logspace(-9, 8, 1000)), logspace(-9, 8, 1000)]; %creates a range of 2000 frequencies to numerically integrate over
SE_part = nan(1, length(omega_list)); %integrand in the spectral efficiency integral
SE_full = nan(1, length(omega_list)); %integrand in the spectral efficiency integral
powerSpectrum_part_uncnd = nan(1, length(omega_list));
powerSpectrum_part_cnd = nan(1, length(omega_list));
powerSpectrum_full_uncnd = nan(1, length(omega_list));
powerSpectrum_full_cnd = nan(1, length(omega_list));
for k = 1:length(omega_list)
    omega = omega_list(k);
    
    powerSpectrum_part_uncnd(k) = det(C'*(((Ar-eye(4)*sqrt(-1)*omega)\(sigma^2/(gamma^2+omega^2)*brbrt + BrBrt))/(Ar'+eye(4)*sqrt(-1)*omega))*C/(2*pi)); %power spectrum unconditioned on input
    powerSpectrum_part_cnd(k) = det(C'*(((Ar-eye(4)*sqrt(-1)*omega)\(BrBrt))/(Ar'+eye(4)*sqrt(-1)*omega))*C/(2*pi)); %power spectrum conditioned on input
    SE_part(k) = log(powerSpectrum_part_uncnd(k)) - log(powerSpectrum_part_cnd(k)); %integrand in the spectral efficiency integral

    powerSpectrum_full_uncnd(k) = det((((Ar-eye(4)*sqrt(-1)*omega)\(sigma^2/(gamma^2+omega^2)*brbrt + BrBrt))/(Ar'+eye(4)*sqrt(-1)*omega))/(2*pi)); %power spectrum unconditioned on input
    powerSpectrum_full_cnd(k) = det((((Ar-eye(4)*sqrt(-1)*omega)\(BrBrt))/(Ar'+eye(4)*sqrt(-1)*omega))/(2*pi)); %power spectrum conditioned on input
    SE_full(k) = log(powerSpectrum_full_uncnd(k)) - log(powerSpectrum_full_cnd(k)); %integrand in the spectral efficiency integral

end

MI_part = 1/2*real(trapz(omega_list, SE_part));
MI_full = 1/2*real(trapz(omega_list, SE_full));


%% Plots
% figure
% plot(omega_list, real(fracPOS))
% title('ACh Partially Observed: SE vs \omega','FontSize',14)
% xlabel('\omega','FontSize',14)
% ylabel('SE','FontSize',14)
% 
% figure
% plot(omega_list, real(fracFOS))
% title('ACh Fully Observed: SE vs \omega','FontSize',14)
% xlabel('\omega','FontSize',14)
% ylabel('SE','FontSize',14)




end
