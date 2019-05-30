% Last Updated: July 20th, 2018

% Function to find the mutual information for ACh

%function [MI_FOS, MI_POS, MI_gap] = MI_ACh(p)

function [SE_FOS, SE_POS_low, SE_POS_high] = MI_ACh(p)

%% Input Parameters
% Constants we must specify
%p; %prob of high input
% N_x = 1000; %number of channels
% M = ???; %number of sources

% Derived Constants
% epsilon = 1/sqrt(N_x);
% D = ??? %diffusion constant; dependent on M
% eta = 2*D/epsilon^2

% Since we have not specified M (and thus D) we will choose eta
eta = 1;

% Measurement Vector
C = [0; 0; 1; 1; 1]; % Measurement vector for 2nd state observable

dt = pi*1e-2; % "nominal" timestep

r = 100;
s = 100;

%% Build the transition matrix A

a0 = zeros(5,5);% rates with no input
a1 = zeros(5,5);% rates proportional to input

% rates with no input
a0(2,1) = 0.01; % transition rate from state 1 to state 2
a0(1,2) = 0.66; % transition rate from state 2 to state 1
a0(4,1) = 3e3; % transition rate from state 1 to state 4
a0(1,4) = 15; % transition rate from state 4 to state 1
a0(3,2) = 5e2; % transition rate from state 2 to state 3
a0(2,3) = 1.5e4; % transition rate from state 3 to state 2
a0(4,3) = 4e3; % transition rate from state 3 to state 4
a0(3,4) = 0.01; % transition rate from state 4 to state 3
a0(5,4) = 2e3; % transition rate from state 4 to state 5
a0(4,5) = 0.01; % transition rate from state 5 to state 4
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
b=sqrt(p*(1-p))*[-pi1*(a1(2,1)+a1(4,1))+pi2*a1(1,2)+pi4*a1(1, 4);
    pi1*a1(2,1)-pi2*(a1(1,2)+a1(3,2))+pi3*a1(2,3);
    pi2*a1(3,2)-pi3*(a1(2,3)+a1(4,3))+pi4*a1(3,4);
    pi3*a1(4,3)-pi4*(a1(1,4)+a1(3,4)+a1(5,4))+pi5*a1(4,5);
    pi4*a1(5,4)-pi5*a1(4,5)];

B = [-sqrt(pi1*abar(2,1)), sqrt(pi2*abar(1,2)), 0, 0, 0, 0, -sqrt(pi1*abar(4,1)), sqrt(pi4*abar(1,4)), 0, 0;
    sqrt(pi1*abar(2,1)), -sqrt(pi2*abar(1,2)), -sqrt(pi2*abar(3,2)), sqrt(pi3*abar(2,3)), 0, 0, 0, 0, 0, 0;
    0, 0, sqrt(pi2*abar(3,2)), -sqrt(pi3*abar(2,3)), -sqrt(pi3*abar(4,3)), sqrt(pi4*abar(3,4)), 0, 0, 0, 0;
    0, 0, 0, 0, sqrt(pi3*abar(4,3)), -sqrt(pi4*abar(3,4)), sqrt(pi1*abar(4,1)), -sqrt(pi4*abar(1,4)), -sqrt(pi4*abar(5,4)), sqrt(pi5*abar(4,5));
    0, 0, 0, 0, 0, 0, 0, 0, sqrt(pi4*abar(5, 4)), -sqrt(pi5*abar(4,5))];

BBt=B*B'; % noise covariance conditioned on input
bbt = b*b';
BbBbt=BBt+eta*(b*b'); % unconditioned noise covariance 

%%% PT commented out everything from this point to the end...

%% Eigenvalues of matrix A

[v_A,d_A, u_A]=eig(abar);
d_A=diag(d_A);

for i = 1:5
    u_A(:,i) = u_A(:,i)/((u_A(:,i)'*v_A(:,i))');
end

% Find the eigenvectors associated with nonzero eigenvalues
eig_A=d_A(find(d_A<-1e-9));
eigV_A=v_A(:,find(d_A<-1e-9));
eigU_A=u_A(:,find(d_A<-1e-9));


%disp(eigU_A'*eigV_A)


beta11_prime = C'*eigV_A(:,1)*eigU_A(:,1)'*BbBbt*eigU_A(:,1)*eigV_A(:,1)'*C;
beta12_prime = C'*eigV_A(:,1)*eigU_A(:,1)'*BbBbt*eigU_A(:,2)*eigV_A(:,2)'*C;
beta13_prime = C'*eigV_A(:,1)*eigU_A(:,1)'*BbBbt*eigU_A(:,3)*eigV_A(:,3)'*C;
beta14_prime = C'*eigV_A(:,1)*eigU_A(:,1)'*BbBbt*eigU_A(:,4)*eigV_A(:,4)'*C;
beta21_prime = C'*eigV_A(:,2)*eigU_A(:,2)'*BbBbt*eigU_A(:,1)*eigV_A(:,1)'*C;
beta22_prime = C'*eigV_A(:,2)*eigU_A(:,2)'*BbBbt*eigU_A(:,2)*eigV_A(:,2)'*C;
beta23_prime = C'*eigV_A(:,2)*eigU_A(:,2)'*BbBbt*eigU_A(:,3)*eigV_A(:,3)'*C;
beta24_prime = C'*eigV_A(:,2)*eigU_A(:,2)'*BbBbt*eigU_A(:,4)*eigV_A(:,4)'*C;
beta31_prime = C'*eigV_A(:,3)*eigU_A(:,3)'*BbBbt*eigU_A(:,1)*eigV_A(:,1)'*C;
beta32_prime = C'*eigV_A(:,3)*eigU_A(:,3)'*BbBbt*eigU_A(:,2)*eigV_A(:,2)'*C;
beta33_prime = C'*eigV_A(:,3)*eigU_A(:,3)'*BbBbt*eigU_A(:,3)*eigV_A(:,3)'*C;
beta34_prime = C'*eigV_A(:,3)*eigU_A(:,3)'*BbBbt*eigU_A(:,4)*eigV_A(:,4)'*C;
beta41_prime = C'*eigV_A(:,4)*eigU_A(:,4)'*BbBbt*eigU_A(:,1)*eigV_A(:,1)'*C;
beta42_prime = C'*eigV_A(:,4)*eigU_A(:,4)'*BbBbt*eigU_A(:,2)*eigV_A(:,2)'*C;
beta43_prime = C'*eigV_A(:,4)*eigU_A(:,4)'*BbBbt*eigU_A(:,3)*eigV_A(:,3)'*C;
beta44_prime = C'*eigV_A(:,4)*eigU_A(:,4)'*BbBbt*eigU_A(:,4)*eigV_A(:,4)'*C;

beta11 = C'*eigV_A(:,1)*eigU_A(:,1)'*BBt*eigU_A(:,1)*eigV_A(:,1)'*C;
beta12 = C'*eigV_A(:,1)*eigU_A(:,1)'*BBt*eigU_A(:,2)*eigV_A(:,2)'*C;
beta13 = C'*eigV_A(:,1)*eigU_A(:,1)'*BBt*eigU_A(:,3)*eigV_A(:,3)'*C;
beta14 = C'*eigV_A(:,1)*eigU_A(:,1)'*BBt*eigU_A(:,4)*eigV_A(:,4)'*C;
beta21 = C'*eigV_A(:,2)*eigU_A(:,2)'*BBt*eigU_A(:,1)*eigV_A(:,1)'*C;
beta22 = C'*eigV_A(:,2)*eigU_A(:,2)'*BBt*eigU_A(:,2)*eigV_A(:,2)'*C;
beta23 = C'*eigV_A(:,2)*eigU_A(:,2)'*BBt*eigU_A(:,3)*eigV_A(:,3)'*C;
beta24 = C'*eigV_A(:,2)*eigU_A(:,2)'*BBt*eigU_A(:,4)*eigV_A(:,4)'*C;
beta31 = C'*eigV_A(:,3)*eigU_A(:,3)'*BBt*eigU_A(:,1)*eigV_A(:,1)'*C;
beta32 = C'*eigV_A(:,3)*eigU_A(:,3)'*BBt*eigU_A(:,2)*eigV_A(:,2)'*C;
beta33 = C'*eigV_A(:,3)*eigU_A(:,3)'*BBt*eigU_A(:,3)*eigV_A(:,3)'*C;
beta34 = C'*eigV_A(:,3)*eigU_A(:,3)'*BBt*eigU_A(:,4)*eigV_A(:,4)'*C;
beta41 = C'*eigV_A(:,4)*eigU_A(:,4)'*BBt*eigU_A(:,1)*eigV_A(:,1)'*C;
beta42 = C'*eigV_A(:,4)*eigU_A(:,4)'*BBt*eigU_A(:,2)*eigV_A(:,2)'*C;
beta43 = C'*eigV_A(:,4)*eigU_A(:,4)'*BBt*eigU_A(:,3)*eigV_A(:,3)'*C;
beta44 = C'*eigV_A(:,4)*eigU_A(:,4)'*BBt*eigU_A(:,4)*eigV_A(:,4)'*C;

SE_POS_low = log((beta11_prime/(eig_A(1)*eig_A(1)) + beta12_prime/(eig_A(1)*eig_A(2))+...
    beta13_prime/(eig_A(1)*eig_A(3)) + beta14_prime/(eig_A(1)*eig_A(4))+...
    beta21_prime/(eig_A(2)*eig_A(1)) + beta22_prime/(eig_A(2)*eig_A(2))+...
    beta23_prime/(eig_A(2)*eig_A(3)) + beta24_prime/(eig_A(2)*eig_A(4))+...
    beta31_prime/(eig_A(3)*eig_A(1)) + beta32_prime/(eig_A(3)*eig_A(2))+...
    beta33_prime/(eig_A(3)*eig_A(3)) + beta34_prime/(eig_A(3)*eig_A(4))+...
    beta41_prime/(eig_A(4)*eig_A(1)) + beta42_prime/(eig_A(4)*eig_A(2))+...
    beta43_prime/(eig_A(4)*eig_A(3)) + beta44_prime/(eig_A(4)*eig_A(4)))/...
    (beta11/(eig_A(1)*eig_A(1)) + beta12/(eig_A(1)*eig_A(2))+...
    beta13/(eig_A(1)*eig_A(3)) + beta14/(eig_A(1)*eig_A(4))+...
    beta21/(eig_A(2)*eig_A(1)) + beta22/(eig_A(2)*eig_A(2))+...
    beta23/(eig_A(2)*eig_A(3)) + beta24/(eig_A(2)*eig_A(4))+...
    beta31/(eig_A(3)*eig_A(1)) + beta32/(eig_A(3)*eig_A(2))+...
    beta33/(eig_A(3)*eig_A(3)) + beta34/(eig_A(3)*eig_A(4))+...
    beta41/(eig_A(4)*eig_A(1)) + beta42/(eig_A(4)*eig_A(2))+...
    beta43/(eig_A(4)*eig_A(3)) + beta44/(eig_A(4)*eig_A(4))));

SE_POS_high = log((beta11_prime + beta12_prime + beta13_prime + beta14_prime +...
    beta21_prime + beta22_prime + beta23_prime + beta24_prime +...
    beta31_prime + beta32_prime + beta33_prime + beta34_prime +...
    beta41_prime + beta42_prime + beta43_prime + beta44_prime)/...
    (beta11 + beta12 + beta13 + beta14 + beta21 + beta22 + beta23 + beta24 +...
    beta31 + beta32 + beta33 + beta34 + beta41 + beta42 + beta43 + beta44)); 



%% Partially Observed System (POS)

omega_list =  linspace(-10000, 10000, 1000); %creates a range of 1000 frequencies to numerically integrate over

fracPOS = nan(1, length(omega_list));
for k = 1:length(omega_list)
    omega = omega_list(k);
    powerSpectrumPOS_uncnd(k) = C'*(((abar+eye(5)*sqrt(-1)*omega)\BbBbt)/(abar'-eye(5)*sqrt(-1)*omega))*C/(2*pi); %power spectrum unconditioned on input
    powerSpectrumPOS_cnd(k) = C'*(((abar+eye(5)*sqrt(-1)*omega)\BBt)/(abar'-eye(5)*sqrt(-1)*omega))*C/(2*pi); %power spectrum conditioned on input
    fracPOS(k) = (log(powerSpectrumPOS_uncnd(k)) - log(powerSpectrumPOS_cnd(k)))*(r+s)/(pi*((r+s)^2 + omega_list(k)^2)); %integrand in the spectral efficiency integral
end

MI_POS = real(trapz(omega_list, fracPOS)); % MI per specified frequency band (-Pi/dt to +Pi/dt)

%% Fully Observed System (FOS)

% Direct method to compute spectral efficiency for fully observed system
% (FOS)

omega_list =  linspace(-10000, 10000, 1000); %creates a range of 1000 frequencies to numerically integrate over


pdet_powerSpectrumFOS_uncnd = nan(1, length(omega_list)); %pseudodeterminant of fully observed, unconditioned power spectrum
pdet_powerSpectrumFOS_cnd = nan(1, length(omega_list)); %pseudodeterminant of fully observed, conditioned power spectrum
fracFOS = nan(1, length(omega_list)); %integrand in the spectral efficiency integral



for k = 1:length(omega_list)
    omega = omega_list(k);
    powerSpectrumFOS_uncnd = (abar+eye(5)*sqrt(-1)*omega)\BbBbt/(abar'-eye(5)*sqrt(-1)*omega)/(2*pi); %fully observed, unconditioned power spectrum
    
    % Computes the pseudo-determinant of the unconditioned power spectrum by finding the product of the nonzero eigenvalues
    [v_num,d_num]=eig(powerSpectrumFOS_uncnd);
    d_num=diag(d_num);
    eig_num=d_num(find(abs(d_num)>1e-12));
    if length(eig_num) == 2
        pdet_powerSpectrumFOS_uncnd(k) = eig_num(1)*eig_num(2);
    end
    if length(eig_num) == 1
        pdet_powerSpectrumFOS_uncnd(k) = eig_num(1);
    end
    if length(eig_num) == 0
        pdet_powerSpectrumFOS_uncnd(k) = 0;
    end
    if length(eig_num) == 3
        pdet_powerSpectrumFOS_uncnd(k) = eig_num(1)*eig_num(2)*eig_num(3);
    end
    if length(eig_num) == 4
        pdet_powerSpectrumFOS_uncnd(k) = eig_num(1)*eig_num(2)*eig_num(3)*eig_num(4);
    end
    if length(eig_num) == 5
        pdet_powerSpectrumFOS_uncnd(k) = eig_num(1)*eig_num(2)*eig_num(3)*eig_num(4)*eig_num(5);
    end
    
    powerSpectrumFOS_cnd = (abar+eye(5)*sqrt(-1)*omega)\BBt/(abar'-eye(5)*sqrt(-1)*omega)/(2*pi); %fully observed, conditioned power spectrum

    % Computes the pseudo-determinant of the unconditioned power spectrum by finding the product of the nonzero eigenvalues
    [v_denom,d_denom]=eig(powerSpectrumFOS_cnd);
    d_denom=diag(d_denom);
    eig_denom=d_denom(find(abs(d_denom)>1e-12)); %Note: eps(10) returns the distance from 10.0 to the next largest double-precision number
    if length(eig_denom) == 2
        pdet_powerSpectrumFOS_cnd(k) = eig_denom(1)*eig_denom(2);
    end
    if length(eig_denom) == 1
        pdet_powerSpectrumFOS_cnd(k) = eig_denom(1);
    end
    if length(eig_denom) == 0
        pdet_powerSpectrumFOS_cnd(k) = 0;
    end
    if length(eig_denom) == 3
        pdet_powerSpectrumFOS_cnd(k) = eig_denom(1)*eig_denom(2)*eig_denom(3);
    end
    if length(eig_denom) == 4
        pdet_powerSpectrumFOS_cnd(k) = eig_denom(1)*eig_denom(2)*eig_denom(3)*eig_denom(4);
    end
    if length(eig_denom) == 5
        pdet_powerSpectrumFOS_cnd(k) = eig_denom(1)*eig_denom(2)*eig_denom(3)*eig_denom(4)*eig_denom(5);
    end 
end

for k = 1:length(omega_list)
    fracFOS(k) = (log(pdet_powerSpectrumFOS_uncnd(k))-log(pdet_powerSpectrumFOS_cnd(k)))*(r+s)/(pi*((r+s)^2 + omega_list(k)^2)); %integrand in the spectral efficiency integral
end

MI_FOS = real(trapz(omega_list, fracFOS)); % MI per specified frequency band (-Pi/dt to +Pi/dt)
MI_gap = MI_FOS - MI_POS; %difference between fully and partially observed


SE_FOS = fracFOS(1);

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
