% Last Updated: Ooctboer 3rd, 2018

% Function to find the mutual information in a three state ring model 

%function [MI_FOS, MI_POS, MI_gap] = MI_3state_ring(p)


function [SE_FOS, SE_POS_low, SE_POS_high, beta1, beta2, beta3, beta4, beta1_prime, beta2_prime, beta3_prime, beta4_prime, x, y, z] = MI_3state_ring(p)

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
C = [0; 1; 0]; % Measurement vector for 2nd state observable

dt = pi*1e0; % "nominal" timestep

r=100;
s=100;

%% Build the transition matrix A

a0 = zeros(3,3);% rates with no input
a1 = zeros(3,3);% rates proportional to input

% rates with no input
a0(2,1) = 0; % transition rate from state 1 to state 2
a0(3,2) = 1; % transition rate from state 2 to state 3
a0(1,3) = 1; % transition rate from state 3 to state 1
a0 = a0-diag(sum(a0,1)); % cols sum to zero

% rates proportional to input 
a1(2,1) = 1; % transition rate from state 1 to state 2
a1(3,2) = 0; % transition rate from state 2 to state 3
a1(1,3) = 0; % transition rate from state 3 to state 1
a1 = a1-diag(sum(a1,1)); % cols sum to zero

% average transition rates (not conditioned on input)
abar = a0 + p*a1;
% note cols of abar sum to zero.


%% Compute Steady-State, b, and B
        
normfactor = abar(3,2)*abar(1,3)+abar(2,1)*abar(1,3)+abar(2,1)*abar(3,2);
Ybar1 = abar(3,2)*abar(1,3)/normfactor;
Ybar2 = abar(1,3)*abar(2,1)/normfactor;
Ybar3 = abar(2,1)*abar(3,2)/normfactor;
% [Ybar1;Ybar2;Ybar3] is the steady-state dist and a right eigenvector of A 
if abs(Ybar1 + Ybar2 + Ybar3 - 1) > 1e-6
    warning('Steady-state distribution does not sum to one')
end

% input vector b captures effects of the input; scale by sqrt(p*(1-p))
b=sqrt(p*(1-p))*[-Ybar1*a1(2,1)+Ybar2*a1(1,2);
    Ybar1*a1(2,1)-Ybar2*(a1(1,2)+a1(3,2))+Ybar3*a1(2,3);
    Ybar2*a1(3,2)-Ybar3*a1(2,3)];

B = [0, 0 , sqrt(Ybar3*abar(1, 3)); 
    sqrt(Ybar1*abar(2,1)), 0, 0;
    0, sqrt(Ybar2*abar(3, 2)), 0];

B=B-diag(sum(B,1)); % cols sum to zero

BBt=B*B'; % noise covariance conditioned on input
BbBbt=BBt+eta*(b*b'); % unconditioned noise covariance



%% Eigenvalues of matrix A

[v_A,d_A, u_A]=eig(abar);
d_A=diag(d_A);

for i = 1:3
    u_A(:,i) = u_A(:,i)/((u_A(:,i)'*v_A(:,i))');
end

% Find the eigenvectors associated with nonzero eigenvalues
eig_A=d_A(find(real(d_A)<-1e-9));
eigV_A=v_A(:,find(real(d_A)<-1e-9));
eigU_A=u_A(:,find(real(d_A)<-1e-9));

%disp(eigU_A'*eigV_A)

% %Force the eigenvalue-vector pair associated with positive imaginary part
% %to be first
% eigV_A_placeholder = eigV_A;
% eigU_A_placeholder = eigU_A;
% eigA_placeholder = eig_A;
% 
% if imag(eig_A(1)) < 0
%     eigV_A(:,2) = eigV_A_placeholder(:,1);
%     eigV_A(:,1) = eigV_A_placeholder(:,2);
%     eigU_A(:,2) = eigU_A_placeholder(:,1);
%     eigU_A(:,1) = eigU_A_placeholder(:,2);
%     eig_A(1) = eigA_placeholder(2);
%     eig_A(2) = eigA_placeholder(1);
% end
% 


%C = eigU_A(:,1);


beta1_prime = (C'*eigV_A(:,1))*(eigU_A(:,1)'*BbBbt*eigU_A(:,1))*(eigV_A(:,1)'*C);
beta2_prime = (C'*eigV_A(:,2))*(eigU_A(:,2)'*BbBbt*eigU_A(:,2))*(eigV_A(:,2)'*C);
beta3_prime = (C'*eigV_A(:,1))*(eigU_A(:,1)'*BbBbt*eigU_A(:,2))*(eigV_A(:,2)'*C);
beta4_prime = (C'*eigV_A(:,2))*(eigU_A(:,2)'*BbBbt*eigU_A(:,1))*(eigV_A(:,1)'*C);

x = eigV_A(:,1)*eigU_A(:,1)';
y = eigV_A(:,1)*eigU_A(:,1)';
z = eigV_A(:,1)*eigU_A(:,1)';

x = x(2,1);
y = y(2,2);
z = z(2,3);

beta1 = (C'*eigV_A(:,1))*(eigU_A(:,1)'*BBt*eigU_A(:,1))*(eigV_A(:,1)'*C);
beta2 = (C'*eigV_A(:,2))*(eigU_A(:,2)'*BBt*eigU_A(:,2))*(eigV_A(:,2)'*C);
beta3 = (C'*eigV_A(:,1))*(eigU_A(:,1)'*BBt*eigU_A(:,2))*(eigV_A(:,2)'*C);
beta4 = (C'*eigV_A(:,2))*(eigU_A(:,2)'*BBt*eigU_A(:,1))*(eigV_A(:,1)'*C);

%These equations below will need to be rederived since they are only valid
%for real eigenvalues.
% SE_POS_low = log((beta1_prime*eig_A(2)^2 + beta2_prime*eig_A(1)^2 + beta3_prime*eig_A(1)*eig_A(2) + beta4_prime*eig_A(2)*eig_A(1))...
%     /(beta1*eig_A(2)^2 + beta2*eig_A(1)^2 + beta3*eig_A(1)*eig_A(2) + beta4*eig_A(2)*eig_A(1)));
SE_POS_high = log((beta1_prime + beta2_prime + beta3_prime + beta4_prime)/(beta1 + beta2 + beta3 + beta4));
% 
% 
% SE_POS_low_alt = log((beta1_prime/eig_A(1)^2 + beta2_prime/eig_A(2)^2 + beta3_prime/(eig_A(1)*eig_A(2)) + beta4_prime/(eig_A(2)*eig_A(1)))...
%     /(beta1/eig_A(1)^2 + beta2/eig_A(2)^2 + beta3/(eig_A(1)*eig_A(2)) + beta4/(eig_A(2)*eig_A(1))));


SE_POS_low = log((beta1_prime/abs(eig_A(1))^2 + beta2_prime/abs(eig_A(2))^2 + beta3_prime/(eig_A(1)*conj(eig_A(2))) + beta4_prime/(eig_A(2)*conj(eig_A(1))))...
     /(beta1/abs(eig_A(1))^2 + beta2/abs(eig_A(2))^2 + beta3/(eig_A(1)*conj(eig_A(2))) + beta4/(eig_A(2)*conj(eig_A(1)))));


%% Partially Observed System (POS)

omega_list =  linspace(-10000, 10000, 1000); %creates a range of 1000 frequencies to numerically integrate over


powerSpectrumPOS_uncnd = nan(1, length(omega_list));
powerSpectrumPOS_cnd = nan(1, length(omega_list));
fracPOS = nan(1, length(omega_list));
for k = 1:length(omega_list)
    omega = omega_list(k);
    powerSpectrumPOS_uncnd(k) = C'*(((abar+eye(3)*sqrt(-1)*omega)\BbBbt)/(abar'-eye(3)*sqrt(-1)*omega))*C/(2*pi); %power spectrum unconditioned on input
    powerSpectrumPOS_cnd(k) = C'*(((abar+eye(3)*sqrt(-1)*omega)\BBt)/(abar'-eye(3)*sqrt(-1)*omega))*C/(2*pi); %power spectrum conditioned on input
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
    powerSpectrumFOS_uncnd = (abar+eye(3)*sqrt(-1)*omega)\BbBbt/(abar'-eye(3)*sqrt(-1)*omega)/(2*pi); %fully observed, unconditioned power spectrum
    
    % Computes the pseudo-determinant of the unconditioned power spectrum by finding the product of the nonzero eigenvalues
    [v_num,d_num]=eig(powerSpectrumFOS_uncnd);
    d_num=diag(d_num);
    eig_num=d_num(find(abs(d_num)>eps(10)));
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
    
    powerSpectrumFOS_cnd = (abar+eye(3)*sqrt(-1)*omega)\BBt/(abar'-eye(3)*sqrt(-1)*omega)/(2*pi); %fully observed, conditioned power spectrum

    % Computes the pseudo-determinant of the unconditioned power spectrum by finding the product of the nonzero eigenvalues
    [v_denom,d_denom]=eig(powerSpectrumFOS_cnd);
    d_denom=diag(d_denom);
    eig_denom=d_denom(find(abs(d_denom)>eps(10))); %Note: eps(10) returns the distance from 10.0 to the next largest double-precision number
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
end

for k = 1:length(omega_list)
    fracFOS(k) = (log(pdet_powerSpectrumFOS_uncnd(k))-log(pdet_powerSpectrumFOS_cnd(k)))*(r+s)/(pi*((r+s)^2 + omega_list(k)^2)); %integrand in the spectral efficiency integral
end

MI_FOS = real(trapz(omega_list, fracFOS)); % MI per specified frequency band (-Pi/dt to +Pi/dt)
MI_gap = MI_FOS - MI_POS; %difference between fully and partially observed

SE_FOS = real(fracFOS(1));


%%

evalsBBt = eig(BBt);
evalsBbBbt = eig(BbBbt);
pdetBBt = evalsBBt(2)*evalsBBt(3);
pdetBbBbt = evalsBbBbt(2)*evalsBbBbt(3);
ratio = pdetBbBbt/pdetBBt;
integrand  = log(ratio);



%ratioPOS = (C'*BbBbt*C)/(C'*BBt*C);

%% Plots
% figure
% plot(omega_list, real(fracPOS))
% title('Three-State Ring Partially Observed: SE vs \omega','FontSize',14)
% xlabel('\omega','FontSize',14)
% ylabel('SE','FontSize',14)
% 
% figure
% plot(omega_list, real(fracFOS))
% title('Three-State Ring Fully Observed: SE vs \omega','FontSize',14)
% xlabel('\omega','FontSize',14)
% ylabel('SE','FontSize',14)


end
