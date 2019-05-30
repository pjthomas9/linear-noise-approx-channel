% Last Updated: July 26th, 2018

% Comment from PT: please rewrite as a function taking input (p) and maybe
% other parameters? from the calling function.

% QUESTION: Which parameters would be good to include in the function
% calling?  Obviously p.  We could also include the measurement vector C,
%

% also please comment as you go!
%

% Function to find the mutual information in a three state chain model 

%function [MI_FOS, MI_POS, MI_gap] = MI_3state_chain(p)
function [SE_FOS, SE_POS_low, SE_POS_high, SE_FOS_exact] = MI_3state_chain(p)

%% Seed
rng(8)


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
C = [0; 0; 1]; %Measurement vector for 3rd state observable

dt = pi*1e0; %choose "nominal" timestep


%% Build the transition matrix A

a0 = zeros(3,3);% rates with no input
a1 = zeros(3,3);% rates proportional to input

% % RANDOM rates with no input
% a0(2,1)=lognrnd(0,1); % transition rate from state 1 to state 2
% a0(1,2)=lognrnd(0,1); % transition rate from state 2 to state 1
% a0(3,2)=lognrnd(0,1); % transition rate from state 2 to state 3
% a0(2,3)=lognrnd(0,1); % transition rate from state 3 to state 2
% a0=a0-diag(sum(a0,1)); % cols sum to zero

% NONRANDOM rates with no input
a0(2,1)=10; % transition rate from state 1 to state 2
a0(1,2)=1; % transition rate from state 2 to state 1
a0(3,2)=1; % transition rate from state 2 to state 3
a0(2,3)=1; % transition rate from state 3 to state 2
a0=a0-diag(sum(a0,1)); % cols sum to zero

% % RANDOM rates proportional to input 
% a1(2,1)=0; % transition rate from state 1 to state 2
% a1(1,2)=lognrnd(0,1); % transition rate from state 2 to state 1
% a1(3,2)=0; % transition rate from state 2 to state 3
% a1(2,3)=0; % transition rate from state 3 to state 2
% a1=a1-diag(sum(a1,1)); % cols sum to zero

% NONRANDOM rates proportional to input 
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

% note [Ybar1;Ybar2;Ybar3] is the ss dist and a right eigenvector of A' 
if abs(pi1 + pi2 + pi3 - 1) > 1e-6
    warning('Steady-state distribution does not sum to one')
end

J12=pi1*abar(2,1); % steady state flux from 1 to 2. Equals ss flux from 2 to 1.
J23=pi2*abar(3,2); % steady state flux from 2 to 3. Equals ss flux from 3 to 2.


% input vector b captures effects of the input; scale by sqrt(p*(1-p))
b=sqrt(p*(1-p))*[-pi1*a1(2,1)+pi2*a1(1,2);
    pi1*a1(2,1)-pi2*(a1(1,2)+a1(3,2))+pi3*a1(2,3);
    pi2*a1(3,2)-pi3*a1(2,3)];

BBt=2*[J12,-J12,0;-J12,J12+J23,-J23;0,-J23,J23]; % noise covariance conditioned on input
BbBbt=BBt+eta*(b*b'); % unconditioned noise covariance 

%% Eigenvalues of matrix A

[v_A,d_A, u_A]=eig(abar);
d_A=diag(d_A);

for i = 1:3
    u_A(:,i) = u_A(:,i)/((u_A(:,i)'*v_A(:,i))');
end

% Find the eigenvectors associated with nonzero eigenvalues
eig_A=d_A(find(d_A<-eps(10)));
eigV_A=v_A(:,find(d_A<-eps(10)));
eigU_A=u_A(:,find(d_A<-eps(10)));

%C = eigU_A(:,1);

%disp(eigU_A'*eigV_A)


beta1_prime = C'*eigV_A(:,1)*eigU_A(:,1)'*BbBbt*eigU_A(:,1)*eigV_A(:,1)'*C;
beta2_prime = C'*eigV_A(:,2)*eigU_A(:,2)'*BbBbt*eigU_A(:,2)*eigV_A(:,2)'*C;
beta3_prime = C'*eigV_A(:,1)*eigU_A(:,1)'*BbBbt*eigU_A(:,2)*eigV_A(:,2)'*C;

beta1 = C'*eigV_A(:,1)*eigU_A(:,1)'*BBt*eigU_A(:,1)*eigV_A(:,1)'*C;
beta2 = C'*eigV_A(:,2)*eigU_A(:,2)'*BBt*eigU_A(:,2)*eigV_A(:,2)'*C;
beta3 = C'*eigV_A(:,1)*eigU_A(:,1)'*BBt*eigU_A(:,2)*eigV_A(:,2)'*C;

SE_POS_low = log((beta1_prime*eig_A(2)^2 + beta2_prime*eig_A(1)^2 + 2*beta3_prime*eig_A(1)*eig_A(2))/(beta1*eig_A(2)^2 + beta2*eig_A(1)^2 + 2*beta3*eig_A(1)*eig_A(2)));
SE_POS_high = log((beta1_prime + beta2_prime + 2*beta3_prime)/(beta1 + beta2 + 2*beta3));


%% Inidividual Terms in the MI equation

% omega = [-fliplr(logspace(-8, 1)), logspace(-8, 1)]; %frequency range
% 
% eta = 1; %Parameter scaling SNR
% 
% C = [0; 0; 1]; %Measurement vector for 3rd state observable
% 
% numTerm1 = 1./(eig_A(1)+1i .* omega) .* 1./(eig_A(1)-1i .* omega) *(C'*eigV_A(:,1)*eigU_A(:,1)'*(BBt+eta*(b*b'))*eigU_A(:,1)*eigV_A(:,1)'*C);
% numTerm2 = 1./(eig_A(1)+1i .* omega) .* 1./(eig_A(2)-1i .* omega) *(C'*eigV_A(:,1)*eigU_A(:,1)'*(BBt+eta*(b*b'))*eigU_A(:,2)*eigV_A(:,2)'*C);
% numTerm3 = 1./(eig_A(2)+1i .* omega) .* 1./(eig_A(1)-1i .* omega) *(C'*eigV_A(:,2)*eigU_A(:,2)'*(BBt+eta*(b*b'))*eigU_A(:,1)*eigV_A(:,1)'*C);
% numTerm4 = 1./(eig_A(2)+1i .* omega) .* 1./(eig_A(2)-1i .* omega) *(C'*eigV_A(:,2)*eigU_A(:,2)'*(BBt+eta*(b*b'))*eigU_A(:,2)*eigV_A(:,2)'*C);
% 
% denomTerm1 = 1./(eig_A(1)+1i .* omega) .* 1./(eig_A(1)-1i .* omega) *(C'*eigV_A(:,1)*eigU_A(:,1)'*(BBt)*eigU_A(:,1)*eigV_A(:,1)'*C);
% denomTerm2 = 1./(eig_A(1)+1i .* omega) .* 1./(eig_A(2)-1i .* omega) *(C'*eigV_A(:,1)*eigU_A(:,1)'*(BBt)*eigU_A(:,2)*eigV_A(:,2)'*C);
% denomTerm3 = 1./(eig_A(2)+1i .* omega) .* 1./(eig_A(1)-1i .* omega) *(C'*eigV_A(:,2)*eigU_A(:,2)'*(BBt)*eigU_A(:,1)*eigV_A(:,1)'*C);
% denomTerm4 = 1./(eig_A(2)+1i .* omega) .* 1./(eig_A(2)-1i .* omega) *(C'*eigV_A(:,2)*eigU_A(:,2)'*(BBt)*eigU_A(:,2)*eigV_A(:,2)'*C);


%% Calculating the Mutual Information (per unit frequency) for Partially Observed

% frac = @(omega) (1./(eig_A(1)+1i .* omega) .* 1./(eig_A(1)-1i .* omega) *(C'*eigV_A(:,1)*eigU_A(:,1)'*(BBt+eta*(b*b'))*eigU_A(:,1)*eigV_A(:,1)'*C) + ...
%     1./(eig_A(1)+1i .* omega) .* 1./(eig_A(2)-1i .* omega) *(C'*eigV_A(:,1)*eigU_A(:,1)'*(BBt+eta*(b*b'))*eigU_A(:,2)*eigV_A(:,2)'*C) + ...
%     1./(eig_A(2)+1i .* omega) .* 1./(eig_A(1)-1i .* omega) *(C'*eigV_A(:,2)*eigU_A(:,2)'*(BBt+eta*(b*b'))*eigU_A(:,1)*eigV_A(:,1)'*C) + ...
%     1./(eig_A(2)+1i .* omega) .* 1./(eig_A(2)-1i .* omega) *(C'*eigV_A(:,2)*eigU_A(:,2)'*(BBt+eta*(b*b'))*eigU_A(:,2)*eigV_A(:,2)'*C))/...
%     (1./(eig_A(1)+1i .* omega) .* 1./(eig_A(1)-1i .* omega) *(C'*eigV_A(:,1)*eigU_A(:,1)'*(BBt)*eigU_A(:,1)*eigV_A(:,1)'*C) + ...
%     1./(eig_A(1)+1i .* omega) .* 1./(eig_A(2)-1i .* omega) *(C'*eigV_A(:,1)*eigU_A(:,1)'*(BBt)*eigU_A(:,2)*eigV_A(:,2)'*C) + ...
%     1./(eig_A(2)+1i .* omega) .* 1./(eig_A(1)-1i .* omega) *(C'*eigV_A(:,2)*eigU_A(:,2)'*(BBt)*eigU_A(:,1)*eigV_A(:,1)'*C) + ...
%     1./(eig_A(2)+1i .* omega) .* 1./(eig_A(2)-1i .* omega) *(C'*eigV_A(:,2)*eigU_A(:,2)'*(BBt)*eigU_A(:,2)*eigV_A(:,2)'*C));

% By diagonalization

% numerator = @(omega) (1./(eig_A(1)+1i .* omega) .* 1./(eig_A(1)-1i .* omega) *(C'*eigV_A(:,1)*eigU_A(:,1)'*(BBt+eta*(b*b'))*eigU_A(:,1)*eigV_A(:,1)'*C) + ...
%     1./(eig_A(1)+1i .* omega) .* 1./(eig_A(2)-1i .* omega) *(C'*eigV_A(:,1)*eigU_A(:,1)'*(BBt+eta*(b*b'))*eigU_A(:,2)*eigV_A(:,2)'*C) + ...
%     1./(eig_A(2)+1i .* omega) .* 1./(eig_A(1)-1i .* omega) *(C'*eigV_A(:,2)*eigU_A(:,2)'*(BBt+eta*(b*b'))*eigU_A(:,1)*eigV_A(:,1)'*C) + ...
%     1./(eig_A(2)+1i .* omega) .* 1./(eig_A(2)-1i .* omega) *(C'*eigV_A(:,2)*eigU_A(:,2)'*(BBt+eta*(b*b'))*eigU_A(:,2)*eigV_A(:,2)'*C));
% 
% 
% denominator = @(omega) ...
%     (1./(eig_A(1)+1i .* omega) .* 1./(eig_A(1)-1i .* omega) *(C'*eigV_A(:,1)*eigU_A(:,1)'*(BBt)*eigU_A(:,1)*eigV_A(:,1)'*C) + ...
%     1./(eig_A(1)+1i .* omega) .* 1./(eig_A(2)-1i .* omega) *(C'*eigV_A(:,1)*eigU_A(:,1)'*(BBt)*eigU_A(:,2)*eigV_A(:,2)'*C) + ...
%     1./(eig_A(2)+1i .* omega) .* 1./(eig_A(1)-1i .* omega) *(C'*eigV_A(:,2)*eigU_A(:,2)'*(BBt)*eigU_A(:,1)*eigV_A(:,1)'*C) + ...
%     1./(eig_A(2)+1i .* omega) .* 1./(eig_A(2)-1i .* omega) *(C'*eigV_A(:,2)*eigU_A(:,2)'*(BBt)*eigU_A(:,2)*eigV_A(:,2)'*C));

% By direct method
omega_list =  [-fliplr(logspace(-8, 2, 1000)), logspace(-8, 2, 1000)] ; %creates a range of 2000 frequencies to numerically integrate over

fracPOS = nan(1, length(omega_list)); %integrand in the spectral efficiency integral
for k = 1:length(omega_list)
    omega = omega_list(k);
    %powerSpectrumPOS_uncnd(k) = C'*(((abar+eye(3)*sqrt(-1)*omega)\BbBbt)/(abar'-eye(3)*sqrt(-1)*omega))*C/(2*pi); %power spectrum unconditioned on input
    %powerSpectrumPOS_cnd(k) = C'*(((abar+eye(3)*sqrt(-1)*omega)\BBt)/(abar'-eye(3)*sqrt(-1)*omega))*C/(2*pi); %power spectrum conditioned on input
    
    powerSpectrumPOS_uncnd(k) = C'*(((abar+eye(3)*sqrt(-1)*omega)\BbBbt)/(abar'-eye(3)*sqrt(-1)*omega))*C/(2*pi); %power spectrum unconditioned on input
    powerSpectrumPOS_cnd(k) = C'*(((abar+eye(3)*sqrt(-1)*omega)\BBt)/(abar'-eye(3)*sqrt(-1)*omega))*C/(2*pi); %power spectrum conditioned on input
    
    ratioPOS(k) = powerSpectrumPOS_uncnd(k)/powerSpectrumPOS_cnd(k);
    SE_diagMethod(k) = (beta1_prime/(eig_A(1)^2 +omega^2) + beta2_prime/(eig_A(2)^2 + omega^2)...
        + 2*beta3_prime*(eig_A(1)*eig_A(2)+omega^2)/((eig_A(1)^2 + omega^2)*(eig_A(2)^2 +omega^2)))...
        /(beta1/(eig_A(1)^2+omega^2) + beta2/(eig_A(2)^2+omega^2)...
        + 2*beta3*(eig_A(1)*eig_A(2)+omega^2)/((eig_A(1)^2 +omega^2)*(eig_A(2)^2 +omega^2)));
    fracPOS(k) = log(powerSpectrumPOS_uncnd(k)) - log(powerSpectrumPOS_cnd(k)); %integrand in the spectral efficiency integral
end

MI_POS = real(trapz(omega_list, fracPOS)); % MI per specified frequency band (-Pi/dt to +Pi/dt)

%% Plots

% Each plot below represents one of the terms in the Spectral Efficiency
% calculation; we plot both the real (blue) and imaginary (red) parts

% plot(omega, real(numTerm1), omega, imag(numTerm1))
% plot(omega, real(numTerm2), omega, imag(numTerm2))
% plot(omega, real(numTerm3), omega, imag(numTerm3))
% plot(omega, real(numTerm4), omega, imag(numTerm4))
% 
% plot(omega, real(denomTerm1), omega, imag(denomTerm1))
% plot(omega, real(denomTerm2), omega, imag(denomTerm2))
% plot(omega, real(denomTerm3), omega, imag(denomTerm3))
% plot(omega, real(denomTerm4), omega, imag(denomTerm4))

%% Fully Observed System (FOS)

% Diagonalization Method
% numeratorFOS = @(omega) (1./(eig_A(1)+1i .* omega) .* 1./(eig_A(1)-1i .* omega) *(eigV_A(:,1)*eigU_A(:,1)'*(BBt+eta*(b*b'))*eigU_A(:,1)*eigV_A(:,1)') + ...
%     1./(eig_A(1)+1i .* omega) .* 1./(eig_A(2)-1i .* omega) *(eigV_A(:,1)*eigU_A(:,1)'*(BBt+eta*(b*b'))*eigU_A(:,2)*eigV_A(:,2)') + ...
%     1./(eig_A(2)+1i .* omega) .* 1./(eig_A(1)-1i .* omega) *(eigV_A(:,2)*eigU_A(:,2)'*(BBt+eta*(b*b'))*eigU_A(:,1)*eigV_A(:,1)') + ...
%     1./(eig_A(2)+1i .* omega) .* 1./(eig_A(2)-1i .* omega) *(eigV_A(:,2)*eigU_A(:,2)'*(BBt+eta*(b*b'))*eigU_A(:,2)*eigV_A(:,2)'));
% 
% denominatorFOS = @(omega) ...
%     (1./(eig_A(1)+1i .* omega) .* 1./(eig_A(1)-1i .* omega) *(eigV_A(:,1)*eigU_A(:,1)'*(BBt)*eigU_A(:,1)*eigV_A(:,1)') + ...
%     1./(eig_A(1)+1i .* omega) .* 1./(eig_A(2)-1i .* omega) *(eigV_A(:,1)*eigU_A(:,1)'*(BBt)*eigU_A(:,2)*eigV_A(:,2)') + ...
%     1./(eig_A(2)+1i .* omega) .* 1./(eig_A(1)-1i .* omega) *(eigV_A(:,2)*eigU_A(:,2)'*(BBt)*eigU_A(:,1)*eigV_A(:,1)') + ...
%     1./(eig_A(2)+1i .* omega) .* 1./(eig_A(2)-1i .* omega) *(eigV_A(:,2)*eigU_A(:,2)'*(BBt)*eigU_A(:,2)*eigV_A(:,2)'));

% % Direct method
% numeratorFOS = @(omega) (abar+eye(3)*sqrt(-1)*omega)\BbBbt/(abar'-eye(3)*sqrt(-1)*omega)/(2*pi);
% 
% denominatorFOS = @(omega) (abar+eye(3)*sqrt(-1)*omega)\BBt/(abar'-eye(3)*sqrt(-1)*omega)/(2*pi);
% 
% omega =  linspace(-pi/dt, pi/dt, 1000);
% 
% for i = 1:length(omega)
%     [v_num,d_num]=eig(numeratorFOS(omega(i)));
%     d_num=diag(d_num);
%     eig_num=d_num(find(d_num>eps(10)));
%     if length(eig_num) == 2
%         pdetnum(i) = eig_num(1)*eig_num(2);
%     elseif length(eig_num) == 1
%         pdetnum(i) = eig_num(1);
%     end
%     
%     [v_denom,d_denom]=eig(denominatorFOS(omega(i)));
%     d_denom=diag(d_denom);
%     eig_denom=d_denom(find(d_denom>eps(10)));
%     if length(eig_denom) == 2
%         pdetdenom(i) = eig_denom(1)*eig_denom(2);
%     elseif length(eig_denom) == 1
%         pdetdenom(i) = eig_denom(1);
%     end
% end
% 
% fracFOS = nan(1, 1000);
% for k = 1:length(omega)
%     fracFOS(k) = log(pdetnum(k))-log(pdetdenom(k));
% end
% 
% MI_FOS = real(trapz(omega, fracFOS));
% MI_gap = MI_FOS - MI_POS;

% Direct method to compute spectral efficiency for fully observed system
% (FOS)

omega_list =  [-fliplr(logspace(-8, 2, 1000)), logspace(-8, 2, 1000)] ; %creates a range of 2000 frequencies to numerically integrate over


pdet_powerSpectrumFOS_uncnd = nan(1, length(omega_list)); %pseudodeterminant of fully observed, unconditioned power spectrum
pdet_powerSpectrumFOS_cnd = nan(1, length(omega_list)); %pseudodeterminant of fully observed, conditioned power spectrum
fracFOS = nan(1, length(omega_list)); %integrand in the spectral efficiency integral



for k = 1:length(omega_list)
    omega = omega_list(k);
    powerSpectrumFOS_uncnd = (abar+eye(3)*sqrt(-1)*omega)\BbBbt/(abar'-eye(3)*sqrt(-1)*omega)/(2*pi);
    [v_num,d_num]=eig(powerSpectrumFOS_uncnd);
    d_num=diag(d_num);
    %eig_num=d_num(find(abs(d_num)>eps(10)));
    eig_num=d_num(find(abs(d_num)>1e-9));
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
    
    powerSpectrumFOS_cnd = (abar+eye(3)*sqrt(-1)*omega)\BBt/(abar'-eye(3)*sqrt(-1)*omega)/(2*pi);

    [v_denom,d_denom]=eig(powerSpectrumFOS_cnd);
    d_denom=diag(d_denom);
    %eig_denom=d_denom(find(abs(d_denom)>eps(10)));
    eig_denom=d_denom(find(abs(d_denom)>1e-9));
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
    fracFOS(k) = log(pdet_powerSpectrumFOS_uncnd(k))-log(pdet_powerSpectrumFOS_cnd(k)); %integrand in the spectral efficiency integral
end


MI_FOS = real(trapz(omega_list, fracFOS)); % MI per specified frequency band (-Pi/dt to +Pi/dt)
MI_gap = MI_FOS - MI_POS; %difference between fully and partially observed

SE_FOS = real(fracFOS(1));
%% Simplification for 3rd state observable


%integrand in the spectral efficiency integral for partially observed (POS)
%integrandPOS = log((2*pi2*abar(3,2) + sqrt(p*(1-p))*((a1(3,2)*pi2)^2 - 2*a1(3,2)*a1(2,3)*pi2*pi3 + (a1(2,3)*pi3)^2))/(2*pi2*abar(3,2)));

SE = (C'*BbBbt*C)/(C'*BBt*C);


SE_highOmegaApprox = log((beta1_prime + beta2_prime + 2*beta3_prime)/(beta1 + beta2 + 2*beta3));


SE_lowOmegaApprox = log((beta1_prime*eig_A(2)^2 + beta2_prime*eig_A(1)^2 + 2*beta3_prime*(eig_A(1)*eig_A(2)))...
    /(beta1*eig_A(2)^2 + beta2*eig_A(1)^2 + 2*beta3*(eig_A(1)*eig_A(2))));

%% Plots
%  
% figure(1)
% plot(omega_list, real(fracPOS))
% title('Three-State Chain Partially Observed: SE vs \omega','FontSize',14)
% xlabel('\omega','FontSize',14)
% ylabel('SE','FontSize',14)
% 
% figure(2)
% plot(omega_list, real(fracFOS))
% title('Three-State Chain Fully Observed: SE vs \omega','FontSize',14)
% xlabel('\omega','FontSize',14)
% ylabel('SE','FontSize',14)
% 
% 
% figure(3)
% plot(omega_list, real(fracFOS)-real(fracPOS))
% title('Three-State Chain Fully Observed: SE difference vs \omega','FontSize',14)
% xlabel('\omega','FontSize',14)
% ylabel('SE difference','FontSize',14)

%%

SE_FOS_exact = log(BbBbt(1,1)*BbBbt(3,3)/(BBt(1,1)*BBt(3,3)));
end
