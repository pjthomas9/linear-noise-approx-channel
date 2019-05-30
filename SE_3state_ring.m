% Last Updated: May 30th, 2019

% Function to find spectral efficiency for three-state ring (e.g. ChR2)

function [SE_full, SE_part_low, SE_part_high] = SE_3state_ring(p)

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

%% Build the transition matrix A

a0 = zeros(3,3);% rates with no input
a1 = zeros(3,3);% rates proportional to input

% rates with no input
a0(2,1) = 0; % transition rate from state 1 to state 2
a0(3,2) = 50; % transition rate from state 2 to state 3
a0(1,3) = 17; % transition rate from state 3 to state 1
a0 = a0-diag(sum(a0,1)); % cols sum to zero

% rates proportional to input 
a1(2,1) = 5e-3; % transition rate from state 1 to state 2
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

%% Partially Observed System (SE_part)

beta1_prime = (C'*eigV_A(:,1))*(eigU_A(:,1)'*BbBbt*eigU_A(:,1))*(eigV_A(:,1)'*C);
beta2_prime = (C'*eigV_A(:,2))*(eigU_A(:,2)'*BbBbt*eigU_A(:,2))*(eigV_A(:,2)'*C);
beta3_prime = (C'*eigV_A(:,1))*(eigU_A(:,1)'*BbBbt*eigU_A(:,2))*(eigV_A(:,2)'*C);
beta4_prime = (C'*eigV_A(:,2))*(eigU_A(:,2)'*BbBbt*eigU_A(:,1))*(eigV_A(:,1)'*C);

beta1 = (C'*eigV_A(:,1))*(eigU_A(:,1)'*BBt*eigU_A(:,1))*(eigV_A(:,1)'*C);
beta2 = (C'*eigV_A(:,2))*(eigU_A(:,2)'*BBt*eigU_A(:,2))*(eigV_A(:,2)'*C);
beta3 = (C'*eigV_A(:,1))*(eigU_A(:,1)'*BBt*eigU_A(:,2))*(eigV_A(:,2)'*C);
beta4 = (C'*eigV_A(:,2))*(eigU_A(:,2)'*BBt*eigU_A(:,1))*(eigV_A(:,1)'*C);

SE_part_high = log((beta1_prime + beta2_prime + beta3_prime + beta4_prime)/(beta1 + beta2 + beta3 + beta4));

SE_part_low = log((beta1_prime/abs(eig_A(1))^2 + beta2_prime/abs(eig_A(2))^2 + beta3_prime/(eig_A(1)*conj(eig_A(2))) + beta4_prime/(eig_A(2)*conj(eig_A(1))))...
     /(beta1/abs(eig_A(1))^2 + beta2/abs(eig_A(2))^2 + beta3/(eig_A(1)*conj(eig_A(2))) + beta4/(eig_A(2)*conj(eig_A(1)))));

%% Fully Observed System (SE_full)

[~,d_BbBbt]=eig(BbBbt);
d_BbBbt=diag(d_BbBbt);
eig_num=d_BbBbt(find(abs(d_BbBbt)>eps(10)));
pdet_BbBbt = eig_num(1)*eig_num(2);

[~,d_BBt]=eig(BBt);
d_BBt=diag(d_BBt);
eig_num=d_BBt(find(abs(d_BBt)>eps(10)));
pdet_BBt = eig_num(1)*eig_num(2);

SE_full = (log(pdet_BbBbt)-log(pdet_BBt));

%% Plots for SE vs omega
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
