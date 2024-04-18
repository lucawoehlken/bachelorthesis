%willems lemma simulation Luca Wöhlken 15/04/24
clc;
clearvars;
%parameters
nu = 2; %initial input sequenze
T = 0.2; %Sample time
L = 100; %Online Data
N = 2*L+1; %Number of samples/ Offline Data
u_t = ones(N,1); %Input Signal
u_n = 0.1*randn(N, 1); %random noise/offline input
t = linspace(0,T*N,100); %Time Vector

%Continuous time system
Ak = [-0.2, -1; 1, 0];
Bk = [1; 0];
Ck = [1, 0];
Dk = 0;

sysk = ss(Ak, Bk, Ck, Dk);

%discrete time system model
sysd = c2d(sysk, T);
[A, B, C, D] = ssdata(sysd);

y_t = lsim(sysd, u_t); %real output signal
y_n = lsim(sysd, u_n); %offline output

%willems lemma
%Hankel matrix
HL_u = H(L, u_n(1:N));
Hnu_y = H(nu, y_n(1:N-L+nu));
Huy = [HL_u; Hnu_y];
%initial condition
u_i = u_t;
u_i = circshift (u_i,nu,1);
u_i(1:nu) = zeros (nu,1);
y_i = zeros(nu, 1);
uy = [u_i(1:L); y_i];
%alpha calculation
alpha = pinv(Huy)*uy;
%online output data
y_wl = H(L, y_n(1:N))*alpha;
y_wl = circshift(y_wl,-nu);

%graph plotts
plot(t(1:L-1), y_t(1:L-1),'LineWidth', 4); %plot of real system
hold on;
plot(t(1:L-1), y_wl(1:L-1), 'LineWidth', 2); %plot of data based system
legend('real system', 'data based system')
hold off;

%hankel-matrix function
function Hankel = H (L,z)

  N = size (z,1);     % one column = one signal element's sequence
  m = size (z,2);     % one row = one non-scalar sequence element
  Hankel = zeros (L*m,N-L+1);   % Initialize array of zeros with correct size

  for i = 1:1:N-L+1   % all columns of the Hankel matrix
    for j = 0:1:L-1   % all (block) rows of the Hankel matrix
      Hankel(j*m+1:j*m+m,i) = z(i+j,:)';
    end
  end

end