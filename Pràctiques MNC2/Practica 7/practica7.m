%% Practical 7

clear; close all;
format long;

% Define the variables given by the statement
m1 = 1; m2 = 1; m3 = 1;
k1 = 1; k2 = sqrt(2); k3 = sqrt(3); k4 = 4;

% Declare B matrix
B = zeros(6, 6);
B(1:3, 4:6) = eye(3);
B(4, 1) = -(k1 + k2)/m1; B(4, 2) = k2/m1;
B(5, 1) = k2/m2; B(5, 2) = -(k2 + k3)/m2; B(5, 3) = k3/m2;
B(6, 2) = k3/m3; B(6, 3) = -(k3 + k4)/m3;

% Initial condition
z0 = [0.281, 0.033, -1.33, 1.12, 0.35, -0.299]';

%% (a) Method 1
% Let $A = e^{B\Delta t}$ so we won't have to compute the exponential of a
% matrix for each iteration. With the reasoning presented in the statement,
% we have that $z_j = Az_{j-1}$

N = 400;
dt = 0.25; % define increase of t
Tm = N*dt; % Total time 

A = expm(dt*B);
jVec = 0:N-1; tVec = dt*jVec;

zj = z0;
x2Vec = zeros(N);
for j = jVec
    x2Vec(j+1) = zj(2); % Store x_2(t)
    zj = A*zj; % Update vector z_j(t)
end

% Plot x2 through time
figure(1)
plot(tVec, x2Vec)
axis([0 Tm -2 2])
title('$x_2(t)$ position', 'Interpreter','latex')
xlabel('Time')
ylabel('$x_2(t)$', 'Interpreter','latex')

%%
% The oscillation of $x_2(t)$ can be clearly seen in Figure 1. Given the
% situation, this is what we ought to obtain. We'd expect similar results
% for the other displacement components of $z(t)$ as well as the
% speeds.

%%
fk = fft(x2Vec); % Compute DFT coefficients with FFT function

% Shift coefficients so that they are centered in omega = 0
fk_shifted = 0*fk;
fk_shifted(1:N/2) = fk(N/2+1:N);
fk_shifted(N/2+1:N) = fk(1:N/2);
omegaVec =2*pi*(-N/2:N/2-1)/Tm;

% Find peak frequencies
freqsF = zeros(1, 6); jj = 1;

% To find the frequencies of the peaks we'll go through all the coefficient
% vector and see when a coefficient is higher than those around it: this
% way we'll find local maxima. We'll store these frequencies.
for ii = 1:length(fk_shifted)
    if ii == 1
        last = 0;
    else
        last = abs(fk_shifted(ii-1));
    end
    
    if ii == length(fk_shifted)
        next = 0;
    else
        next = abs(fk_shifted(ii+1));
    end
    
    if abs(fk_shifted(ii)) > last && abs(fk_shifted(ii)) > next
        freqsF(jj) = (-N/2-1 + ii)*2*pi/Tm;
        jj = jj + 1;
    end
end

% Plot the power spectrum
% Mark in red the normal modes of oscillation (METHOD 1)
figure(2)
semilogy(omegaVec, abs(fk_shifted/N), 'k-o')
hold on
for ii = 1:length(freqsF)
    xline(freqsF(ii), 'r')
end
hold off
axis([0 8 0 1])

%%
% In Figure 2, we can see the power spectrum $|\tilde{f}_k|$ as a function 
% of $\omega _k$. We can also observe the normal modes of oscillation of 
% the system for $x_2(t)$. This frequencies are displayed above (only the 3 
% positive ones since the other three will be symmetric with respect to 
% vertical axis). The values of these peak frequency are displayed above 
% (it will be interesting to compare them to those obtained with the next 
% method).

%% (b) Method 2

eigB = eig(B); % Find eigenvalues of B
freqsB = imag(eigB); % Take imaginary part of the eigenvalues

% Plot the power spectrum
% Mark in red the normal modes of oscillation (METHOD 2)
figure(3)
semilogy(2*pi*(-N/2:N/2-1)/Tm, abs(fk_shifted/N),'k-o')
hold on
for ii = 1:length(freqsB)
    xline(freqsB(ii), 'r');
end
hold off
axis([0 8 0 1])

%%
% Figure 3 is a repetition of Figure 2 but the red vertical lines now 
% represent the normal modes of oscillation obtained with the second 
% method. We can clearly see that they coincide with the peak frequencies 
% of the power spectrum. Moreover, at first sight, it doesn't look like
% they are different from those obtained with the first method.

%%
disp(freqsF(freqsF > 0)) % Frequencies obtained with Method 1
disp(flip(freqsB(freqsB > 0)')) % Frequencies obtained with Method 2

%%
% The frequencies obtained with both methods, even though they are
% graphically located in the peaks of the power spectrum of $x_2(t)$, are a
% little different. It is very likely a numerical error in the
% frequencies computed with method 1: when computing the DFT we are already
% assuming we'll only have $N$ possible frequencies. Since we are making
% $\omega _k$ a discrete value, it will be more probable to make any
% numerical error (we'll have an error due to the method and an error due 
% to the computations). However, this differences are at most approximately 
% $0.03$ (a little lower in reality, and it happens for the lowest 
% characteristic frequency), which is about half the distance 
% $|\omega_{k+1}-\omega_k|$. Hence, this makes sense given the 
% circumstances.
%
% Since we are using $N=400$, it may not lead to such incorrect results. In 
% light of this, though, it will probably be best to use the frequencies 
% obtained with Method 2 (they may accumulate less numerical error).
