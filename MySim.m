clc;clear all;
M0a = 1; M0b = 1;
m0 = [0; 0; M0a; 0; 0; M0b];
dt = 1e-5; % ms
B1 = 1 * 1e-6; % T for cw
gamma = 42.58* 2 * pi * 1e6; % Mhz = 1e6 Hz
dwa = 0;
dwb = 9.4 * 2 * pi;
t = 1;
N = ceil(t / dt);
t_eff = 0.5 * N;

% cw
B1_array = ones(N, 1); % cw 
B1_array(1:t_eff,1) = B1 * B1_array(1:t_eff,1);
B1_array(t_eff:N,1) = 0 * B1_array(t_eff:N,1);

%theta = atan(w1/dw)
theta = atan(B1*gamma * 1e-3 / dwb);
% spin_lock


% w1 = gamma * B1;
w1_array = B1_array * gamma;
T1s = 0.77;
T2s = 0.033;

T1w = 1.5;
T2w = 0.06;

R1b = 1/T1s; R1a = 1/T1w; R2b = 1/T2s; R2a = 1/T2w;

fb = 0.01;
kb = 100; % 100 Hz
% ksw = k; kws = f * k;

M_data = zeros(6, N);
% M(:, 1) = m0;
M_data(:, 1) = m0;
for ii = 1:(N-1)
    M = M_data(:, ii);
    % dmxs = -dw_s * M(2) - R2s*M(1) - ksw*M(1) +kws*M(4);
    % dmys = dw_s * M(1) + B1_array(ii)*M(3) - R2s*M(2) - ksw*M(2) + kws*M(5);
    % dmzs = -B1_array(ii) * M(2) / gamma - R1s * (M(3) - M0s) - ksw*M(3) + kws * M(6);
    % dmxw = -dw * M(5) - R2w*M(4) + ksw*M(1) - kws*M(4);
    % dmyw = dw * M(4) + B1_array(ii)*M(6) - R2w*M(5) + ksw*M(2) - kws*M(5);
    % dmzw = -B1_array(ii) * M(5) / gamma - R1w * (M(6) - M0s) + ksw*M(3) - kws*M(6);

    La = [- R2a -dwa 0;
        +dwa -R2a +w1_array(ii);
        0 -w1_array(ii) -R1a];
    Lb = [-R2b -dwb 0;
        +dwb -R2b +w1_array(ii);
        0 -w1_array(ii) -R1b];
    Kb = eye(3) * kb;
    Ka = fb * Kb;
    A = [La-Ka Kb;
        Ka Lb-Kb];
    % A = A / 1000;
    C = [0;0;R1a*M0a;0;0;R1b*M0b];
    % C = C / 1000;

    M = A * M + C;

    if ii == t_eff/dt
        M
    end
    M_data(:, ii+1) = M;
end

figure
plot3(M_data(1,:), M_data(2,:), M_data(3,:));
% plot3(M_data(4,:), M_data(5,:), M_data(6,:));