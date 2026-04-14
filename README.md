# ADC-Project
%% =====================  PART 1  ==========================

%%% User Input %%%
M = input('enter the number of signals:');
Length = input('enter the length of the largest signal: ');
fprintf('Enter the signals as in a row\n');
Result = zeros(M, Length);
t = 0:Length-1;

for i = 1: M
    Result(i,:) = input("");
end

NumOfN = rank(Result);

%%% Gram-Schmidt method %%%

S = Result';
[length, M] = size(S);
Q = zeros(length, M);

for i = 1:M
    v = S(:,i);
    for j = 1:i-1
        v = v - (Q(:,j)' * v) * Q(:,j);
    end
    if norm(v) > 1e-10   % <-- guard against zero vectors (linearly dependent)
        Q(:,i) = v / norm(v);
    end
end

N = Q';   % Transpose of Q
display(N)

%%% Plot the M input signals %%%
figure;
for i = 1: M
    subplot(M,1,i);
    stem(Result(i,:), 'filled');
    title(['signal ', num2str(i)]);
end

%%% plot the N basis functions %%%
figure;
for i = 1: NumOfN
    subplot(NumOfN,1,i);
    stairs([t t(end)+1], [N(i,:) N(i,end)], 'LineWidth', 2);
    title(['basis function ', num2str(i)]);
end

%%% Plot the constellation diagram of the signals %%%

%In Case of 2D plot
if NumOfN <= 2
 figure;
 plot(Result(:,1), Result(:,2), 'X', 'MarkerSize', 10, 'LineWidth', 2);
 grid on;
 xlabel('\phi_1');
 ylabel('\phi_2');
 title('Constellation Diagram (2D)');

%In Case of 3D plot
else
 figure;
 plot3(Result(:,1), Result(:,2), Result(:,3), 'X', 'MarkerSize', 10, 'LineWidth', 2);
 grid on;
 hold on;
 xlabel('\phi_1');
 ylabel('\phi_2');
 zlabel('\phi_3');
 title('3D Constellation Diagram');
 view(3);
 hold off;
end

%%% Calculate the energy of each symbols using the constellation diagram %%%

E = sum(Result.^2, 2);

fprintf('\nEnergY of signal:\n');
for i = 1:M
    fprintf('Signal %d energy = %.2f\n', i, E(i));
end

%% ===================== PART 2 ==========================

Tb = 1;
fs = 100;
dt = Tb / fs;
Nb = 5000;
EbNo_dB = -10:2:10;

% Unipolar NRZ symbols
s1_uni = ones(1, fs);
s2_uni = zeros(1, fs);

% Manchester symbols
s1_man = [ones(1, fs/2), -ones(1, fs/2)];
s2_man = [-ones(1, fs/2), ones(1, fs/2)];

% Gram-Schmidt for Unipolar NRZ
E1_uni = sum(s1_uni.^2) * dt;          
phi1_uni = s1_uni / sqrt(E1_uni);      % basis function

s11_uni = sum(s1_uni .* phi1_uni) * dt; 
s21_uni = sum(s2_uni .* phi1_uni) * dt; 

Eb_uni = (s11_uni^2 + s21_uni^2) / 2;  % average energy per bit

% Gram-Schmidt for Manchester 
E1_man = sum(s1_man.^2) * dt;          
phi1_man = s1_man / sqrt(E1_man);

s11_man = sum(s1_man .* phi1_man) * dt; 
s21_man = sum(s2_man .* phi1_man) * dt; 

Eb_man = (s11_man^2 + s21_man^2) / 2;  % average energy = 1

% Constellation diagram 
figure;
subplot(1,2,1);
plot([s21_uni s11_uni], [0 0], 'bX', 'MarkerSize', 14, 'LineWidth', 2);
xline((s11_uni + s21_uni)/2, 'r--', 'Threshold');
xlabel('\phi_1');
title('Unipolar NRZ - Constellation (no noise)');
grid on;

subplot(1,2,2);
plot([s21_man s11_man], [0 0], 'bX', 'MarkerSize', 14, 'LineWidth', 2);
xline((s11_man + s21_man)/2, 'r--', 'Threshold');
xlabel('\phi_1');
title('Manchester - Constellation (no noise)');
grid on;

% ---- Generate random bits ----
bits = rand(1, Nb) > 0.5;

% Optimal thresholds
threshold_uni = (s11_uni + s21_uni) / 2;  % = 0.5
threshold_man = (s11_man + s21_man) / 2;  % = 0

% Pre-allocate BER arrays
BER_uni_sim   = zeros(1, numel(EbNo_dB));
BER_man_sim   = zeros(1, numel(EbNo_dB));
BER_uni_theory = zeros(1, numel(EbNo_dB));
BER_man_theory = zeros(1, numel(EbNo_dB));

% ---- Main simulation loop ----
for k = 1:numel(EbNo_dB)

    EbNo_lin = 10^(EbNo_dB(k) / 10);

    % Noise variance on the correlator output = No/2
    No_uni = Eb_uni / EbNo_lin;
    No_man = Eb_man / EbNo_lin;

    sigma_uni = sqrt(No_uni / 2);
    sigma_man = sqrt(No_man / 2);

    % Correlator outputs: constellation point + Gaussian noise
    x_uni = zeros(1, Nb);
    x_man = zeros(1, Nb);

    for i = 1:Nb
        if bits(i) == 1
            x_uni(i) = s11_uni + sigma_uni * randn();
            x_man(i) = s11_man + sigma_man * randn();
        else
            x_uni(i) = s21_uni + sigma_uni * randn();
            x_man(i) = s21_man + sigma_man * randn();
        end
    end

    % Plot received constellation at -10, 0, +10 dB
    if ismember(EbNo_dB(k), [-10 0 10])
        figure;
        subplot(1,2,1);
        plot(x_uni, zeros(1,Nb), '.', 'MarkerSize', 2);
        xline(threshold_uni, 'r--', 'Threshold');
        xlabel('\phi_1');
        title(['Unipolar NRZ Received  Eb/No = ' num2str(EbNo_dB(k)) ' dB']);
        grid on;

        subplot(1,2,2);
        plot(x_man, zeros(1,Nb), '.', 'MarkerSize', 2);
        xline(threshold_man, 'r--', 'Threshold');
        xlabel('\phi_1');
        title(['Manchester Received  Eb/No = ' num2str(EbNo_dB(k)) ' dB']);
        grid on;
    end

    % Decision
    decoded_uni = x_uni >= threshold_uni;
    decoded_man = x_man >= threshold_man;

    % Simulated BER
    BER_uni_sim(k) = sum(decoded_uni ~= bits) / Nb;
    BER_man_sim(k) = sum(decoded_man ~= bits) / Nb;

    % Unipolar NRZ:  s11=1, s21=0  ->  Pe = (1/2)*erfc( sqrt(Eb/No / 2) )
    % Manchester:    s11=1, s21=-1 ->  Pe = (1/2)*erfc( sqrt(Eb/No) )

    BER_uni_theory(k) = 0.5 * erfc(sqrt(EbNo_lin / 2));
    BER_man_theory(k) = 0.5 * erfc(sqrt(EbNo_lin));
end

% ---- BER plot ----
figure;
semilogy(EbNo_dB, BER_uni_sim,   'ro-',  'LineWidth', 1.5, 'DisplayName', 'Unipolar NRZ (Simulated)');
hold on;
semilogy(EbNo_dB, BER_uni_theory, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Unipolar NRZ (Theoretical)');
semilogy(EbNo_dB, BER_man_sim,   'bs-',  'LineWidth', 1.5, 'DisplayName', 'Manchester (Simulated)');
semilogy(EbNo_dB, BER_man_theory, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Manchester (Theoretical)');
grid on;
xlabel('E_b/N_o (dB)');
ylabel('Bit Error Rate');
title('BER vs Eb/No - Unipolar NRZ and Manchester');
legend('Location', 'southwest');
hold off;
