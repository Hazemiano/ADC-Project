clear all;
M = input('Enter the number of signals: ');
Length = input('Enter the length of the largest signal: '); 
S = zeros(M, Length);
t = 0:Length-1; 
t_c = t - floor(Length/2); 

%%% Input Signals Menu %%%
for i = 1:M
    fprintf('\n type of Signal %d\n', i);
    fprintf('1. Custom Array (Any length, e.g., [1 2 -3])\n');
    fprintf('2. Sin Wave\n');
    fprintf('3. Cos Wave\n');
    fprintf('4. Rectangular Pulse (Rect)\n');
    fprintf('5. Triangular Pulse\n');
    fprintf('6. Sinc\n');
    
    choice = input('Select the signal type (1-6): ');
    
    switch choice
        case 1
            % Zero-Padding
            fprintf('Enter the signal as a row vector:\n');
            user_input = input('');
            user_input = user_input(:)';
            len_temp = length(user_input);
            
            if len_temp <= Length
                 % Zero-Padding
                S(i, 1:len_temp) = user_input;
            else
              
                S(i,:) = user_input(1:Length);
                disp('Warning: Input signal was longer than the max length you are selected');
            end
            
     case 2
    
    period=input('Enter period');
    
 
    
    temp=zeros(1,Length);  
    n =0:period-1; 
    temp(1:period) =sin(2*pi*n/period);  
    S(i,:) =temp;

case 3
    
    period=input('Enter period ');
  
    
    temp =zeros(1,Length);
    n= 0:period-1;
    temp(1:period) =cos(2*pi*n/period);
    S(i,:) =temp;
    
            
       case 4

        dur = input('Enter pulse duration (samples): ');

         if dur < 1
             disp('Invalid duration! Using 1.');
             dur = 1;
         end

        rect_wave = zeros(1,Length);
        rect_wave(1:dur) =1;   
        S(i,:) =rect_wave;
      case 5

    period = input('nter period');

    
    temp = zeros(1, Length);

    n=0:period-1;

    half=floor(period/2);

    tri_one_period = zeros(1, period);

    tri_one_period(1:half)=linspace(0, 1, half);
    tri_one_period(half+1:end)=linspace(1, 0, period-half);
    temp(1:period)=tri_one_period;
    S(i,:)=temp;
            
        case 6
            % Sinc
            sinc_wave = zeros(1, Length);
            for k = 1:Length
                x = t_c(k) / (Length/8); 
                if x == 0
                    sinc_wave(k) = 1; 
                else
                    sinc_wave(k) = sin(pi * x) / (pi * x);
                end
            end
            S(i,:) = sinc_wave;
            
        otherwise
            disp('Invalid choice! Defaulting to zeros.');
    end
end

%%% Gram–Schmidt method %%% 
S_mat = S';  
[len, M_S] = size(S_mat);  
Q = zeros(len, M_S);%% empty array
counter=0;

for i = 1:M_S
    g = S_mat(:,i);
    for j = 1:counter
       Sij = sum(Q(:,j) .* g); %% s21 , s31 ,s32..
        g = g - (Sij * Q(:,j)); %%g2,g3.. 
    end
    energy_g = sum(g.^2);
    if energy_g > 1e-10   %  guard against zero vectors 
        counter = counter + 1; 
        Q(:, counter) = g / sqrt(energy_g);
    end
end

N = Q(:, 1:counter)';   
disp('The Basis Functions are:');
display(N);

%%% Plot the M input signals %%%
axisxyz = zeros(M, counter);
for i = 1:M
    for j = 1:counter
        axisxyz(i,j) = sum(S(i,:)' .* Q(:,j)); 
    end
end

t_pulse = [0:Length-1; 1:Length];
t_pulse = t_pulse(:);

figure;
for i = 1: M
    subplot(M,1,i);
    s_plot = [S(i,:); S(i,:)];
    s_plot = s_plot(:);
    plot(t_pulse, s_plot, 'LineWidth', 1.5); 
    title(['Signal ', num2str(i)]);
    grid on;
end 

%%% plot the N basis functions %%%
figure;
for i = 1: counter
    subplot(counter,1,i);
    n_plot = [N(i,:); N(i,:)];
    n_plot = n_plot(:);
    plot(t_pulse, n_plot, 'LineWidth', 1.5); 
    title(['Basis Function ', num2str(i)]);
    grid on;
end 

%%% Plot the constellation diagram of the signals %%%
if counter == 1
    % In Case of 1D plot (1 Basis Function)
    figure;
    plot(axisxyz(:,1), zeros(M,1), 'X', 'MarkerSize', 10, 'LineWidth', 2);
    grid on;
    xlabel('\phi_1');
    title('Constellation Diagram (1D)');
    
elseif counter == 2
    % In Case of 2D plot (2 Basis Functions)
    figure;
    plot(axisxyz(:,1), axisxyz(:,2), 'X', 'MarkerSize', 10, 'LineWidth', 2);
    grid on;
    xlabel('\phi_1');
    ylabel('\phi_2');
    title('Constellation Diagram (2D)');
    
else
    % In Case of 3D plot (3 Basis Functions or more)
    figure;
    plot3(axisxyz(:,1), axisxyz(:,2), axisxyz(:,3), 'X', 'MarkerSize', 10, 'LineWidth', 2);
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
E = sum(axisxyz.^2, 2);
fprintf('\nEnergy of signal:\n');
for i = 1:M
    fprintf('Signal %d energy = %.2f\n', i, E(i));
end


%%% Unipolar NRZ encoding %%%

%%% Part 2 - BER %%%
Nb = 100000;
Tb = 1;
fs = 100;
ts = 1/fs;

b = rand(1, Nb) > 0.5;

% build unipolar NRZ
UPNRZ = zeros(1, Nb*fs);
for i = 1:Nb
    UPNRZ((i-1)*fs+1 : i*fs) = b(i);
end

% build manchester
Manchester = zeros(1, Nb*fs);
for i = 1:Nb
    idx = (i-1)*fs+1 : i*fs;
    if b(i) == 1
        Manchester(idx) = [ones(1,fs/2) -ones(1,fs/2)];
    else
        Manchester(idx) = [-ones(1,fs/2) ones(1,fs/2)];
    end
end

% plot waveforms
t = 0:ts:Nb*Tb-ts;
figure;
subplot(2,1,1);
plot(t, UPNRZ);
axis([0 10 -0.5 1.5]);
title('Unipolar NRZ');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(t, Manchester);
axis([0 10 -1.5 1.5]);
title('Manchester Code');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% basis functions
phi1_NRZ = ones(1,fs) / sqrt(Tb);
phi2_NRZ = [ones(1,fs/2) -ones(1,fs/2)] / sqrt(Tb);
phi1_MAN = [ones(1,fs/2) -ones(1,fs/2)] / sqrt(Tb);
phi2_MAN = ones(1,fs) / sqrt(Tb);

% Eb for each scheme
Eb_NRZ = 1 * Tb;
Eb_MAN = 1 * Tb;

SNR_dB = -10:2:10;
snr = 10.^(SNR_dB/10);

BER_UPNRZ = zeros(1, numel(SNR_dB));
BER_Manc  = zeros(1, numel(SNR_dB));

x1_NRZ = zeros(1, Nb);
x2_NRZ = zeros(1, Nb);
x1_MAN = zeros(1, Nb);
x2_MAN = zeros(1, Nb);

for k = 1:numel(SNR_dB)

    % noise std for each scheme
    No_NRZ    = Eb_NRZ / snr(k);
    sigma_NRZ = sqrt(No_NRZ * fs / 2);

    No_MAN    = Eb_MAN / snr(k);
    sigma_MAN = sqrt(No_MAN * fs / 2);

    % unipolar NRZ
    received = UPNRZ + sigma_NRZ * randn(size(UPNRZ));
    decoded  = zeros(1, Nb);
    for i = 1:Nb
        idx = (i-1)*fs+1 : i*fs;
        r = sum(received(idx) .* phi1_NRZ) * ts;
        if r >= sqrt(Tb)/2
            decoded(i) = 1;
        end
        if k == numel(SNR_dB)
            x1_NRZ(i) = sum(received(idx) .* phi1_NRZ) * ts;
            x2_NRZ(i) = sum(received(idx) .* phi2_NRZ) * ts;
        end
    end
    BER_UPNRZ(k) = sum(decoded ~= b) / Nb;

    % manchester
    received = Manchester + sigma_MAN * randn(size(Manchester));
    decoded  = zeros(1, Nb);
    for i = 1:Nb
        idx = (i-1)*fs+1 : i*fs;
        r = sum(received(idx) .* phi1_MAN) * ts;
        if r >= 0
            decoded(i) = 1;
        end
        if k == numel(SNR_dB)
            x1_MAN(i) = sum(received(idx) .* phi1_MAN) * ts;
            x2_MAN(i) = sum(received(idx) .* phi2_MAN) * ts;
        end
    end
    BER_Manc(k) = sum(decoded ~= b) / Nb;

end

% theoretical BER
BER_UPNRZ_theory = 0.5 * erfc(sqrt(snr/4));
BER_Manc_theory  = 0.5 * erfc(sqrt(snr));

% plot BER
figure;
semilogy(SNR_dB, BER_UPNRZ,        'b-o', 'DisplayName', 'Unipolar NRZ (Simulated)');
hold on;
semilogy(SNR_dB, BER_UPNRZ_theory, 'b--', 'DisplayName', 'Unipolar NRZ (Theoretical)');
semilogy(SNR_dB, BER_Manc,         'r-o', 'DisplayName', 'Manchester (Simulated)');
semilogy(SNR_dB, BER_Manc_theory,  'r--', 'DisplayName', 'Manchester (Theoretical)');
hold off;
legend('Location', 'southwest');
title('BER vs Eb/No');
xlabel('Eb/No (dB)');
ylabel('BER');
grid on;

% constellation with noise at Eb/No = 10 dB
figure;
subplot(1,2,1);
plot(x1_NRZ(b==1), x2_NRZ(b==1), '.b', 'MarkerSize', 3, 'DisplayName', 'bit=1');
hold on;
plot(x1_NRZ(b==0), x2_NRZ(b==0), '.r', 'MarkerSize', 3, 'DisplayName', 'bit=0');
hold off;
legend;
grid on;
axis equal;
xlabel('\phi_1');
ylabel('\phi_2');
title('Unipolar NRZ Constellation (with noise)');

subplot(1,2,2);
plot(x1_MAN(b==1), x2_MAN(b==1), '.b', 'MarkerSize', 3, 'DisplayName', 'bit=1');
hold on;
plot(x1_MAN(b==0), x2_MAN(b==0), '.r', 'MarkerSize', 3, 'DisplayName', 'bit=0');
hold off;
legend;
grid on;
axis equal;
xlabel('\phi_1');
ylabel('\phi_2');
title('Manchester Constellation (with noise)');
