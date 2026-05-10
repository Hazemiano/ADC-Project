
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
              
                S(i, :) = user_input(1:Length);
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
            % Rectangular Pulse
            rect = floor(Length/4); 
            rect_wave = zeros(1, Length);
            rect_wave(abs(t_c) <= rect) = 1;
            S(i,:) = rect_wave;
            
        case 5
            % Triangular Pulse
            tri = 1 - abs(t_c) /(Length/4);
            tri(tri < 0) = 0; 
            S(i,:) =tri;
            
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
% Parameters
Nb = 100000;     % Number of bits
Tb = 1;
fs = 100;
ts = 1/fs;

% Generate ONE random bitstream
b = rand(1, Nb) > 0.5;

UPNRZ = repelem(b, fs);

t = 0:ts:NbTb-ts;

figure;
plot(t, UPNRZ);
axis([0 10 -0.5 1.5]);   % Show only first 10 bits (important!)
title('Unipolar NRZ');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% manchester
Manchester = zeros(1, Nbfs);

for i = 1:Nb
    idx = (i-1)fs + 1 : ifs;

    if b(i) == 1
        Manchester(idx) = [ones(1, fs/2) -ones(1, fs/2)];
    else
        Manchester(idx) = [-ones(1, fs/2) ones(1, fs/2)];
    end
end

figure;
plot(t, Manchester);
axis([0 10 -1.5 1.5]);   % Show first 10 bits only
title('Manchester Code');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
   %

SPOWER = 1;

SNR = [-10:2:10]; %in Db
snr = 10.^(0.1.SNR);

%for I = 1:length(snr)
   % noise = 1 / sqrt(2) (randn(1, 10000) + 1i * randn(1, 10000));
    %u =Manchester+noise .* snr(I);
%end
SPOWER=1;
for k = 1:length(SNR)

    N0=SPOWER/snr(k);
    sigma=sqrt(N0/2);
    noise=sigma/sqrt(2)(randn(size(Manchester))+1irandn(size(Manchester)));
    Part2Result= Manchester+noise;
    subplot(3,4,k);   

    scatter(real(Part2Result), imag(Part2Result), 10, 'filled');
    grid on;
    axis equal;

    xlabel('In-phase');
    ylabel('Quadrature');

    title(['Manchester, SNR = ', num2str(SNR(k)), ' dB']);
end

SPOWER=1;

for k = 1:length(SNR)

    N0=SPOWER/snr(k);
    sigma=sqrt(N0/2);
    noise=sigma/sqrt(2)(randn(size(UPNRZ))+1irandn(size(UPNRZ)));
    Part3Result= UPNRZ+noise;
    subplot(3,4,k);   

    scatter(real(Part2Result), imag(Part2Result), 10, 'filled');
    grid on;
    axis equal;

    xlabel('In-phase');
    ylabel('Quadrature');

    title(['Manchester, SNR = ', num2str(SNR(k)), ' dB']);
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

