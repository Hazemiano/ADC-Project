
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
