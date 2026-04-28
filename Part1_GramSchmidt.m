M = input('enter the number of signals: ');
Length = input('enter the length of the largest signal: ');
fprintf('Enter the signals as in arow\n'); 
Result = zeros(M, Length);
t = 0:Length-1; 
for i = 1: M
    Result(i,:)=input("");
end
%%% Gram–Schmidt method %%% 
S = Result';   
[length, M] = size(S);  
Q = zeros(length, M);
NumOfN=0;
for i = 1:M
    g = S(:,i);
    for j = 1:NumOfN
       Sij = sum(Q(:,j) .* g); %% s21 , s31 ,s32..
        g = g - (Sij * Q(:,j)); %%g2,g3.. 
    end
    energy_g = sum(g.^2);
    if energy_g > 1e-10   % <-- guard against zero vectors 
        NumOfN = NumOfN + 1; 
        Q(:, NumOfN) = g / sqrt(energy_g);
    end
end
N = Q(:, 1:NumOfN)';   
display(N);

%%% Plot the M input signals %%%
Coords = zeros(M, NumOfN);
for i = 1:M
    for j = 1:NumOfN
        Coords(i,j) = sum(Result(i,:)' .* Q(:,j)); 
    end
end
t_pulse = [0:Length-1; 1:Length];
t_pulse = t_pulse(:);

figure;
for i = 1: M
    subplot(M,1,i);
    s_plot = [Result(i,:); Result(i,:)];
    s_plot = s_plot(:);
    plot(t_pulse, s_plot, 'LineWidth', 1.5); 
    title(['signal ', num2str(i)]);
    grid on;
    xlim([0 Length]); 
end 

%%% plot the N basis functions %%%
figure;
for i = 1: NumOfN
    subplot(NumOfN,1,i);
    n_plot = [N(i,:); N(i,:)];
    n_plot = n_plot(:);
    plot(t_pulse, n_plot, 'LineWidth', 1.5); 
    title(['basis function ', num2str(i)]);
    grid on;
    xlim([0 Length]);
end 

%%% Plot the constellation diagram of the signals %%%
%In Case of 2D plot
if NumOfN <= 2
 figure;
 plot(Coords(:,1), Coords(:,2), 'X', 'MarkerSize', 10, 'LineWidth', 2);
 grid on;
 xlabel('\phi_1');
 ylabel('\phi_2');
 title('Constellation Diagram (2D)');
%In Case of 3D plot
else
 figure;
 plot3(Coords(:,1), Coords(:,2), Coords(:,3), 'X', 'MarkerSize', 10, 'LineWidth', 2);
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
E = sum(Coords.^2, 2);
fprintf('\nEnergy of signal:\n');
for i = 1:M
    fprintf('Signal %d energy = %.2f\n', i, E(i));
end
