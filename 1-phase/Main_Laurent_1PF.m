%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code by Juan S. Giraldo - jnse@ieee.org
% 02/2022 - EEMCS, University of Twente, NL
% 3Ph-Power flow using Laurent series expansion
%Please cite reference [A] when publishing results based on this algorithm.
%[A] Juan S. Giraldo, Oscar D. Montoya, Pedro P. Vergara, and Federico Milano, "A Fixed-Point Current Injection Power Flow forElectric Distribution Systems using Laurent Series", Electric Power Systems Research
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;

%% Import data


Vb = 11; % kV phase base voltage
Sbase = 1000.0; % kVA 1phase base apparent power


System_Data_Lines = readmatrix('Lines_136.csv');
System_Data_Nodes = readmatrix('Nodes_136.csv');
% System_Data_Nodes = readmatrix('Nodes_33_laurent.csv');
% System_Data_Lines = readmatrix('Lines_33_laurent.csv');

System_Data_Lines(:,3:4) = System_Data_Lines(:,3:4)*1/(Vb^2*1000/Sbase);%*1/(Vb^2*1000/Sbase);
System_Data_Lines(:,5) = System_Data_Lines(:,5)*(Vb^2*1000/Sbase);%*(Vb^2*1000/Sbase);


%% Make Ybus
%(Must be ordered - 1, 2, ... nb. and bus 1 must be the slack!)

RX = 1.0;

[Yss, Ysd, Yds, Ydd] = makeYbus_1f(System_Data_Lines, System_Data_Nodes, RX);

%% Power Flow


non_slack = find(System_Data_Nodes(:,2) ~= 1);
nb = length(System_Data_Nodes(:,1));

Vs = 1*exp(1j*0.0*pi/180);

%
num_points = 1;
mag_vector = 1.0;
angle_vector = 0.0;




num_load_models = 1;
alpha_P_vector = 1.0;
alpha_I_vector = 0.0;
alpha_Z_vector = 0.0;
lambda = 1.0;


K = 100;
tolerance = zeros(K,length(mag_vector));
Delta = zeros(length(mag_vector),1);
Time = zeros(length(mag_vector),1);
Iter = zeros(length(mag_vector),1);
min_V = zeros(length(mag_vector),1);


S = zeros(3*(nb-1),1);




DELTA = zeros(num_load_models,1);
ITER = zeros(num_load_models,1);
TIME = zeros(num_load_models,1);


alpha_P = repmat(alpha_P_vector(1) * ones(nb-1,1),1,1);
alpha_I = repmat(alpha_I_vector(1) * ones(nb-1,1),1,1);
alpha_Z = repmat(alpha_Z_vector(1) * ones(nb-1,1),1,1);

S = lambda(1)*(System_Data_Nodes(2:end,3) + 1j*System_Data_Nodes(2:end,4)) / Sbase;% + 0.00001/ Sbase;

mag_0 = mag_vector(1);
angle_0 = angle_vector(1);
Vo = mag_0*exp(1j*(angle_0*pi/180 + 0));

V_0 = repmat((Vo), (nb-1), 1);  % Flat start
V_0_0 = V_0;  % Flat start



k = 0;
epsilon = 1e-6;
tol = 100;

Vr = real(V_0);
Vi = imag(V_0);
B = sparse(diag(alpha_Z.*conj(S)) + Ydd);
C = sparse(Yds*Vs + alpha_I.*conj(S));

%% Iterative process
tic
while k <= K
    
    k = k + 1;
    
    
    
    A = sparse(diag(alpha_P.*conj(V_0.^(-2)).*conj(S)));
    D = sparse(2*(alpha_P.*conj(V_0.^(-1))).*conj(S));
    
    
    M11 = real(A) - real(B);
    M12 = imag(A) + imag(B);
    M21 = imag(A) - imag(B);
    M22 = -real(A) - real(B);
    N1 = real(C) + real(D);
    N2 = imag(C) + imag(D);
    
    M = [M11, M12; M21, M22];
    N = [N1;N2];
    V = M\N;
    
    V_0 = (V(1:1*(nb-1),1) + 1j*V(1*(nb-1)+1:end));
    V_0(isinf(V_0) | isnan(V_0)) = 4.0;
    
    
    Sd = -conj((Yds*Vs + Ydd*V_0).*conj(V_0)); %Sd([ncarga_a*3-2;ncarga_b*3-1;ncarga_c*3])
    tol = max(abs(Sd - (alpha_P + alpha_I.*V_0 + alpha_Z.*abs(V_0).^2).*S));
    %         tolerance(k,n) = tol;
    
    if tol <= epsilon
        break;
    end
    
end
time = toc;

Ss = conj(diag(conj(Vs))*(Yss*Vs + Ysd*V_0));
Delta(1,1) = 1/(max(abs(V_0_0))/min(abs(V_0_0-V_0)));
Iter(1,1) = k;
Time(1,1) = time*1000; % miliseconds
min_V(1,1) = min(abs(V_0(abs(V_0)>0.0)));

DELTA(1) = mean(Delta);
ITER(1) = mean(Iter);
TIME(1) = mean(Time);
%% Results

% Results = [DELTA,ITER,TIME];
Results = [Delta,Iter,Time];

% writematrix(Results,'Results_DS_start.csv')




