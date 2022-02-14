%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code by Juan S. Giraldo - jnse@ieee.org - University of Twente, NL
% 3Ph-Power flow using Laurent series expansion
%Please cite reference [A] when publishing results based on this algorithm.
%[A] Juan S. Giraldo, Oscar D. Montoya, Pedro P. Vergara, and Federico Milano, "A Fixed-Point Current Injection Power Flow forElectric Distribution Systems using Laurent Series", Electric Power Systems Research
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

%% Import data


Vb = 4.16/sqrt(3); % kV phase base voltage
Sbase = 1000.0/3; % kVA 3phase base apparent power


System_Data_Nodes = readmatrix('Nodes_34_3Ph.csv');
System_Data_Lines = readmatrix('Lines_34_3Ph.csv');


if max(System_Data_Nodes(:,1) ~= length(System_Data_Nodes(:,1)))
    [System_Data_Nodes, System_Data_Lines] = renumerate(System_Data_Nodes, System_Data_Lines);
end

System_Data_Lines(:,3:end) = System_Data_Lines(:,3:end)*1/(Vb^2*1000/Sbase);

%% Make Ybus
%(Must be ordered - 1, 2, ... nb. and bus 1 must be the slack!)

[Yss, Ysd, Yds, Ydd, Zr] = Ybus_3F_func(System_Data_Lines);

%% Power Flow


non_slack = find(System_Data_Nodes(:,2) ~= 1);
nb = length(System_Data_Nodes(:,1));
ncarga_a = find(System_Data_Nodes(2:end,3));
ncarga_b = find(System_Data_Nodes(2:end,4));
ncarga_c = find(System_Data_Nodes(2:end,5));

a = 1*exp(1j*2*pi/3);
Vs = [1;a^2;a];

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


alpha_P = repmat(alpha_P_vector(1) * ones(nb-1,1),3,1);
alpha_I = repmat(alpha_I_vector(1) * ones(nb-1,1),3,1);
alpha_Z = repmat(alpha_Z_vector(1) * ones(nb-1,1),3,1);

SA = lambda(1)*(System_Data_Nodes(2:end,3) + 1j*System_Data_Nodes(2:end,6)) / Sbase + 0.00001/ Sbase;
SB = lambda(1)*(System_Data_Nodes(2:end,4) + 1j*System_Data_Nodes(2:end,7)) / Sbase + 0.00001/ Sbase;
SC = lambda(1)*(System_Data_Nodes(2:end,5) + 1j*System_Data_Nodes(2:end,8)) / Sbase + 0.00001/ Sbase;
S(3*System_Data_Nodes(2:nb,1)-5) = SA;
S(3*System_Data_Nodes(2:nb,1)-4) = SB;
S(3*System_Data_Nodes(2:nb,1)-3) = SC;

for n = 1:num_points
    
    
    mag_0 = mag_vector(n);
    angle_0 = angle_vector(n);
    Vo = [mag_0*exp(1j*(angle_0*pi/180 + 0)); mag_0*exp(1j*(angle_0*pi/180 - 2*pi/3));mag_0*exp(1j*(angle_0*pi/180 + 2*pi/3))];
    
    V_0 = repmat((Vo), (nb-1), 1);  % Flat start
    V_0_0 = V_0;  % Flat start
    
    T = repmat(([1;1;1]), (nb-1), 1);
    
    
    
    
    k = 0;
    epsilon = 1e-6;
    tol = 100;
    
    Vr = real(V_0);
    Vi = imag(V_0);
    B = sparse(diag(alpha_Z.*conj(S)) + Ydd);
    C = sparse(Yds*Vs + alpha_I.*conj(T).*conj(S));
    
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
        
        V_0 = (V(1:3*(nb-1),1) + 1j*V(3*(nb-1)+1:end));
        V_0(isinf(V_0) | isnan(V_0)) = 4.0;
        
        
        Sd = -conj((Yds*Vs + Ydd*V_0).*conj(V_0)); %Sd([ncarga_a*3-2;ncarga_b*3-1;ncarga_c*3])
        tol = max(abs(Sd([ncarga_a*3-2;ncarga_b*3-1;ncarga_c*3]) - (alpha_P([ncarga_a*3-2;ncarga_b*3-1;ncarga_c*3]) + (alpha_I([ncarga_a*3-2;ncarga_b*3-1;ncarga_c*3]).*(T([ncarga_a*3-2;ncarga_b*3-1;ncarga_c*3])).*(V_0([ncarga_a*3-2;ncarga_b*3-1;ncarga_c*3]))) + (alpha_Z([ncarga_a*3-2;ncarga_b*3-1;ncarga_c*3]).*abs(V_0([ncarga_a*3-2;ncarga_b*3-1;ncarga_c*3])).^2)).*S([ncarga_a*3-2;ncarga_b*3-1;ncarga_c*3])));
        %         tolerance(k,n) = tol;
        
        if tol <= epsilon
            break;
        end
        
    end
    time = toc;
    
    Ss = conj(diag(conj(Vs))*(Yss*Vs + Ysd*V_0));
    Delta(n,1) = 1/(max(abs(V_0_0))/min(abs(V_0_0-V_0)));
    Iter(n,1) = k;
    Time(n,1) = time*1000; % miliseconds
    min_V(n,1) = min(abs(V_0(V_0>0.0)));
end


DELTA(1) = mean(Delta);
ITER(1) = mean(Iter);
TIME(1) = mean(Time);
%% Results

Results = [Delta,Iter,Time];





