function [Yss, Ysd, Yds, Ydd] = makeYbus_1f(branch_info, bus_info, RX)

nb = length(bus_info(:,1));  % number of buses
nl = length(branch_info(:,1));   % number of lines

sl = bus_info(bus_info(:,2) == 1);     % Slack node(s)

%
% %% for each branch, compute the elements of the branch admittance matrix where
% %%
% %%      | Is |   | Yss  Ysd |   | Vs |
% %%      |    | = |          | * |    |
% %%      |-Id |   | Yds  Ydd |   | Vd |
% %%


stat = branch_info(:, 6);%  %% ones at in-service branches
Ys = stat./((branch_info(:, 3)*RX + 1j * branch_info(:, 4))) ;% %% series admittance
Bc = stat.*branch_info(:,5);  % line charging susceptance
tap = stat.*branch_info(:, 7);  % default tap ratio = 1

Ytt = Ys + 1j * Bc / 2;
Yff = Ytt./ (tap);
Yft = - Ys./ (tap);
Ytf = Yft;


%% build connection matrices
f = branch_info(:, 1);  %% list of "from" buses
t = branch_info(:, 2);  %% list of "to" buses
%% connection matrix for line & from buses
Cf = sparse(1:nl,f,ones(nl,1),nl,nb);
%% connection matrix for line & to buses
Ct = sparse(1:nl,t,ones(nl,1),nl,nb);


%% build Yf and Yt such that Yf * V is the vector of complex branch currents injected
%% at each branch's "from" bus, and Yt is the same for the "to" bus end
i = repmat(1:nl,1,2);  %% double set of row indices

Yf = sparse(i',[f;t],[Yff;Yft]);
Yt = sparse(i',[f;t],[Ytf;Ytt]);

%% build Ybus
Ybus = Cf' * Yf + Ct' * Yt;        % Full Ybus

Yss = Ybus(1,1);
Ysd = Ybus(1,2:end);
Yds = Ybus(2:end,1);
Ydd = Ybus(2:end,2:end);



end

