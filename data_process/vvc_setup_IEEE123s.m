%% Distribution Modeling Setup File
load IEEE123_system_data3.mat 

% Variables:
% Pload Qload Pslack Qslack Ppv Ybus Y00 Y0L YL0 YLL V Vpu nodes load_nodes pv_nodes Sphase Vphase Iphase Vbase Sbase Ybase

% Load and voltage data setup
Vmag = abs(V(2:end,:));
Vmag_pu = abs(Vpu(2:end,:));
Pg = Ppv/Sphase;
Pload_pu = Pload/Sphase;
Qload_pu = Qload/Sphase;
Pslack_pu = Pslack/Sphase;
Qslack_pu = Qslack/Sphase;
v0 = Vphase;
v0_pu = 1;

% Incidence matrix setup
[N,T] = size(Pload);
[Mtilde, y] = nodeEdgeIncidence(Ybus);
[temp, L] = size(Mtilde);
M = Mtilde(2:end,:);
invM = inv(M);

% Line impedance
z = 1./y;
r = real(z);
x = imag(z);

% r to x ratio
z_ratio = r./x;

% Linear approximation
P = Pload_pu;
Q = Qload_pu;
R = invM'*diag(r)*invM;
X = invM'*diag(x)*invM;

% Index of nodes with loads
slack_idx = find(nodes == 150);
temp_node = nodes;
temp_node(slack_idx) = [];

noload_idx = find(P(:,1) == 0);
load_idx = find(P(:,1) ~= 0);

pv_idx = find(ismember(temp_node,pv_nodes));
no_pv_idx = find(~ismember(temp_node,pv_nodes));