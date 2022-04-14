clear all, close all, clc

run('vvc_setup_IEEE123s.m')
load vvc_srate_2_25.mat

Sphase_shaohui = 5.1834e+7/1000;

Qg_pred_cvar_qonly_sampling_reg = csvread('qg_pred_cvar_qonly_sampling_reg2.csv',0,0)/Sphase_shaohui;
Qg_pred_cvar_qonly_reg = csvread('qg_pred_cvar_qonly_reg2.csv',0,0)/Sphase_shaohui;

delete_time = csvread('delete_list.csv',0,0);

Q(:,delete_time) = [];
Qg(:,delete_time) = [];
P(:,delete_time) = [];
% Qg_pred_cvar_qonly_reg(:,delete_time) = [];

Tnew = T - length(delete_time);

Qg_cvar_qonly_sampling_reg = zeros(N,Tnew); Qg_cvar_qonly_sampling_reg(pv_idx,:) = Qg_pred_cvar_qonly_sampling_reg;
Qg_cvar_qonly_reg = zeros(N,Tnew); Qg_cvar_qonly_reg(pv_idx,:) = Qg_pred_cvar_qonly_reg;

v_opt = zeros(122,Tnew);
v_cvar = v_opt;
v_cvar_qonly = v_opt;

w = ones(122,1);
YLLi = inv(YLL);

Q_opt = Qg - Q;
Q_cvar_qonly_sampling = Qg_cvar_qonly_sampling_reg - Q;
Q_cvar_qonly = Qg_cvar_qonly_reg - Q;
for t = 1:Tnew
    if rem(t,100) == 0
       t
    end
    x_opt = [-P(:,t); Q_opt(:,t)];
    x_cvar = [-P(:,t); Q_cvar_qonly_sampling(:,t)];
    x_cvar_qonly = [-P(:,t); Q_cvar_qonly(:,t)];
    
    v_opt(:,t) = FPM(w, x_opt, YLLi);
    v_cvar(:,t) = FPM(w, x_cvar, YLLi);
    v_cvar_qonly(:,t) = FPM(w, x_cvar_qonly, YLLi);
end

vmag_opt = abs(v_opt);
vmag_cvar_qonly_sampling = abs(v_cvar);
vmag_cvar_qonly = abs(v_cvar_qonly);

dv_opt = vmag_opt - ones(N,Tnew);
dv_cvar_qonly_sampling = vmag_cvar_qonly_sampling - ones(N,Tnew);
dv_cvar_qonly = vmag_cvar_qonly - ones(N,Tnew);

figure, plot(vmag_opt'), ylabel('|v| (pu)'), title('Optimal')
figure, plot(vmag_cvar_qonly_sampling'), ylabel('|v| (pu)'), title('CVaR')
figure, plot(vmag_cvar_qonly'), ylabel('|v| (pu)'), title('CVaR - qonly')

figure
subplot(1,3,1), hist(vec(dv_opt),100), title('Optimal'), xlim([-0.06 0])
subplot(1,3,2), hist(vec(dv_cvar_qonly_sampling),100), title('CVaR - qonly sampling')
subplot(1,3,3), hist(vec(dv_cvar_qonly),100), title('CVaR - qonly')

filename = 'dv_opt2.txt';
csvwrite(filename,dv_opt);

filename = 'dv_cvar_qonly_sampling2.txt';
csvwrite(filename,dv_cvar_qonly_sampling);

filename = 'dv_cvar_qonly2.txt';
csvwrite(filename,dv_cvar_qonly);


%%
function v = FPM(w0,xwye,YLLi)
    jay = sqrt(-1);
    err = 10;
    v_minus = w0;
    while err > 1e-12
        cV = diag(conj(v_minus));
        cV_inv = inv(cV);
        term1 = YLLi*cV_inv;
        term2 = -jay*YLLi*cV_inv;
        Mwye = [term1, term2];
        v_plus = Mwye*xwye+w0;
        err = sum(abs(v_minus-v_plus));
        v_minus = v_plus;
    end
    v = v_minus;
end