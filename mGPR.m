function [predict_FPCA_force_score,predict_FPCA_force_covariance,displacement_allcondition_mean,force_allcondition_mean,displacement_allcondition_centered,force_allcondition_centered,FPCA_dis_mean,FPCA_dis,FPCA_force_mean,FPCA_force]=mGPR(displacement_allcondition,force_allcondition,displacementmean,displacementFPCAbasis,forcemean,forceFPCAbasis)
load('all_UAVpose35_210612');
displacement_allcondition_mean = displacementmean';
force_allcondition_mean = forcemean';%this mean is important for following steps, centering
%plot(displacement_allcondition_mean);
displacement_allcondition_centered = displacement_allcondition - displacement_allcondition_mean;
force_allcondition_centered = force_allcondition - force_allcondition_mean;
nn=20;mm=12;% the truncation size in FPCA, POD 30
FPCA_dis_mean = displacement_allcondition_centered;
FPCA_dis = FPCA_dis_mean*displacementFPCAbasis(:,1:nn);
FPCA_force_mean = force_allcondition_centered;
FPCA_force = FPCA_force_mean*forceFPCAbasis(:,1:mm);
% datset 4:UAV
xtr=all_UAVpose35(1:34,:);
ytr=FPCA_force(1:34,:);
sn=0.0001;
xt=all_UAVpose35(35,:);
yt=FPCA_force(35,:);
%%
covfunc = @covSEard;
para_init = @SE_init;
n_parameter = para_init(xtr, ytr, 1);
dk = cell(n_parameter,1);
k = @(hyp,x,z) covfunc(hyp,x,z);
for i = 1:n_parameter
    dk{i} = @(hyp,x,z) covfunc(hyp,x,z,i);
end
%------------------------------------------------------------------
options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'display','off');
nlml_gp= Inf;
numinit = 100;
funcGP = @mvgp_solve_gpml_diagnal_gradient;
for j=1:numinit                 % multiple output regression
    kernel_gp = para_init(xtr, ytr);
    %------------------------------------------------------------------
    [diag_Omega_gp,non_diag_Omega_gp] = Omega_init(xtr, ytr);
    %------------------------------------------------------------------
    % param_gp = [log([sn; kernel_gp]);diag_Omega_gp];
    param_gp = [log([sn; kernel_gp]);log(diag_Omega_gp)];
    %------------------------------------------------------------------
    % Optimization
    [~,nlml_gp_new] = fminunc(@(w) funcGP(w,xtr,ytr,k,dk),param_gp,options);
    if (nlml_gp_new < nlml_gp)
        param_gp_final = param_gp;
        nlml_gp = nlml_gp_new;
    end
end
[w_gp_final,nlml_gp_final] = fminunc(@(w) funcGP(w,xtr,ytr,k,dk),param_gp_final,options);
%------------------------------------------------------------------
[prediction,covariance_pre] = mGPR_prediction_diagnal(w_gp_final,xtr,ytr,k,xt);
%------------------------------------------------------------------
predict_FPCA_force_score = prediction(1:size(yt,2));
predict_FPCA_force_covariance = covariance_pre;
