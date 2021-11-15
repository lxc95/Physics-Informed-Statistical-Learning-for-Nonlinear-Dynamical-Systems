function [predict_FPCA_displacement_score,full_predict_displacement_filed,covariance_FPCA_displacement,full_predict_displacement_POD,covariance_POD_displacement]=Prediction(r,B_estimated,covariance_LinearReg,predict_FPCA_force_score,predict_FPCA_force_covariance,displacement_allcondition_mean,force_allcondition_mean,displacement_allcondition_centered,force_allcondition_centered,FPCA_dis_mean,FPCA_dis,FPCA_force_mean,FPCA_force,displacement_allcondition,force_allcondition,displacementmean,displacementFPCAbasis,forcemean,forceFPCAbasis)
nn=20;mm=12;% the truncation size in FPCA, POD 30
predict_FPCA_displacement_score = predict_FPCA_force_score*B_estimated;
full_predict_displacement_filed = displacement_allcondition_mean' + displacementFPCAbasis(:,1:nn)*predict_FPCA_displacement_score';
covariance_FPCA_displacement = B_estimated'*predict_FPCA_force_covariance*B_estimated + covariance_LinearReg;
%convert concatenated data to seperated data
full_predict_displacement_POD = zeros(r,400);
for i = 1:r
    full_predict_displacement_POD(i,:)=full_predict_displacement_filed((400*(i-1)+1):(400*(i-1)+400));
end
%now firstly just consider
%covariance(full_predict_displacement_POD(:,t),full_predict_displacement_POD(:,t))!!!!
%which is at the same time
t_step_predict = 400;
basis_special = displacementFPCAbasis(t_step_predict:400:(r*400),1:nn);
covariance_POD_displacement = basis_special*covariance_FPCA_displacement*basis_special';
