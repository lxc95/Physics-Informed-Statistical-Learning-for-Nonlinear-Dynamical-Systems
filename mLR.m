function [B_estimated,covariance_LinearReg]=mLR(FPCA_force,FPCA_dis)
num_pre = 35;
B_estimated = inv((FPCA_force(1:34,:)'*FPCA_force(1:34,:)))*FPCA_force(1:34,:)'*FPCA_dis(1:34,:);
covariance_LinearReg = 1/(34-1)*(FPCA_dis(1:34,:) - FPCA_force(1:34,:)*B_estimated)'*(FPCA_dis(1:34,:) - FPCA_force(1:34,:)*B_estimated);%note to use 34 training data
