function [out,Var_co] = mGPR_prediction_diagnal(w,x,y,k,xt)
n = size(x,1); d = size(y,2);
para_sigma2 = exp(w(1));
para_diagonal_Lambda = exp(w(end-d+1:end));%enforce the diagonal elements to be positive
% para_Lambda = [para_diagonal_Lambda(:)];
% para_kernel = exp(w(2:end-d));
[L_Lambda,Lambda] = vec2mat_diag([para_diagonal_Lambda(:);zeros((d-1)*d/2,1)],d);
K11 = k(w(2:end-d),x,x)+eye(n)*para_sigma2;
[L_kernel,p] = chol(K11,'lower');

vv = L_kernel\y;            vv_Lambda = L_Lambda\y';
alpha = L_kernel'\vv;       alpha_Lambda = L_Lambda'\vv_Lambda;

% e = d*n/2*log(2*pi)+ d*sum(log(diag(L_kernel)))+ n*sum(log(diag(L_Lambda)))...
%     + trace(0.5*alpha*alpha_Lambda);
%
K22 = k(w(2 : end-d),xt,xt);
K21 = k(w(2 : end-d),xt,x);

% Solve the mean
Eft = K21*alpha;

% Solve the variance
v = L_kernel\K21';
v_temp = L_kernel'\v;
v_new = K21*v_temp;
%----------------------------------------------------------------------
Varft_total = kron(diag(K22) - diag(v_new),diag(Lambda));
Var_co = kron(K22 - v_new,Lambda);
Varft = reshape(Varft_total,d,size(xt,1))';
out = [Eft Varft p];


% % Solve the variance
% v = L_kernel\K21';
% %----------------------------------------------------------------------
% Varft_total = kron(diag(K22) - sum(v.^2,1)',diag(Lambda));
% Varft = reshape(Varft_total,d,size(xt,1))';
% out = [Eft Varft*1e7 p];
%----------------------------------------------------------------------
end