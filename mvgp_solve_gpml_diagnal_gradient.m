function [e,eg] = mvgp_solve_gpml_diagnal_gradient(w,x,y,k,xt)
n = size(x,1); d = size(y,2);
para_sigma2 = exp(w(1));
para_diagonal_Lambda = exp(w(end-d+1:end));%enforce the diagonal elements to be positive
para_Lambda = para_diagonal_Lambda(:);
para_kernel = exp(w(2:end-d));
[L_Lambda,Lambda] = vec2mat_diag([para_diagonal_Lambda(:);zeros((d-1)*d/2,1)],d);
K11 = k(w(2:end-d),x,x)+eye(n)*para_sigma2;
[L_kernel,p] = chol(K11,'lower');
if p>0
    
    % Add jitter and try again
    jitter = 1e-9*diag(rand(n,1));
    [L_kernel,p] = chol(K11+jitter,'lower');
    
    % Still no luck
    if p>0, e=nan; return; end
    
end
vv = L_kernel\y;            vv_Lambda = L_Lambda\y';
alpha = L_kernel'\vv;       alpha_Lambda = L_Lambda'\vv_Lambda;

e = d*n/2*log(2*pi)+ d*sum(log(diag(L_kernel)))+ n*sum(log(diag(L_Lambda)))...
    + trace(0.5*alpha*alpha_Lambda);

num_Lambda = numel(para_Lambda); num_kernel = numel(para_kernel);
dk = xt;
% Allocate space
eg = zeros(1+num_Lambda+num_kernel,1);
% Derivative w.r.t. sigma2
invK = L_kernel'\(L_kernel\eye(n));
invLambda = L_Lambda'\(L_Lambda\eye(d));
eg(1) = d*0.5*trace(invK) -  0.5*trace(alpha*invLambda*alpha');
% The params in the kernel
for j= 1:num_kernel
    %             dK = dk{j}(bsxfun(@minus,x(:),x(:)'),para_kernel);
    dK = dk{j}(w(2 : end-d),x,x);
    
    eg(j+1) = d*0.5*trace(invK*dK) - 0.5*trace(alpha*invLambda*alpha'*dK);
end
% Return derivatives
%eg = eg.*exp(w);%how to explain the exp()
eg = eg;

% The parameters in the Omega(diagonal and non-diagonal)
dLambda = zeros(d,d);
for i = 1:d
    for j = 1:i
        E = zeros(d,d);
        if i==j
            E(i,j) = para_diagonal_Lambda(i);
        else
            E(i,j) = 1;
        end
        dLambda(i,j) = 0.5*n*trace(invLambda*(E*L_Lambda'+L_Lambda*E)) - 0.5*trace(...
            alpha_Lambda*invK*alpha_Lambda'*(E*L_Lambda'+L_Lambda*E));
    end
end

% reshape dOmega in a vector,eg_Omega
eg_Lambda = cell(d,1);
for l= 1:d
    eg_Lambda{l} = diag(dLambda,1-l);
end

eg_Lambda = cell2mat(eg_Lambda);

% throw all the parameter in eg_Omega into eg
for i = 1 :num_Lambda
    eg(i+1+num_kernel) = eg_Lambda(i);
end

end