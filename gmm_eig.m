%Fit GMM and obtain eigen values
function eig_val = gmm_eig(X,K,Covtype,num_iter)

rng(1);
[pi,mu,sigma,~,converged] = fit_gmm(X,K,Covtype,num_iter);

% [sigma,alpha,beta] = FitMGGD(X');

if(~isempty(sigma) && ~isempty(pi))

    V = reshape(sigma,[2,8]);
    W = repmat(pi,[4,1]);W = reshape(W,[2,8]);
    c_mt = W.*V;
    c_mt = [sum(c_mt(:,1:2:end),2),sum(c_mt(:,2:2:end),2)];
    [~,D] = eig(c_mt);
    eig_val = [D(1,1),D(2,2)];
else
    eig_val = [0,0];
end
% eig_val = [0,0];