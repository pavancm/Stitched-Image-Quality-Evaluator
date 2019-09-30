function [pi_l, mu, sigma, l, converged] = fit_gmm(data,K,Covtype,num_iter) 
%GMM
%Covtype = 1 for diagonal, 2 for full covariance
warning off;

converged = false;

D = size(data,2);
N = size(data,1);
Pl_x = zeros(N,K);
Px_l = zeros(N,K);

%k means
[idx,mu] = kmeans(data,K,'MaxIter',10);
sigma = repmat(diag(var(data)),[1,1,K]);
pi_l = ones(1,K,'like',data)/K ;

l = zeros(num_iter,1);
l_old = -Inf;

tol = 1e-6;
prob_thresh = 1e-8;
for iter = 1:num_iter
    %E step
    log_Px_l = logdensity(data,mu,sigma,pi_l,Covtype);
    
    
    if(isempty(log_Px_l))
        pi_l = [];mu = [];sigma = [];l = [];
        break;
    end
    
    max_Px_l = max(log_Px_l,[],2);
    Pl_x = exp(bsxfun(@minus, log_Px_l, max_Px_l));
    logpdf = log(sum(Pl_x,2)) + max_Px_l;
    %Q value
    ll = sum(logpdf);
    
    Pl_x = bsxfun(@rdivide,Pl_x,sum(Pl_x,2));
    %remove small posteriors
    Pl_x(Pl_x<prob_thresh) = 0;
    Pl_x = bsxfun(@rdivide,Pl_x,sum(Pl_x,2));
    
    %M Step
    sigma = zeros(D,D,K);
    Pl_x_sum = sum(Pl_x,1);
    for k = 1:K
%         mu(k,:) = Pl_x(:,k)' * data / Pl_x_sum(k);
        mu(k,:)=0;
        mean_sub = bsxfun(@minus,data,mu(k,:));
        le = bsxfun(@times,mean_sub,sqrt(Pl_x(:,k)));
        sigma(:,:,k) = ((le'*le/Pl_x_sum(k)));
    end
    pi_l = Pl_x_sum/sum(Pl_x_sum);
    
    l_diff = ll - l_old;
    if l_diff >=0 && l_diff < tol*abs(ll)
        converged = true;
        break; 
    end;
    l_old = ll;
    
end