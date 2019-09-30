function   [log_lh,mahalaD,re_init]=logdensity(X, mu, Sigma, p, CovType)

    re_init = 0;
    log_prior = log(p);
    [n,d]=size(X);
    k=size(mu,1);
    log_lh = zeros(n,k);
    mahalaD = zeros(n,k);
    logDetSigma = -Inf;
    for j = 1:k 
        if CovType == 2 %full covariacne
            % compute the log determinant of covariance
            [L,f] = chol(Sigma(:,:,j) );
            diagL = diag(L);
            logDetSigma = 2*sum(log(diagL));
        else %diagonal covariance
            L = sqrt(diag(Sigma(:,:,j))); % a vector
            logDetSigma = sum( log(diag(Sigma(:,:,j))) );
     
        end
        
        Xcentered = bsxfun(@minus, X, mu(j,:));
       
        if(isempty(L))
            log_lh = [];
            break
        end
        
        if(size(Xcentered,2)~=size(L,1))
            log_lh = [];
            re_init = 1; %initialize once more
            break;
        end
        
        if CovType == 2
            xRinv = Xcentered /L ;
        else
            xRinv = bsxfun(@rdivide,Xcentered ,L');
        end
        mahalaD(:,j) = sum(xRinv.^2, 2);

        log_lh(:,j) = -0.5 * mahalaD(:,j) +...
            (-0.5 *logDetSigma + log_prior(j)) - d*log(2*pi)/2;
        %get the loglikelihood for each point with each component
        %log_lh is a N by K matrix, log_lh(i,j) is log \alpha_j(x_i|\theta_j)
    end

   
