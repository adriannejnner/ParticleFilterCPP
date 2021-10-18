function variance = bisection(n,t,threshold,weight,Control_tumour_growth,PFTumourVolume)
    % start bisection for search variance
    init_variance_low = 0;
    init_variance_upper = 2500;
    init_variance = 2500;
    ESS_t_temp = n;
    tol = 1;
    weight_t_temp = zeros(1,n);
    
    while tol > 0.001
    %for j = 1:30
       init_variance = (init_variance_upper + init_variance_low)/2;
        
        for i = 1:n
           % update the weight
           weight_t_temp(1,i) = log(weight(t-1,i)) + log_normpdf(Control_tumour_growth(t),PFTumourVolume(t,i),init_variance);
       end
       
       weight_t_temp(1,:) = exp(weight_t_temp(1,:) - logsumexp(weight_t_temp(1,:),2));
       ESS_t_temp = 1/sum(weight_t_temp(1,:).^2);
       ESS_t_temp
       tol = abs(ESS_t_temp - n*threshold);
       
       % bisection searching for variance
       
       if  ESS_t_temp < n*threshold
           init_variance_low = init_variance;
       else
           init_variance_upper = init_variance;
       end
       
    end
    
    variance = init_variance;
end