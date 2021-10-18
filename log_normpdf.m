function result = log_normpdf(y,x,sigma)
% y is scalar
% x is scalar
% sigma = 1

result = log(1/(sigma*sqrt(2*pi))) - 0.5*((y-x)/sigma)^2;

end