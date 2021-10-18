function C = weightedcov(Y, w)
ctrl = isvector(w) & isreal(w) & ~any(isnan(w)) & ~any(isinf(w)) & all(w > 0);
% if ctrl
w = w(:) / sum(w);                                                              % w is column vector
% else
%   error('Check w: it needs be a vector of real positive numbers with no infinite or nan values!')
% end
% ctrl = isreal(Y) & ~any(isnan(Y)) & ~any(isinf(Y)) & (size(size(Y), 2) == 2);
% if ~ctrl
%   error('Check Y: it needs be a 2D matrix of real numbers with no infinite or nan values!')
% end
% ctrl = length(w) == size(Y, 1);
% if ~ctrl
%   error('size(Y, 1) has to be equal to length(w)!')
% end
[T, N] = size(Y);                                                                 
C = Y - repmat(w' * Y, T, 1);                                                     
C = C' * (C .* repmat(w, 1, N));                                                 
C = 0.5 * (C + C');                                                   
end