function [ x_movm ] = TriangMovMean(x,np_win)
%TriangMovMean calculates the moving mean of x with triangular weights
% Input Arguments:
%   x: original series
%   np_win: window size (number of points)

assert(size(x,2) == 1,'Error./n x must be a column vector')

np_h = length(x); %length of history

%weights of moving mean
triag_win_wts = triang(np_win);

%create convolution matrixt
cov_mat_trg_wt = zeros(np_h);
%develop convolution matrix as diagonal elements
for k = 1:np_win
    n_diag = np_h-abs(k-floor(np_win/2)-1);
    cov_mat_trg_wt = cov_mat_trg_wt + diag(triag_win_wts(k)*ones(n_diag,1),k-floor(np_win/2)-1);
end
cov_mat_trg_wt = cov_mat_trg_wt./repmat(sum(cov_mat_trg_wt,2),1,np_h); %adjust sum of weights to equal to one
cov_mat_trg_wt = sparse(cov_mat_trg_wt); %convert it to a sparse matrix to save computation 

x_movm = cov_mat_trg_wt*x;

end