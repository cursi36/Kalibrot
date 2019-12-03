function y = Mycombvec(M)
%MyCOMBVEC Create all combinations of vectors.
% input: M = matrix of vectors [v1;v2;....;vn]
% Each vector has dimension m;
%
% output: y = nxm^n Matrix of combinations in 
%  
% Francesco Cursi, 08-12-19


[rows,~] = size(M);

if length(M) == 0
  y = [];
else
  y = M(1,:);
  for i=2:rows
    z = M(i,:);
    y = [copy_blocked(y,size(z,2)); copy_interleaved(z,size(y,2))];
end
end

%=========================================================
function b = copy_blocked(m,n)

[mr,mc] = size(m);
b = zeros(mr,mc*n);
ind = 1:mc;
for i=[0:(n-1)]*mc
  b(:,ind+i) = m;
end
%=========================================================

function b = copy_interleaved(m,n)

[mr,mc] = size(m);
b = zeros(mr*n,mc);
ind = 1:mr;
for i=[0:(n-1)]*mr
  b(ind+i,:) = m;
end
b = reshape(b,mr,n*mc);
