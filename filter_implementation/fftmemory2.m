function A = fftmemory2(A, B,option)
% option: 'same';truncated version or 'full';full convolution 
if nargin<3||isempty(option)
    option = 'same';
end
A = single(A);
B = single(B);
dims = 1:3;
ifun = @(m,n) ceil((n-1)/2)+(1:m);
subs(1:3) = {':'};
l = [0,0,0];

for dim=dims
    m = size(A,dim);
    n = size(B,dim);
    subs{dim} = ifun(m,n);
    l(dim) = m+n-1;
end

A = fftn(A,l);
B = fftn(B,l);
A =A.*B;
B = [];
A = ifftn(A);
if option == 'same'
    A = real(A(subs{:}));
end
end



