function A = fftmemory(A, B,option)
% option: 'same';truncated version or 'full';full convolution 
if nargin<3||isempty(option)
    option = 'same';
end
    
A = single(A);
B = single(B);
dims = 1:3;
ifun = @(m,n) ceil((n-1)/2)+(1:m);
subs(1:3) = {':'};

for dim=dims
    m = size(A,dim);
    n = size(B,dim);
    l = m+n-1;
    A = fft(A,l,dim);
    B = fft(B,l,dim);
    subs{dim} = ifun(m,n);
end

A =A.*B;
B = [];
for dim=dims
    A = ifft(A,[],dim);
end 

if option == 'same'
    A = real(A(subs{:}));
end
end



