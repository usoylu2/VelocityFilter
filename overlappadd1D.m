function X = overlappadd1D(A, B, L)
A = transpose(squeeze(A));  
B = transpose(squeeze(B));
N1=length(A);
M=length(B);
A=[A zeros(1,mod(-N1,L))];
N2=length(A);
B=[B zeros(1,L-1)];
B=fft(B,L+M-1);
S=N2/L;
index=1:L;
index2=1:L+M-1;
X=zeros(1,N2 + M-1);
for stage=1:S
    xm=[A(index) zeros(1,M-1)];		% Selecting sequence to process
    xm=fft(xm,L+M-1);
    xm=xm.*B;
    xm=ifft(xm);
    X(index2) = X(index2) + xm;
    index=stage*L+1:(stage+1)*L;
    index2=stage*L+1:((stage+1)*L + M-1);
end
i=1:N1+M-1;
X= X(i);
X = reshape(X,1,1,N1+M-1);

dims = 3;
ifun = @(m,n) ceil((n-1)/2)+(1:m);
subs(1:dims) = {':'};
subs{3} = ifun(N1,M);
X = X(subs{:});
end
