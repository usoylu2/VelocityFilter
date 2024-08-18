function phi = overlapadd(bubble,filter,kx,kz,kt,option)
% option: 'same';truncated version or 'full';full convolution 
% kx,kz,kt: number of seperations for each dimensions

bubble = single(bubble);
filter = single(filter);    
if nargin<6||isempty(option)
    option = 'same';
end

ifun = @(m,n) ceil((n-1)/2)+(1:m);
subs(1:3) = {':'};
dims = 1:3;
for dim=dims
    m = size(bubble,dim);
    n = size(filter,dim);
    subs{dim} = ifun(m,n);
end

%Overlapp - Add 
M = size(bubble);
N = size(filter);
phi = rand(M+N-1, 'single');
for p=1:kx
    for j=1:kz
        for i=1:kt
            tStart = tic;
            temp = fftmemory2(bubble((1:M(1)/kx) +(p-1)*M(1)/kx,(1:M(2)/kz) +(j-1)*M(2)/kz,(1:M(3)/kt) +(i-1)*M(3)/kt), filter ,'full');
            phi((1:M(1)/kx+ N(1)-1) +(p-1)*M(1)/kx,(1:M(2)/kz+ N(2)-1) +(j-1)*M(2)/kz,(1:M(3)/kt + N(3)-1)+M(3)/kt*(i-1)) = phi((1:M(1)/kx+ N(1)-1) +(p-1)*M(1)/kx,(1:M(2)/kz+ N(2)-1) +(j-1)*M(2)/kz,(1:M(3)/kt + N(3)-1)+M(3)/kt*(i-1))+temp;
            temp = [];
            disp(['Expected Run Time:',num2str(toc(tStart)*kx*kz*kt),'seconds']);           
        end
    end
end
%End

if option == 'same'
    phi = real(phi(subs{:}));
end
end
