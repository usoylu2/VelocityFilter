function phi_c = overlapaddefficient(bubble,filter,kx,kz,kt,option)
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
phi_c = rand(M, 'single');

for p=1:kx
    for j=1:kz
        for i=1:kt
            tStart = tic;
            temp = fftmemory2(bubble((1:M(1)/kx) +(p-1)*M(1)/kx,(1:M(2)/kz) +(j-1)*M(2)/kz,(1:M(3)/kt) +(i-1)*M(3)/kt), filter ,'full');
            index1 = (1:M(1)/kx+ N(1)-1) +(p-1)*M(1)/kx;
            index2 = (1:M(2)/kz+ N(2)-1) +(j-1)*M(2)/kz;
            index3 = (1:M(3)/kt + N(3)-1)+M(3)/kt*(i-1);
            c1 = intersect(index1, subs{1});
            c2 = intersect(index2, subs{2});
            c3 = intersect(index3, subs{3});
            c1new = c1 - subs{1}(1)+1;
            c2new = c2 - subs{2}(1)+1;
            c3new = c3 - subs{3}(1)+1;
            c1x = c1 - index1(1)+1;
            c2x = c2 - index2(1)+1;
            c3x = c3 - index3(1)+1;
            phi_c(c1new,c2new,c3new) = phi_c(c1new,c2new,c3new) + temp(c1x,c2x,c3x);
            temp = [];
            disp(['Expected Run Time:',num2str(toc(tStart)*kx*kz*kt),'seconds']);           
        end
    end
end
%End

end