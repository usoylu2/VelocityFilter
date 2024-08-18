function [] = kamoimager(data,option)
% options are original or envelope

if nargin<2||isempty(option)
    option = 'original';
end

if option == 'original'
    a = size(data);
    for i=1:a(3)
        imagesc((data(:,:,i)));
        colorbar;
        title(i);
        pause(0.01);
    end
end

if option == 'envelope'
    a = size(data);
    for i=1:a(3)
        imagesc(abs(hilbert((data(:,:,i)))));
        title(i);
        pause(0.01);
    end
end
end

