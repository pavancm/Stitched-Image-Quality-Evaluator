
function [x,y] = extract_patch(image,block_size)
%function to extract patch

ind_y = 0;
x=[];y=[];x_val=[];y_val=[];c_val=[];

e=1;
while(ind_y + block_size <= size(image,1))    
    ind_x = 0;
    while(ind_x + block_size<= size(image,2))
        block_im = image(ind_y+1:ind_y+block_size,ind_x+1:ind_x+block_size,:);
        
        e = e + 1;
        x = [x;ind_y+1 ind_y+block_size];
        y = [y;ind_x+1 ind_x+block_size];
        ind_x = ind_x + block_size;
    end
    block_im = image(ind_y+1:ind_y+block_size,ind_x+1:size(image,2),:);
    
    %border case, if the border patch is less than half the size of block_size it is rejected
    
    if(size(block_im,1)>=.5*block_size && size(block_im,2)>=.5*block_size)
        e = e + 1;
        x = [x;ind_y+1 ind_y+block_size];
        y = [y;ind_x+1 size(image,2)];
        
    end
    ind_y = ind_y + block_size;
end

ind_x = 0;
while(ind_x + block_size<= size(image,2))
    block_im = image(ind_y+1:size(image,1),ind_x+1:ind_x+block_size,:);
    
    %border case, if the border patch is less than half the size of block_size it is rejected
    
    if(size(block_im,1)>=0.5*block_size && size(block_im,2)>=0.5*block_size)
        e = e + 1;
        x = [x;ind_y+1 size(image,1)];
        y = [y;ind_x+1 ind_x+block_size];
    end
    ind_x = ind_x + block_size;
end

block_im = image(ind_y+1:size(image,1),ind_x+1:size(image,2),:);

%border case, if the border patch is less than half the size of
    %block_size it is rejected
    
if(size(block_im,1)>=0.5*block_size && size(block_im,2)>=0.5*block_size)
    e = e + 1;
    x = [x;ind_y+1 size(image,1)];
    y = [y;ind_x+1 size(image,2)];
end

end