function grad_im = find_neighbor(quant_im,q)

switch(q)
    case 1
        grad_im{1} = quant_im(1:size(quant_im,1),1:size(quant_im,2)-1); %Immediate Horizontal 
        grad_im{2} = quant_im(1:size(quant_im,1),2:size(quant_im,2));
    case 2
        grad_im{1} = quant_im(1:size(quant_im,1),1:size(quant_im,2)-2);%Alternate Horizontal 
        grad_im{2} = quant_im(1:size(quant_im,1),3:size(quant_im,2));
    case 3
        grad_im{1} = quant_im(1:size(quant_im,1)-1,1:size(quant_im,2));%Immediate Vertical
        grad_im{2} = quant_im(2:size(quant_im,1),1:size(quant_im,2));
    case 4
        grad_im{1} = quant_im(1:size(quant_im,1)-2,1:size(quant_im,2));%Alternate Vertical
        grad_im{2} = quant_im(3:size(quant_im,1),1:size(quant_im,2));
    case 5
        grad_im{1} = quant_im(1:size(quant_im,1)-1,2:size(quant_im,2));%Immediate Diagonal
        grad_im{2} = quant_im(2:size(quant_im,1),1:size(quant_im,2)-1);
    case 6
        grad_im{1} = quant_im(1:size(quant_im,1)-2,3:size(quant_im,2));%Alternate Diagonal
        grad_im{2} = quant_im(3:size(quant_im,1),1:size(quant_im,2)-2);
    case 7
        grad_im{1} = quant_im(1:size(quant_im,1)-1,1:size(quant_im,2)-1);%Immediate cross diagonal
        grad_im{2} = quant_im(2:size(quant_im,1),2:size(quant_im,2));
    case 8
        grad_im{1} = quant_im(1:size(quant_im,1)-2,1:size(quant_im,2)-2);%Alternate cross diagonal
        grad_im{2} = quant_im(3:size(quant_im,1),3:size(quant_im,2));
end

end