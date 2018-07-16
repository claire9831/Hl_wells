% Created by Rui Yan, April 2018
% This function creates a mask for the well image
function I = mask(Image)
I = Image;
se = strel('disk', 7);
Ie = imerode(I, se);
Iobr = imreconstruct(Ie, I);
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), ...
imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
bw = im2bw(Iobrcbr, graythresh(Iobrcbr)+0.001);

L = bwlabel(bw);
[Ny,Nx] = size(bw);
whichlabeliswell = L(round(Ny/2),round(Nx/2));
wellmask = L == whichlabeliswell;
%I(~wellmask)= min(I(:));
I(~wellmask)= 0;
return;
