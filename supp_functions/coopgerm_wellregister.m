function [xwell,ywell] = coopgerm_wellregister(subIm)

archetmp = 3e4*ones(509,509,'uint16'); % template for the ideal well

wellmask = false(509);
inds = sub2ind([509 509],[111 111 399 399],[111 399 111 399]);

wellmask(inds)=1;
wellmask = imdilate(wellmask,strel('disk',110,0));
% wellinner = imdilate(wellmask,strel('disk',106-37,0));

archetmp(~wellmask) = 2100; archetmp(wellmask) = 1200;
[Nytmp,Nxtmp] = size(archetmp);

imcorr = normxcorr2e(archetmp,subIm,'same');

% param(1) = scaling, param(2) = angle, param(3) =

lb = [0.9 % lower bounds on parameters
    -10*pi/180
    0
    0];   
    
ub = [1.1 % upper bounds on parameters
    10*pi/180
    Nxtmp
    Nytmp];

[~,xymaxcorr] = max(imcorr(:));
[ymaxcorr,xmaxcorr]=ind2sub(size(imcorr),xymaxcorr);

params0 = [1; 0; mod(xmaxcorr,294); mod(ymaxcorr,294)]; % starting guesses for the parameters


optimparams = fmincon(@coopgerm_anglescalingx0y0,params0,[],[],[],[],lb,ub,...
    [],[],...
    imcorr,Nxtmp,Nytmp);

scaling = optimparams(1);
angle = optimparams(2);
x0 = optimparams(3);
y0 = optimparams(4);

% calculate lattice vectors

dxh = scaling*294*cos(angle);
dyh = scaling*294*sin(angle);

dxv = -scaling*294*sin(angle);
dyv = scaling*294*cos(angle);

horizshift = -2:6;
vertshift = -2:6;

[horizshift_g,vertshift_g]=meshgrid(horizshift,vertshift);
xs = x0+horizshift_g*dxh+vertshift_g*dxv;
ys = y0+horizshift_g*dyh+vertshift_g*dyv;

xwell = colfilt(xs,[2 2],'sliding',@mean); 
ywell = colfilt(ys,[2 2],'sliding',@mean);

if size(xwell,1)==size(xs,1),
    xwell(end,:) = [];
    ywell(end,:) = [];
end;

if size(xwell,2)==size(xs,2),
    xwell(:,end) = [];
    ywell(:,end) = [];
end;


is_ok = xwell>0 & xwell<=size(imcorr,2) & ywell>0 & ywell<size(imcorr,1);

xwell = xwell(is_ok); ywell = ywell(is_ok);

if ~iscolumn(xwell),
    xwell = xwell(:);
    ywell = ywell(:);
end;

return;

function mfit = coopgerm_anglescalingx0y0(params,imcorr,Nxtmp,Nytmp)

scaling = params(1);
angle = params(2);
x0 = params(3);
y0 = params(4);

% calculate lattice vectors

dxh = scaling*294*cos(angle);
dyh = scaling*294*sin(angle);

dxv = -scaling*294*sin(angle);
dyv = scaling*294*cos(angle);

horizshift = 0:5;
vertshift = 0:5;

[horizshift_g,vertshift_g]=meshgrid(horizshift,vertshift);
xs = x0+horizshift_g*dxh+vertshift_g*dxv;
ys = y0+horizshift_g*dyh+vertshift_g*dyv;

is_ok = xs >= 0.5*Nxtmp & xs <= size(imcorr,2)-0.5*Nxtmp ...
    & ys>=0.5*Nytmp & ys <= size(imcorr,1)-0.5*Nytmp;

xtofit = xs(is_ok);
ytofit = ys(is_ok);

if ~iscolumn(xtofit),
    xtofit = xtofit(:);
    ytofit = ytofit(:);
end;

mfit = -mean(interp2(imcorr,xtofit,ytofit,'*linear'));

return;