function [profile] = contactfind(croppedImg, x, y, r,n, verbose)
    %if particle(n).g2 < g2thresh
    debug= false;
    padding = 3;
    maskwidth_fraction = 1/10;
    %if uint16(y+r+padding)<size(Gimg, 1)&uint16(x+r+padding)<size(Gimg,2)&uint16(y1-r1-padding)>1&uint16(x1-r1-padding)>1
%     croppedImg = (Gimg(uint16(y-r-padding):uint16(y+r+padding),uint16(x-r-padding):uint16(x+r+padding)));

    [a,b] = size(croppedImg);
    [X, Y] = meshgrid( (1:b)-r-padding, (1:a)-r-padding);
    R = sqrt(X.^2 + Y.^2);
    mask_1 = (r-r*maskwidth_fraction <R&R<r-1); % smooth 1 px around the radius
    maskedImg = uint8(double(croppedImg).*mask_1);
    test = croppedImg(mask_1);
    
    [hist_counts, d]= imhist(test, 256);
    L = otsuthresh(hist_counts);
    BW = imbinarize(maskedImg,L);
    
    [centers,radii] = imfindcircles(BW,[r-8 r+4]);
    
    if isempty(radii)
        %BWn =BW;
        mask_1 = (r-r*maskwidth_fraction<R & R<r-2);
        %radii = r1;
        finrad = r-1;
        flag = false;
        inds_slice = find(BW == 1);
        avgY = mean(Y(inds_slice));
        avgX = mean(X(inds_slice));
        [X, Y] = meshgrid( (1:b)-r-padding+avgX/4, (1:a)-r+avgY/4-padding);

        R = sqrt(X.^2 + Y.^2);
        theta = atan2(Y, X);
        %mask_1 = (r1-r/2<R & R<r1-2);
        BWn = logical(BW.*mask_1);
            
        [centersn,radiin] = imfindcircles(BWn,[int16(r - 8) int16(r+1)]);
        maskedImg = uint8(double(croppedImg).*mask_1);
        xn = -avgX/4 +x;
        yn = -avgY/4 +y;
            
    elseif ~isempty(radii)
    dif = 10;
    i=1;
%     figuredisp = figure;
%     axdisp = subplot(1,1,1, 'Parent', figuredisp);
    while length(radii)>0
        if dif<5 | i>9 %note that i and dif limits are set arbitrarily
            flag = true;
            %

            mask_1 = (r-r*maskwidth_fraction<R & R<r-1);
            r = radii(1);
            finrad = radii(1)-1;
            BWn = logical(BW.*mask_1);
            [centers,radii] = imfindcircles(BWn,[int16(radii(1)- 4) int16(radii(1)+1)]);
            maskedImg = uint8(double(croppedImg).*mask_1);
            

            radii = [];
        else
            
%           
            dif = radii(1)-1-(r-r/2);
            finrad= radii(1)-1;
            slicemask = (R>radii(1)-1 & R<radii(1));
            BWslice = logical(BW.*slicemask);
            inds_slice = find(BW == 1);
            avgY = mean(Y(inds_slice));
            avgX = mean(X(inds_slice));

            

            [X, Y] = meshgrid( (1:b)-r-padding+avgX/4, (1:a)-r+avgY/4-padding);

            R = sqrt(X.^2 + Y.^2);
            theta = atan2(Y, X);
            mask_1 = (r-r*maskwidth_fraction<R & R<radii(1)-1);
            BWn = logical(BW.*mask_1);

            

            [centersn,radiin] = imfindcircles(BWn,[int16(radii(1)- 8) int16(radii(1)+1)]);
            
            
            maskedImg = uint8(double(croppedImg).*mask_1);
            
            if debug
            figurepanel = figure;
            ax = subplot(2, 2, 1, 'Parent', figurepanel);
            imshow(maskedImg, 'Parent', ax)
            title(num2str(n))
            hold on;
            plot(ax,centers(1), centers(2),'x')
            
            ax3 = subplot(2, 2, 3, 'Parent', figurepanel);
            imshow(BW, 'Parent', ax3)
            hold on;
            plot(ax3,avgX+a/2, avgY+b/2,'x')
            
            ax4 = subplot(2, 2, 4, 'Parent', figurepanel);
            imshow(BWn, 'Parent', ax4)
            
            ax2 = subplot(2,2,2, 'Parent', figurepanel);
              imshow(maskedImg, 'Parent', ax2)
               hold on;

            if ~isempty(centersn)
                
                plot(ax2,centersn(1), centersn(2), 'x')
%                 
                
            end   
            end
            i=i+1;
            BW = BWn;
            if r - r/2 > r-4
                mask_1 = (r-r/2<R & R<radii(1)-4);
            end
            centers = centersn;
            radii = radiin;
            
            
            xn = -avgX/4 +x;
            yn= -avgY/4 +y;
        end
%        plot(axdisp, 1:length(cdisp), cdisp)
    end
    
    end
    
    
rn= finrad+1;
newcoords = [xn, yn, rn];
%newcoords = [x,y, r];

[a,b] = size(croppedImg);
[X, Y] = meshgrid( (1:b)-r-padding, (1:a)-r-padding);
R = sqrt(X.^2 + Y.^2);
m2 = (rn/2<R&R<rn-1);

maskedImg = uint8(double(croppedImg).*m2);
values = double(croppedImg(m2));
theta = atan2(Y, X);
angles = theta(m2);
combo = [angles, values];
profile = sortrows(combo);

profile = smoothdata(profile,1,'sgolay', length(profile)/25);

profX = [profile(1:uint16(length(profile)/4),1)+(2*pi), profile(1:uint16(length(profile)/4), 2)];
profile = [profile; profX];
%end