function contactDetection(directory, fileNames, boundaryType,frameidind,verbose)
% %UNTITLED4 Summary of this function goes here
% %   Detailed explanation goes here
% 
% 


% directory = '/eno/cllee3/DATA/240226/run1/'
% %topDirectory = '/Users/carmenlee/Desktop/20150731reprocesseduniaxial/'
% % %topDirectory = './DATA/test/Step09/'
% fileNames = '75Hz*.jpg' %image format and regex
% frameidind = 15
% %
% files = dir([directory,fileNames])
% boundaryType = "annulus"; %if airtable use "airtable" if annulus use "annulus"
% radiusRange = [40, 57];
% %radiusRange = [45, 78]; %airtable
% 
% verbose = false;
% % comment out this section if running in PeGSDiskPrep.m
%   clear particle;
% % % %
% % % %
% % % % %directory = '/mnt/ncsudrive/c/cllee3/MATLAB/PEGS-master/DATA/';
% directory = '/eno/cllee3/DATA/jekollme/20160711/Steps/step09/'
% % %directory = '/eno/cllee3/DATA/230428/run2/warpedimg/'
% % %directory = './DATA/warpedimg/'
% % directory = './DATA/'
% fileNames = '*0532.jpg'
% boundaryType ="airtable"
% verbose = false
% % %filename = ''
%files = dir([directory, filename]);
% Img = imread([directory,files(1).name]);%particle image
% 
% 
% Rimg = Img(:,:,1);
% Gimg = Img(:,:,2);
% bckgnd = poly2mask([2730,2724,3651,3661],[212,212,6099,6100], length(Gimg), length(Gimg));
% %Gimgfine = histeq(Gimg);
% 
% % Gimgfine = im2double(Gimg);
% % histogram(Gimgfine)
% %  Gimgfine = imadjust(Gimgfine, [13/255, 39/255]);
% %adjusting image contrast
% 
% Gimg=imsubtract(Gimg,Rimg./20);
% Gimg = im2double(Gimg);
% Gimg = Gimg.*~bckgnd;
% %Gimg = Gimg.*(Gimg > 0);
% %Gimg = imadjust(Gimg,[0, 0.6]);
% 
% % figure;
% % imshow(Gimgfine)
% %read in particle position information and put into 'particle' structure
% centersfile = dir([directory, 'centers_tracked.txt'])
% pData = dlmread([directory,centersfile.name],',', 1,0); %Read Position data from centers file
% frame =1;
% data = find(pData(:,1) == frame);
%         
%    N = size(data,1);
%         
%         particle(1:N) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[], 'nonContact', [], 'contactPos', [], 'contactInt', [], 'edge', 0);
%         for n = 1:N %Bookkeeping
%             particle(n).id= pData(data(n),2);
%             particle(n).x = pData(data(n),3); %-xoffset;
%             particle(n).y = pData(data(n),4); %-yoffset;
%             particle(n).r = pData(data(n),5);
%             particle(n).edge = pData(data(n), 6);
%             
%         end
% 
% verbose = true;
% fsigma = 141; %photoelastic stress coefficient
% g2cal = 145; %Calibration Value for the g^2 method, can be computed by joG2cal.m
% dtol = 30; % How far away can the outlines of 2 particles be to still be considered Neighbours
% 
% contactG2Threshold = 1; %sum of g2 in a contact area larger than this determines a valid contact
% CR = 15;
% % % % % % % % % 

%% 
 particleNumber1 = 1566;
 particleNumber2 = 1639;
%verbose = false;
  
% width = 18;
padding =1;

if boundaryType == "airtable"
    minpeakheight = 0.35;
    minpeakprominence = 0.10;
    minpeakprom_main = 0.20;
    pxPerMeter = 0.0077 / 74;
    fsigma = 100; %photoelastic stress coefficient
    g2cal = 100; %Calibration Value for the g^2 method, can be computed by joG2cal.m
    dtol = 15; % How far away can the outlines of 2 particles be to still be considered Neighbours

    contactG2Threshold = 5; %sum of g2 in a contact area larger than this determines a valid contact
    CR = 10;
    imadjust_limits = [0, 1]
    fineimadjust_limits = [0, 100/255]
    rednormal = 2
elseif boundaryType == "annulus"
    directory = [directory, 'warpedimg/']
    
    minpeakheight = 0.10;
    minpeakprominence = 0.02;
    minpeakprom_main = 0.025;
    pxPerMeter = 0.015/939;
    fsigma = 141; %photoelastic stress coefficient
    g2cal = 145; %Calibration Value for the g^2 method, can be computed by joG2cal.m
    dtol = 30; % How far away can the outlines of 2 particles be to still be considered Neighbours

    contactG2Threshold = 1; %sum of g2 in a contact area larger than this determines a valid contact
    CR = 15; 
    imadjust_limits = [0, 0.6]
    fineimadjust_limits = [0/255, 30/255]%[13/255, 39/255]
    rednormal = 8
end

[directory, fileNames]
files = dir([directory, fileNames(1:end-4), '.tif']);

centersfile = dir([directory, 'centers_tracked.txt']);
pData = dlmread([directory,centersfile.name],',', 1,0); %Read Position data from centers file
%frame =1
length(files)
for imgnumb = 1:length(files)
%for imgnumb = 344:length(files)-344
    clear particle
    imgnumb = imgnumb
    frame = str2num(files(imgnumb).name(frameidind:frameidind+3))
    data = find(pData(:,1) == frame);
    if ~isempty(data)
        

        N = size(data,1);

        particle(1:N) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[], 'nonContact', [], 'contactPos', [], 'contactInt', [], 'edge', 0);
        for n = 1:N %Bookkeeping
            particle(n).id= pData(data(n),2);
            particle(n).x = pData(data(n),3); %-xoffset;
            particle(n).y = pData(data(n),4); %-yoffset;
            particle(n).r = round(pData(data(n),5));
            particle(n).edge = pData(data(n), 6);
            particle(n).rm = particle(n).r*pxPerMeter;
            particle(n).fsigma = fsigma;
        end


    Img = imread([directory,files(imgnumb).name]);%particle image
    Rimg = Img(:,:,1);
    Gimg = Img(:,:,2); %force image
   
    
    
    if boundaryType == "airtable"
        %Gimgp=imsubtract(Gimg,Rimg.*0.5);
        Gimgp= im2double(Gimg);
        %Gimg = Gimg-0.5*Rimg;
        Gimgd = Gimgp.*(Gimgp > 0);
        %Gimgd=imsubtract(Gimg,Rimg.*0.5);
        %Gimgd= im2double(Gimgd);
        stretchlim(Gimgp)
        Gimgd = imadjust(Gimgd,stretchlim(Gimgd));
        figure;
        imshow(Gimgp)
        viscircles([pData(data,3), pData(data,4)], round(pData(data,5)))
        title('Gimgp')
        Gimgp = imadjust(Gimgp,[0,0.6]);

%         figure
%         imshow(Gimgd)
    end

 if boundaryType == "annulus"
   bckgnd = poly2mask([2731,2719,3643,3666],[212,212,6099,6100], length(Gimg), length(Gimg));
   Gimg = inpaintCoherent(Gimg,bckgnd,'SmoothingFactor',5,'Radius',15);
 
%    figure;
%    imshow(imadjust(Gimg, [0,0.15]));
%    hold on
%    plot([2731,2719,3643,3666, 2731],[212,212,6099,6100, 212],'b','LineWidth',1)
%    drawnow;
   Gimg=imsubtract(Gimg,Rimg./rednormal);
   sigma  = 50
    G = fspecial('gaussian', 3*sigma+1, sigma);
    yb = imfilter(imcomplement(Rimg), G, 'replicate');
%     figure;
%     imshow(yb)
%     drawnow;
    Gimg = bsxfun(@minus, Gimg,yb*.09);
    
%     figure;
%     imshow(Gimg)
    %sigma = 80; % choosen by visual inspection ephraim
    
   Gimg= im2double(Gimg);
%Gimg = Gimg-0.5*Rimg;
    Gimg = Gimg.*(Gimg > 0);
    
    
    Gimgd = imadjust(Gimg,imadjust_limits);
%     figure;
%     imshow(Gimgd)
%     
%     title('Gimgd')
    %Gimgp = imadjust(Gimg,stretchlim(Gimg));
    Gimgp = Gimg;
    
%     figure;
%     imshow(Gimgp)
%     title('Gimgp')
%     drawnow;
 end
    %Gimg = Gimg.*(Gimg > 0);
    
for n=1:N
    
    %create a circular mask
    % => Find a better way yo do this masking!

    r = particle(n).r*1.0;
    if round(particle(n).y+r)<size(Gimg, 1)&&round(particle(n).x+r)<size(Gimg,2)&&round(particle(n).y-r)>1&&round(particle(n).x-r)>1

        mask = abs(-r:r);
        mask = mask.^2 + mask.^2';
        mask1 = (sqrt(mask) <= r);

        %This crops out a particle
        cropXstart = round(particle(n).x-r);
        cropXstop = round(particle(n).x-r)+ size(mask1,1)-1;
        cropYstart = round(particle(n).y-r);
        cropYstop = round(particle(n).y-r)+ size(mask1,2)-1;
        cimg = im2double(Gimgd(cropYstart:cropYstop, cropXstart:cropXstop));
        %imshow(cimg)
        particleImg = cimg.*mask1;
        %imshow(particleImg)
        particle(n).forceImage=particleImg;

        %create a circular mask with a radius that is one pixel smaller
        %for cropping out the relevant gradient

        mask2 = double(sqrt(mask) <= r-1);

        %Compute G^2 for each particle
        [gx,gy] = gradient(particleImg);
        g2 = (gx.^2 + gy.^2).*mask2;
        particle(n).g2 = sum(sum(g2));
        particle(n).f = particle(n).g2/g2cal;
    else
        error('badimage!!')
        
    end
end



Gimgfine = Gimgp;

% imshow(Gimgfine)

 
xmat = zeros([N,1]);
ymat = zeros([N,1]); % Preallocation
rmat = zeros([N,1]);
Gimgfine = imadjust(Gimgp, fineimadjust_limits);
% figure;
%    imshow(Gimgfine)
for l = 1:N
    
    xmat(l) = particle(l).x;
    ymat(l) = particle(l).y; %Pulls data from particle structure
    rmat(l) = particle(l).r;
    

end

%histogram(rmat)
rmats = rmat; %Saves our radius matrix for later

dmat = pdist2([xmat,ymat],[xmat,ymat]); %Creates a distance matrix for particle center locations
rmat = rmat + rmat'; %Makes a combination of radii for each particle

friendmat = dmat < (rmat + dtol) & dmat~=0; %Logical "friend" matrix

friendmat = triu(friendmat); %Only examine the upper triangle portion (no repeats)
[f1, f2] = find(friendmat == 1); %Creates an index of particles that are considered touching
%% 




mask = abs(-CR:CR);
mask = mask.^2 + mask.^2';
mask = double(sqrt(mask) <= CR-1);
%% 

for l = 1:length(f1)
    %if f1(l) == particleNumber1 & f2(l) == particleNumber2
    x1 = particle(f1(l)).x;
    y1 = particle(f1(l)).y;
    r1 = particle(f1(l)).r;
    x2 = particle(f2(l)).x;
    y2 = particle(f2(l)).y;
    r2 = particle(f2(l)).r;
    %if uint16(y1+r1+padding)<size(Gimg, 1)&uint16(x1+r1+padding)<size(Gimg,2)&uint16(y1-r1-padding)>1&uint16(x1-r1-padding)>1&uint16(y2+r2+padding)<size(Gimg, 1)&uint16(x2+r2+padding)<size(Gimg,2)&uint16(y2-r2-padding)>1&uint16(x2-r2-padding)>1
    
    
    
    contactXp1 = round(x1 + (r1-1 - CR) * cos(atan2(y2-y1,x2-x1)));
    contactYp1 = round(y1 + (r1 -1- CR) * sin(atan2(y2-y1,x2-x1)));
    
    contactXp2 = round(x1 + (r1-1 + CR + dmat(f1(l),f2(l)) - rmat(f1(l),f2(l))) * cos(atan2(y2-y1,x2-x1))); %%note to carm, we might want to filter by general g2 first
    contactYp2 = round(y1 + (r1-1 + CR + dmat(f1(l),f2(l)) - rmat(f1(l),f2(l))) * sin(atan2(y2-y1,x2-x1)));
    contactImg = im2double(imcrop(Gimgd,[contactXp1-CR contactYp1-CR CR*2 CR*2]));
    
    contactImg = contactImg.*mask;
    
    [gx,gy] = gradient(contactImg);
    g2 = (gx.^2 + gy.^2);
    contactG2p1 = sum(sum(g2));
    contactIp1 = sum(sum(contactImg));
    
    contactImg = im2double(imcrop(Gimgd,[contactXp2-CR contactYp2-CR CR*2 CR*2]));
%     if l==2107
%         n
%         imshow(contactImg)
%         figure
%         imshow(mask)
%         drawnow;
%         n
%     end
%     l
    contactImg = contactImg.*mask;
    contactG2p2 = gradientcalculator(contactImg);
    
    contactIp2 = sum(sum(contactImg));



    
    %[m1,n1,l1] = size(contactG2p1)
    %[m2,n2,l2] = size(contactG2p2)
    %if we declare our contact valid
%      viscircles(ax1, [contactXp1; contactYp1]', CR,'EdgeColor','w');           
%      viscircles(ax1,[contactXp2; contactYp2]', CR,'EdgeColor','w');
%      text(ax1,contactXp1, contactYp1,num2str(contactG2p1),'Color','r');
%      text(ax1, contactXp2, contactYp2,num2str(contactG2p2),'Color','r');
    if(contactG2p1 > contactG2Threshold && contactG2p2 > contactG2Threshold)
        
        %Plot contact area
       
        %this is a valid contact, remember it
        particle(f1(l)).z= particle(f1(l)).z +1; %increase coordination number
        particle(f1(l)).contactG2s(particle(f1(l)).z)=contactG2p1; %remember the g2 value of the current contact area
        particle(f1(l)).contactIs(particle(f1(l)).z)=contactIp1; %changes to color
        particle(f1(l)).color(particle(f1(l)).z)='r'; %changes to color
        particle(f1(l)).neighbours(particle(f1(l)).z) = particle(f2(l)).id; %particle m is now noted as a neigbour in the particle l datastructure
        particle(f1(l)).betas(particle(f1(l)).z) = atan2(y2-y1,x2-x1); %the contact angle to particle m is now noted in the particle l datastructure
        particle(f2(l)).z= particle(f2(l)).z +1; %increase coordination number
        particle(f2(l)).contactG2s(particle(f2(l)).z)=contactG2p2; %remember the g2 value of the current contact area
        particle(f2(l)).contactIs(particle(f2(l)).z)=contactIp2;
        particle(f2(l)).color(particle(f2(l)).z)='r'; %changes to color
        particle(f2(l)).neighbours(particle(f2(l)).z) = particle(f1(l)).id; %particle m is now noted as a neigbour in the particle l datastructure
        particle(f2(l)).betas(particle(f2(l)).z) = atan2(y1-y2,x1-x2);
  
    else %we try the more refined method of contact detection
        
        contactImg = im2double(imcrop(Gimgfine,[contactXp1-CR contactYp1-CR CR*2 CR*2]));
    
        contactImg = contactImg.*mask;
        contactG2p1=gradientcalculator(contactImg);
        
        contactIp1 = sum(sum(contactImg));

        contactImg = im2double(imcrop(Gimgfine,[contactXp2-CR contactYp2-CR CR*2 CR*2]));
        contactImg = contactImg.*mask;
        contactG2p2 = gradientcalculator(contactImg);
        
        contactIp2 = sum(sum(contactImg));
%         text(ax1,contactXp1, contactYp1+10,num2str(contactG2p1),'Color','b');
%         text(ax1, contactXp2, contactYp2+10,num2str(contactG2p2),'Color','b');
%         if contactG2p1> contactG2Threshold && contactG2p2 > contactG2Threshold
%             particle(f1(l)).z= particle(f1(l)).z +1; %increase coordination number
%             particle(f1(l)).contactG2s(particle(f1(l)).z)=contactG2p1; %remember the g2 value of the current contact area
%             particle(f1(l)).contactIs(particle(f1(l)).z)=contactIp1; %changes to color
%             particle(f1(l)).color(particle(f1(l)).z)='w'; %changes to color
%             particle(f1(l)).neighbours(particle(f1(l)).z) = particle(f2(l)).id; %particle m is now noted as a neigbour in the particle l datastructure
%             particle(f1(l)).betas(particle(f1(l)).z) = atan2(y2-y1,x2-x1); %the contact angle to particle m is now noted in the particle l datastructure
%             particle(f2(l)).z= particle(f2(l)).z +1; %increase coordination number
%             particle(f2(l)).contactG2s(particle(f2(l)).z)=contactG2p2; %remember the g2 value of the current contact area
%             particle(f2(l)).contactIs(particle(f2(l)).z)=contactIp2;
%             particle(f2(l)).color(particle(f2(l)).z)='w'; %changes to color
%             particle(f2(l)).neighbours(particle(f2(l)).z) = particle(f1(l)).id; %particle m is now noted as a neigbour in the particle l datastructure
%             particle(f2(l)).betas(particle(f2(l)).z) = atan2(y1-y2,x1-x2);
%         else

        croppedImg = (Gimgfine(uint16(y1-r1-padding):uint16(y1+r1+padding),uint16(x1-r1-padding):uint16(x1+r1+padding)));
        %croppedImg = im2uint8(croppedImg);
        %croppedImg = imadjust(croppedImg, [0.2, 175/255]);
        
        %croppedImg = imadjust(croppedImg, [0, 0.15]);
        [profile] = contactfind(croppedImg, x1, y1, r1-1, f1(l), verbose);
        
        if any(profile(:,1)>minpeakheight)
        [pkints, locints] = findpeaks(profile(:,2),profile(:,1), "MinPeakHeight", minpeakheight,'MinPeakProminence', minpeakprominence, 'MaxPeakWidth', pi);
        else
            pkints = [];
            locints = [];
        end
        [pks,locs] =findpeaks(profile(:,2), profile(:,1),'MinPeakProminence',minpeakprom_main,'MaxPeakWidth', pi);
        pks = [pks;pkints];
        locs = [locs;locints];
        inds = find(locs>(pi));
        locs(inds) = locs(inds)-2*pi;
        
        croppedImg2 = (Gimgfine(uint32(y2-r2-padding):uint32(y2+r2+padding),uint32(x2-r2-padding):uint32(x2+r2+padding)));
        [profile2] = contactfind(croppedImg2, x2, y2, r2-1, f2(l), verbose);


        [pks2,locs2] =findpeaks(profile2(:,2), profile2(:,1),'MinPeakProminence',minpeakprom_main,'MaxPeakWidth', pi);
        if any(profile2(:,1)>minpeakheight)
            [pkints2, locints2] = findpeaks(profile2(:,2),profile2(:,1), "MinPeakHeight", minpeakheight,'MinPeakProminence', 0.02, 'MaxPeakWidth', pi);
        else
            pkints2 = [];
            locints2 = [];
        end
        pks2 = [pks2;pkints2];
        locs2 = [locs2;locints2];
        inds2 = find(locs2>(pi));
        locs2(inds2) = locs2(inds2)-2*pi;
        

        
    %%Merging peaks that are probably too close together
    
    
        newlocs = [];
        newpks = [];
        len =1;
        

        while len <= length(locs)
            matching = find((locs > locs(len)-pi/6) & (locs< locs(len)+pi/6));
            if length(matching) >= 2
                if matching(1)==len
                newlocs = [newlocs,(locs(matching(1))+locs(matching(2)))/2];
                newpks = [newpks,pks(matching(1))];
                end
             elseif length(matching) ==1
               
               newlocs = [newlocs,locs(matching(1))];
               newpks = [newpks,pks(matching(1))];
             end
             len = len+1;
        end
        newlocs = newlocs';

        newlocs2 = [];
        newpks2 = [];
        len =1;
        while len <= length(locs2)
            matching2 = find((locs2 > locs2(len)-pi/6) & (locs2< locs2(len)+pi/6));
            if length(matching2) >= 2
             if matching2(1)==len
        
                newlocs2 = [newlocs2,(locs2(matching2(1))+locs2(matching2(2)))/2];
                  
                newpks2 = [newpks2,pks2(matching2(1))];
             elseif pks2(matching2(1)-pks2(matching2(2) < 0.001))
                 newlocs2 = [newlocs2,locs2(matching2(1))];
             end
             elseif length(matching2) ==1
%                if locs2(matching2(1)) > 0
%                newlocs2(matching2(1)) = locs2(matching2(1));
%                else
%                    newlocs2(matching2(1)) = locs2(matching2(1));
%                end
                newlocs2 = [newlocs2,locs2(matching2(1))];
               newpks2 = [newpks2,pks2(matching2(1))];
            end
            len = len+1;
        end
        
        
        newlocs2 = (newlocs2)';
        
        invnewlocs2 = zeros(size(newlocs2));
        for angle=1:length(newlocs2)
            if newlocs2(angle)>0
                invnewlocs2(angle) = newlocs2(angle)-pi;
            elseif newlocs2(angle)<0
                invnewlocs2(angle) = newlocs2(angle)+pi;
            else
                invnewlocs2(angle) = newlocs2(angle);
            end
        end
        nominalAngle = atan2(y2-y1, x2-x1);
        if ismember(1,ismembertol([-pi, pi], nominalAngle, 0.2))
            
            if nominalAngle<0
                nominalAngle = [nominalAngle, nominalAngle+2*pi];
            else
                nominalAngle = [nominalAngle, nominalAngle-2*pi];
            end
        end
        angle = ismembertol(newlocs, nominalAngle, 0.25);
        angle2 = ismembertol(invnewlocs2, nominalAngle, 0.25);
        

        nominalAngle = atan2(y2-y1, x2-x1);
         if ismember(1, angle) && ismember(1, angle2)

            if x1 > x2
                xcoords = [x2,x1];
            else
                xcoords = [x1,x2];
            end
            if y1 > y2
                ycoords = [y2,y1];
            else
                ycoords = [y1,y2];
            end
            
            if nominalAngle>0 && nominalAngle <pi/2
                rotAngle = -(pi/2-nominalAngle);
            elseif nominalAngle>pi/2 && nominalAngle < pi
                rotAngle = -(pi/2-nominalAngle);
            elseif nominalAngle<0 && nominalAngle >-pi/2
                rotAngle = pi/2+nominalAngle;
            elseif nominalAngle<-pi/2 && nominalAngle > -pi
                rotAngle = nominalAngle + pi/2;
            end
            


            rotatedImage = imrotate(Gimgfine(uint16(ycoords(1)-r1):uint16(ycoords(2)+r1),uint16(xcoords(1)-r1):uint16(xcoords(2)+r1)), rotAngle*180/pi);
            
            
            
            rotimgsize = size(rotatedImage);
            subImage = rotatedImage(uint16(rotimgsize(1)*2/5):uint16(rotimgsize(1)/5)*3, uint16(rotimgsize(2)*2/5):uint16(rotimgsize(2)*3/5));
            intProfile = mean(subImage, 2);
            
            
            
%           
%             if any(intProfile(:,1)>0.15)
%                 [peaks, loc] = findpeaks(profile2(:,2),profile2(:,1), "MinPeakHeight",0.15,'MinPeakProminence', 0.05, 'MaxPeakWidth', pi);
%             else
%                 peaks = [];
%                 loc = [];
%             end
            [peaks, loc] = findpeaks(intProfile, 'MinPeakHeight',0.15, 'MinPeakProminence', 0.05);
            
            [gmag, gdir] = imgradient((subImage));
            dirProfile = mean(gmag, 2);
            

            if verbose
            if f1(l) ==particleNumber1 && f2(l) ==particleNumber2
            figurerot = figure;
            axes1 = subplot(1,6, 1, 'Parent', figurerot);
            axes2 =subplot(1,6, 2, 'Parent', figurerot);
            axes3 =subplot(1,6, 3, 'Parent', figurerot);
            axes4 =subplot(1,6, 4, 'Parent', figurerot);
            
            imshow(Gimgfine(ycoords(1)-r1:ycoords(2)+r1,xcoords(1)-r1:xcoords(2)+r1), 'Parent', axes1)
            hold on;
            title(axes1,num2str(nominalAngle))
            imshow(rotatedImage, 'Parent', axes2);
            title(num2str(rotAngle), 'Parent', axes2)
            imshow(subImage, 'Parent', axes3);
            plot(axes4,intProfile); hold on;
            plot(axes4,loc, peaks, 'o');

            
            axes5 =subplot(1,6, 5, 'Parent', figurerot);

            imshow(uint8(gmag), 'Parent', axes5);
            
            axes6 =subplot(1,6, 6, 'Parent', figurerot);
            
            plot(axes6, dirProfile)
            end
            
            end
            if length(peaks) >=2 
                if any(intProfile(loc(1):loc(2))<0.001)
                    valid = true;
                    'a';
                else 
                    valid =true;
                    'b';
                end
            elseif length(peaks) ==1
                if any(dirProfile>0.1)
                    valid = true;
                    'c';
                else 
                    valid = false ;
                    'd';
                end
            else
                valid = false;
                'e';
            end
%             if particleNumber2 == f2(l) && particleNumber1 == f1(l)
%                 error()
%             end
            valid =true;
            if (contactG2p2 > 0.33 && contactG2p1 >0.5) || (contactG2p2 > 0.5 && contactG2p1 >0.33) || valid == true
            particle(f1(l)).z= particle(f1(l)).z +1; %increase coordination number
            particle(f1(l)).contactG2s(particle(f1(l)).z)=contactG2p1; %remember the g2 value of the current contact area
            particle(f1(l)).contactIs(particle(f1(l)).z)=contactIp1;
            particle(f1(l)).color(particle(f1(l)).z)='y';
            particle(f1(l)).neighbours(particle(f1(l)).z) = particle(f2(l)).id; %particle m is now noted as a neigbour in the particle l datastructure
            particle(f1(l)).betas(particle(f1(l)).z) = nominalAngle; %the contact angle to particle m is now noted in the particle l datastructure
            particle(f2(l)).z= particle(f2(l)).z +1; %increase coordination number
            particle(f2(l)).contactG2s(particle(f2(l)).z)=contactG2p2; %remember the g2 value of the current contact area
            particle(f2(l)).contactIs(particle(f2(l)).z)=contactIp2;
            particle(f2(l)).color(particle(f2(l)).z)='y';
            particle(f2(l)).neighbours(particle(f2(l)).z) = particle(f1(l)).id; %particle m is now noted as a neigbour in the particle l datastructure
            particle(f2(l)).betas(particle(f2(l)).z) = nominalAngle - pi;
            else
                s = length(particle(f1(l)).nonContact);
                w = length(particle(f2(l)).nonContact);
                particle(f1(l)).nonContact(s+1) = f2(l);
                particle(f2(l)).nonContact(w+1) = f1(l);
            end
         else
                s = length(particle(f1(l)).nonContact);
                w = length(particle(f2(l)).nonContact);
                particle(f1(l)).nonContact(s+1) = f2(l);
                particle(f2(l)).nonContact(w+1) = f1(l);
        end
        
        if f1(l) == particleNumber1 && f2(l) == particleNumber2
        figure1 = figure;
        subplot(2,2,1, 'Parent', figure1)

       
        
        imshow(croppedImg2);
        title([num2str(f2(l)),',', num2str(particle(f2(l)).g2)]);
        
        hold on;
        contactLoc2 = [r2+padding+r2.*cos(newlocs2), r2+padding+r2.*sin(newlocs2)];
        
        plot(r2+padding+r2.*cos(nominalAngle-pi), r2+padding+r2.*sin(nominalAngle-pi),'yo')
        if ~isempty(contactLoc2)
        plot(contactLoc2(:,1), contactLoc2(:,2), 'ro');
        plot(r2+padding, r2+padding, 'gx')
        end
        subplot(2,2,3, 'Parent', figure1);
        plot(profile2(:,1), profile2(:,2));
        title(num2str(f2(l)))
        hold on;
        plot(newlocs2, newpks2, 'o');
        
        subplot(2, 2, 2, 'Parent', figure1);
        
        
        %mask_1 = (r-r*maskwidth_fraction<R & R<r-1);
        imshow(croppedImg);
        %maskedImg = uint8(double(croppedImg).*mask_1);
        title([num2str(f1(l)),',', num2str(particle(f1(l)).g2)]);
        hold on;
        contactLoc = [r1+padding+r1.*cos(newlocs), r1+padding+r1.*sin(newlocs)];
        hold on;
        plot(r1+padding+r1.*cos(nominalAngle), r1+padding+r1.*sin(nominalAngle),'yo', 'MarkerSize', 10)
        if ~isempty(contactLoc)
        plot(contactLoc(:,1), contactLoc(:,2), 'ro');
        plot(r1+padding, r1+padding, 'gx')
        plot(r1+r1*cos(0), r1+r1*sin(0), 'bo')
        end
        subplot(2,2,4, 'Parent', figure1)
        plot(profile(:,1), profile(:,2));
        hold on;
        plot(newlocs, newpks, 'o');
        end
        
%         subplot(2,3,6, 'Parent', figure1)
%         plot(c);
%         hold on;
%         plot(loc, peaks, 'o');
        
        
        
        if isempty(particle(f1(l)).contactPos)
        particle(f1(l)).contactPos = newlocs;
        particle(f1(l)).contactInt = newpks;
        end
    if isempty(particle(f2(l)).contactPos)
        particle(f2(l)).contactPos = newlocs2;
        particle(f2(l)).contactInt = newpks2;
    end
        %end
    end
    

    %end
    %end
    
end
%% 

%Check if any of the walls is a neighbour as well

if boundaryType == "airtable"
circs = [[particle.y]', [particle.x]', [particle.r]']; %Makes a circs matrix from old matrices

rightwall = max(circs(:,2) + circs(:,3));
leftwall = min(circs(:,2) - circs(:,3)); %Finds our theorhetical wall locations
topwall = min(circs(:,1) - circs(:,3));
bottomwall = max(circs(:,1) + circs(:,3));

rwi = find(circs(:,2) + circs(:,3) + dtol*1.5 >= rightwall);
lwi = find(circs(:,2) - circs(:,3) - dtol*1.5 <= leftwall); %Indexes based on particles that would be considered to be touching the wall
bwi = find(circs(:,1) + circs(:,3) + dtol*1.5 >= bottomwall);
twi = find(circs(:,1) - circs(:,3) - dtol*1.5 <= topwall);

for l = 1:length(lwi) %Runs through each index to check for contacts via gradients
    x = circs(lwi(l),2);
    y = circs(lwi(l),1);
    r = circs(lwi(l),3);
    
    contactX = x-(r-CR);
    contactY = y;
    
    contactImg = im2double(imcrop(Gimgd,[contactX-CR contactY-CR CR*2 CR*2]));
    contactImg = contactImg.*mask;
    contactG2 = gradientcalculator(contactImg);
    
    
    if(contactG2 > contactG2Threshold)
        cI = sum(sum(contactImg));

        particle(lwi(l)).z= particle(lwi(l)).z +1; %increase coordination number
        particle(lwi(l)).contactG2s(particle(lwi(l)).z)=contactG2;
        particle(lwi(l)).contactIs(particle(lwi(l)).z)=cI;
        particle(lwi(l)).neighbours(particle(lwi(l)).z) = -1; %the wall is now noted as a neigbour in the particle l datastructure
        particle(lwi(l)).betas(particle(lwi(l)).z) = pi; %the contact angle to the wall is now noted in the particle l datastructure
        particle(lwi(l)).color(particle(lwi(l)).z)='g';
    else
        croppedImg = (Gimgfine(uint16(y-r-padding):uint16(y+r+padding),uint16(x-r-padding):uint16(x+r+padding)));
        [profile] = contactfind(croppedImg, x, y, r-1, lwi(l), verbose);
        
        if any(profile(:,1)>minpeakheight)
        [pkints, locints] = findpeaks(profile(:,2),profile(:,1), "MinPeakHeight", minpeakheight,'MinPeakProminence', minpeakprominence, 'MaxPeakWidth', pi);
        else
            pkints = [];
            locints = [];
        end
        [pks,locs] =findpeaks(profile(:,2), profile(:,1),'MinPeakProminence',minpeakprom_main,'MaxPeakWidth', pi);
        pks = [pks;pkints];
        locs = [locs;locints];
        inds = find(locs>(pi));
        locs(inds) = locs(inds)-2*pi;

        newlocs = [];
        newpks = [];
        len =1;
        

        while len <= length(locs)
            matching = find((locs > locs(len)-pi/6) & (locs< locs(len)+pi/6));
            if length(matching) >= 2
                if matching(1)==len
                newlocs = [newlocs,(locs(matching(1))+locs(matching(2)))/2];
                newpks = [newpks,pks(matching(1))];
                end
             elseif length(matching) ==1
               
               newlocs = [newlocs,locs(matching(1))];
               newpks = [newpks,pks(matching(1))];
             end
             len = len+1;
        end
        newlocs = newlocs';


        nominalAngle = pi;
        angle = ismembertol(newlocs, nominalAngle, 0.25);
        if ismember(1, angle)
            particle(lwi(l)).z= particle(lwi(l)).z +1; %increase coordination number
        particle(lwi(l)).contactG2s(particle(lwi(l)).z)=contactG2;
        %particle(lwi(l)).contactIs(particle(lwi(l)).z)=cI;
        particle(lwi(l)).neighbours(particle(lwi(l)).z) = -1; %the wall is now noted as a neigbour in the particle l datastructure
        particle(lwi(l)).betas(particle(lwi(l)).z) = pi; %the contact angle to the wall is now noted in the particle l datastructure
        particle(lwi(l)).color(particle(lwi(l)).z)='b';
        end
    end
end

for l = 1:length(rwi)
    x = circs(rwi(l),2);
    y = circs(rwi(l),1);
    r = circs(rwi(l),3);
    
    contactX = x+(r-CR);
    contactY = y;
    
    contactImg = im2double(imcrop(Gimgd,[contactX-CR contactY-CR CR*2 CR*2]));
    contactImg = contactImg.*mask;
    contactG2 = gradientcalculator(contactImg);
    
    
    if(contactG2 > contactG2Threshold)
        cI = sum(sum(contactImg));
        %if(cI > contactIThreshold)
        %this is a valid contact, remember it
%         if(verbose)
%             text(contactX,contactY,num2str(contactG2),'Color','w');
%             viscircles([contactX; contactY]', CR,'EdgeColor','w');
%         end
        particle(rwi(l)).z= particle(rwi(l)).z +1; %increase coordination number
        particle(rwi(l)).contactG2s(particle(rwi(l)).z)=contactG2;
        %particle(rwi(l)).contactIs(particle(rwi(l)).z)=cI;
        particle(rwi(l)).neighbours(particle(rwi(l)).z) = -2; %the wall is now noted as a neigbour in the particle l datastructure
        particle(rwi(l)).betas(particle(rwi(l)).z) = 0; %the contact angle to the wall is now noted in the particle l datastructure
        particle(rwi(l)).color(particle(rwi(l)).z)='g';
    else
        croppedImg = (Gimgfine(uint16(y-r-padding):uint16(y+r+padding),uint16(x-r-padding):uint16(x+r+padding)));
        [profile] = contactfind(croppedImg, x, y, r-1, rwi(l), verbose);
        
        if any(profile(:,1)>minpeakheight)
        [pkints, locints] = findpeaks(profile(:,2),profile(:,1), "MinPeakHeight", minpeakheight,'MinPeakProminence', minpeakprominence, 'MaxPeakWidth', pi);
        else
            pkints = [];
            locints = [];
        end
        [pks,locs] =findpeaks(profile(:,2), profile(:,1),'MinPeakProminence',minpeakprom_main,'MaxPeakWidth', pi);
        pks = [pks;pkints];
        locs = [locs;locints];
        inds = find(locs>(pi));
        locs(inds) = locs(inds)-2*pi;

        newlocs = [];
        newpks = [];
        len =1;
        

        while len <= length(locs)
            matching = find((locs > locs(len)-pi/6) & (locs< locs(len)+pi/6));
            if length(matching) >= 2
                if matching(1)==len
                newlocs = [newlocs,(locs(matching(1))+locs(matching(2)))/2];
                newpks = [newpks,pks(matching(1))];
                end
             elseif length(matching) ==1
               
               newlocs = [newlocs,locs(matching(1))];
               newpks = [newpks,pks(matching(1))];
             end
             len = len+1;
        end
        newlocs = newlocs';


        nominalAngle = 0;
        angle = ismembertol(newlocs, nominalAngle, 0.25);
        if ismember(1, angle)
            particle(rwi(l)).z= particle(rwi(l)).z +1; %increase coordination number
        particle(rwi(l)).contactG2s(particle(rwi(l)).z)=contactG2;
        %particle(rwi(l)).contactIs(particle(rwi(l)).z)=cI;
        particle(rwi(l)).neighbours(particle(rwi(l)).z) = -2; %the wall is now noted as a neigbour in the particle l datastructure
        particle(rwi(l)).betas(particle(rwi(l)).z) = 0; %the contact angle to the wall is now noted in the particle l datastructure
        particle(rwi(l)).color(particle(rwi(l)).z)='b';
        end
    end
end

for l = 1:length(twi)
    x = circs(twi(l),2);
    y = circs(twi(l),1);
    r = circs(twi(l),3);
    
    contactX = x;
    contactY = y-(r-CR);
    
    contactImg = im2double(imcrop(Gimgd,[contactX-CR contactY-CR CR*2 CR*2]));
    contactImg = contactImg.*mask;
    
    contactG2 = gradientcalculator(contactImg);
    
    if(contactG2 > contactG2Threshold)
        cI = sum(sum(contactImg));
        %if(cI > contactIThreshold)
        %this is a valid contact, remember it
%         if(verbose)
%             text(contactX,contactY,num2str(contactG2),'Color','w');
%             viscircles([contactX; contactY]', CR,'EdgeColor','w');
%         end
        particle(twi(l)).z= particle(twi(l)).z +1; %increase coordination number
        particle(twi(l)).contactG2s(particle(twi(l)).z)=contactG2;
        %particle(twi(l)).contactIs(particle(twi(l)).z)=cI;
        particle(twi(l)).neighbours(particle(twi(l)).z) = -3; %the wall is now noted as a neigbour in the particle l datastructure
        particle(twi(l)).betas(particle(twi(l)).z) = -pi/2; %the contact angle to the wall is now noted in the particle l datastructure
        particle(twi(l)).color(particle(twi(l)).z)='g';
   else
        croppedImg = (Gimgfine(uint16(y-r-padding):uint16(y+r+padding),uint16(x-r-padding):uint16(x+r+padding)));
        [profile] = contactfind(croppedImg, x, y, r-1, twi(l), verbose);
        
        if any(profile(:,1)>minpeakheight)
        [pkints, locints] = findpeaks(profile(:,2),profile(:,1), "MinPeakHeight", minpeakheight,'MinPeakProminence', minpeakprominence, 'MaxPeakWidth', pi);
        else
            pkints = [];
            locints = [];
        end
        [pks,locs] =findpeaks(profile(:,2), profile(:,1),'MinPeakProminence',minpeakprom_main,'MaxPeakWidth', pi);
        pks = [pks;pkints];
        locs = [locs;locints];
        inds = find(locs>(pi));
        locs(inds) = locs(inds)-2*pi;

        newlocs = [];
        newpks = [];
        len =1;
        

        while len <= length(locs)
            matching = find((locs > locs(len)-pi/6) & (locs< locs(len)+pi/6));
            if length(matching) >= 2
                if matching(1)==len
                newlocs = [newlocs,(locs(matching(1))+locs(matching(2)))/2];
                newpks = [newpks,pks(matching(1))];
                end
             elseif length(matching) ==1
               
               newlocs = [newlocs,locs(matching(1))];
               newpks = [newpks,pks(matching(1))];
             end
             len = len+1;
        end
        newlocs = newlocs';


        nominalAngle = -pi/2;
        angle = ismembertol(newlocs, nominalAngle, 0.25);
        if ismember(1, angle)
            particle(twi(l)).z= particle(twi(l)).z +1; %increase coordination number
        particle(twi(l)).contactG2s(particle(twi(l)).z)=contactG2;
        %particle(twi(l)).contactIs(particle(twi(l)).z)=cI;
        particle(twi(l)).neighbours(particle(twi(l)).z) = -3; %the wall is now noted as a neigbour in the particle l datastructure
        particle(twi(l)).betas(particle(twi(l)).z) = -pi/2; %the contact angle to the wall is now noted in the particle l datastructure
        particle(twi(l)).color(particle(twi(l)).z)='b';
        end 
    end
end

for l = 1:length(bwi)
    x = circs(bwi(l),2);
    y = circs(bwi(l),1);
    r = circs(bwi(l),3);
    
    contactX = x;
    contactY = y+(r-CR);
    
    contactImg = im2double(imcrop(Gimgd,[contactX-CR contactY-CR CR*2 CR*2]));
    contactImg = contactImg.*mask;
    
    [gx,gy] = gradient(contactImg);
    g2 = (gx.^2 + gy.^2);
    contactG2 = sum(sum(g2));
    
    if(contactG2 > contactG2Threshold)
        cI = sum(sum(contactImg));
        %if(cI > contactIThreshold)
        %this is a valid contact, remember it
%         if(verbose)
%             text(contactX,contactY,num2str(contactG2),'Color','w');
%             viscircles([contactX; contactY]', CR,'EdgeColor','w');
%         end
        particle(bwi(l)).z= particle(bwi(l)).z +1; %increase coordination number
        particle(bwi(l)).contactG2s(particle(bwi(l)).z)=contactG2;
        %particle(bwi(l)).contactIs(particle(bwi(l)).z)=cI;
        particle(bwi(l)).neighbours(particle(bwi(l)).z) = -4; %the wall is now noted as a neigbour in the particle l datastructure
        particle(bwi(l)).betas(particle(bwi(l)).z) = pi/2; %the contact angle to the wall is now noted in the particle l datastructure
        particle(bwi(l)).color(particle(bwi(l)).z)='g';
        else
        croppedImg = (Gimgfine(uint16(y-r-padding):uint16(y+r+padding),uint16(x-r-padding):uint16(x+r+padding)));
        [profile] = contactfind(croppedImg, x, y, r-1, bwi(l), verbose);
        
        if any(profile(:,1)>minpeakheight)
        [pkints, locints] = findpeaks(profile(:,2),profile(:,1), "MinPeakHeight", minpeakheight,'MinPeakProminence', minpeakprominence, 'MaxPeakWidth', pi);
        else
            pkints = [];
            locints = [];
        end
        [pks,locs] =findpeaks(profile(:,2), profile(:,1),'MinPeakProminence',minpeakprom_main,'MaxPeakWidth', pi);
        pks = [pks;pkints];
        locs = [locs;locints];
        inds = find(locs>(pi));
        locs(inds) = locs(inds)-2*pi;

        newlocs = [];
        newpks = [];
        len =1;
        

        while len <= length(locs)
            matching = find((locs > locs(len)-pi/6) & (locs< locs(len)+pi/6));
            if length(matching) >= 2
                if matching(1)==len
                newlocs = [newlocs,(locs(matching(1))+locs(matching(2)))/2];
                newpks = [newpks,pks(matching(1))];
                end
             elseif length(matching) ==1
               
               newlocs = [newlocs,locs(matching(1))];
               newpks = [newpks,pks(matching(1))];
             end
             len = len+1;
        end
        newlocs = newlocs';


        nominalAngle = pi/2;
        angle = ismembertol(newlocs, nominalAngle, 0.25);
        if ismember(1, angle)
           particle(bwi(l)).z= particle(bwi(l)).z +1; %increase coordination number
        particle(bwi(l)).contactG2s(particle(bwi(l)).z)=contactG2;
        %particle(bwi(l)).contactIs(particle(bwi(l)).z)=cI;
        particle(bwi(l)).neighbours(particle(bwi(l)).z) = -4; %the wall is now noted as a neigbour in the particle l datastructure
        particle(bwi(l)).betas(particle(bwi(l)).z) = pi/2; %the contact angle to the wall is now noted in the particle l datastructure
        particle(bwi(l)).color(particle(bwi(l)).z)='b';
        end 
    end
end

end

%%
%%now let's go back through everyone and look for single contacts

contactnumbers = [particle.z];
% figure;
% histogram(contactnumbers)
contactIntensities = zeros(N, 1);
if boundaryType == "airtable"
    iterate = 0;
elseif boundaryType == "annulus"
    iterate =1;
end
for rounding =1:iterate
for n=1:N
    

    if particle(n).z ==1 
        
        contactIntensities(n) = particle(n).contactIs ;
        noncontactcandidates = particle(n).nonContact;
        pop = zeros(length(noncontactcandidates), 1);
    end
    if particle(n).z == 1 && particle(n).contactIs > 20 &&particle(n).edge == 0 % based on high intensity,look for high intensity at ohter location
        y1 = particle(n).y;
        x1 = particle(n).x;
        r1 = particle(n).r;
           
        %for debugging
%         croppedImg = (Gimgfine(uint16(y1-r1-padding):uint16(y1+r1+padding),uint16(x1-r1-padding):uint16(x1+r1+padding)));
%         figure1 = figure;
        
%         subplot(2,1,1, 'Parent', figure1)
%         imshow(Gimgfine(uint16(y1-r1-50):uint16(y1+r1+50),uint16(x1-r1-100):uint16(x1+r1+100)));
%         
%         hold on;
%         plot(r1*cos(particle(n).betas(1))+100+r1,100+r1+r1*sin(particle(n).betas(1)), 'o');
%         title(num2str(particle(n).id))
%         profile = contactfind(croppedImg, x1, y1, r1, n, verbose);
        
        newlocs = particle(n).contactPos; %these are the peaks that we found previously
        newpks = particle(n).contactInt;

        angle = ismembertol(newlocs, particle(n).betas(1), 0.25); %rule out the contact we've detected already
        newneigh = find(angle==0); %potential neighbours are these ones
        
        for k = 1:length(newneigh)
            if newpks(newneigh(k))>0.09 %likely candidate
                %search for particle neighbour
                arrow = newlocs(newneigh(k));
                
                
                for i=1:length(noncontactcandidates)
                ind = noncontactcandidates(i);
                x = particle(ind).x;
                y = particle(ind).y;
                r = particle(ind).r;
                
                nominalAngle = atan2(y-y1, x-x1);
                angle = ismembertol(arrow, nominalAngle, 0.3);
                if angle ==1
                    contactXp2 = x + (r - CR) * cos(atan2(y1-y,x1-x));
                    contactYp2 = y + (r - CR) * sin(atan2(y1-y,x1-x));
                    contactImg = im2double(imcrop(Gimgfine,[contactXp2-CR contactYp2-CR CR*2 CR*2]));
    
                     contactImg = contactImg.*mask;
%                      figure;
%                      imshow(contactImg);
                     contactIp2 = sum(sum(contactImg));
                     [gx,gy] = gradient(contactImg);
                     g2 = (gx.^2 + gy.^2);
                       contactG2p2 = sum(sum(g2));

                     contactXp1 = x1 + (r1 - CR) * cos(atan2(y-y1,x-x1));
                    contactYp1 = y1 + (r1 - CR) * sin(atan2(y-y1,x-x1));
                    contactImg1 = im2double(imcrop(Gimgfine,[contactXp1-CR contactYp1-CR CR*2 CR*2]));
    
                     contactImg1 = contactImg1.*mask;
%                      figure;
%                      imshow(contactImg);
                     contactIp1 = sum(sum(contactImg1));
                     [gx,gy] = gradient(contactImg);
                     g2 = (gx.^2 + gy.^2);
                     contactG2p1 = sum(sum(g2));
                     if contactIp2 > 75 || contactIp1 >75
                        'yep';
                        alreadyneighbours = [particle(n).neighbours];
                        if ~any(alreadyneighbours == ind)
                        particle(n).z= particle(n).z +1; %increase coordination number
                        particle(n).contactG2s(particle(n).z)=contactG2p1; %remember the g2 value of the current contact areas
                        particle(n).contactIs(particle(n).z)=contactIp1;
                        particle(n).color(particle(n).z)='w';
                        particle(n).neighbours(particle(n).z) = particle(ind).id; %particle m is now noted as a neigbour in the particle l datastructure
                        particle(n).betas(particle(n).z) = nominalAngle; %the contact angle to particle m is now noted in the particle l datastructure
                        particle(ind).z= particle(ind).z +1; %increase coordination number
                        particle(ind).contactG2s(particle(ind).z)=contactG2p2; %remember the g2 value of the current contact area
                        particle(ind).contactIs(particle(ind).z)=contactIp2;
                        particle(ind).color(particle(ind).z)='w';
                        particle(ind).neighbours(particle(ind).z) = particle(n).id; %particle m is now noted as a neigbour in the particle l datastructure
                        particle(ind).betas(particle(ind).z) = nominalAngle - pi;
                        pop(i) = i;
                        end
                     end
                end
                end

            end
        end
        

%         elseif particle(n).z ==1 && particle(n).contactIs <4 && particle(n).edge ==0
%             y1 = particle(n).y;
%             x1 = particle(n).x;
%             r1 = particle(n).r;
%         
%              croppedImg = (Gimgfine(uint16(y1-r1-padding):uint16(y1+r1+padding),uint16(x1-r1-padding):uint16(x1+r1+padding)));
%              figure1 = figure;
% 
%              subplot(2,1,1, 'Parent', figure1)
%              imshow(Gimgfine(uint16(y1-r1-50):uint16(y1+r1+50),uint16(x1-r1-50):uint16(x1+r1+50)));
% 
%              hold on;
%              plot(r1*cos(particle(n).betas(1))+50+r1,50+r1+r1*sin(particle(n).betas(1)), 'o');
%              title(num2str(particle(n).id))
%              subplot(2,1,2, 'Parent', figure1)
%         profile = contactfind(croppedImg, x1, y1, r1, n, verbose);
%         newlocs = particle(n).contactPos;
%         newpks = particle(n).contactInt;
%         
%         plot(profile(:,1), profile(:,2));
%         hold on;
%         plot(newlocs, newpks, 'o');
%         
%         angle = ismembertol(newlocs, particle(n).betas(1), 0.25);
%         newneigh = find(angle==0);
%         drawnow
    if ~isempty(pop)
    index = particle(n).nonContact(nonzeros(pop)); %check if n also is inside of the neighbour particles non-Contact
    for t=1:length(index)
        neighbourind = find(particle(index(t)).nonContact ==n);
        recip = nonzeros(neighbourind);
        particle(index(t)).nonContact(recip)=[];
    end
    particle(n).nonContact(nonzeros(pop)) = [];
    
    end
    end
    

end
end
   
if verbose
h3 = figure(20);
hAx1 = subplot(1,1,1,'Parent', h3);
imshow(Gimgfine, 'Parent', hAx1);
hold (hAx1, 'on');
    for n = 1:length(particle)
        particle(n).id
        %viscircles([particle(n).x, particle(n).y], particle(n).r, 'EdgeColor', particle(n).color);
        z = particle(n).z; %get particle coordination number
        if (z>0) %if the particle does have contacts
            for m = 1:z %for each contact
                %draw contact lines
                lineX(1)=particle(n).x;
                lineY(1)=particle(n).y;
                lineX(2) = lineX(1) + particle(n).r * cos(particle(n).betas(m));
                lineY(2) = lineY(1) + particle(n).r * sin(particle(n).betas(m));
%                 cX = lineX(1) + (particle(n).r-CR) * cos(particle(n).betas(m));
%                 cY = lineY(1) + (particle(n).r-CR) * sin(particle(n).betas(m));
                plot(hAx1, lineX, lineY,particle(n).color(m),'LineWidth',2);
                %text(hAx1,lineX(1),lineY(1),num2str(n),'Color','r')
            end
        end
        text(hAx1, particle(n).x, particle(n).y, num2str(n), 'Color', 'r')
    end
drawnow;
end
[num2str(sum([particle.z])), 'detected'   ]
save([directory, files(imgnumb).name(1:end-4),'_preprocessing.mat'],'particle')
    end
   
end

function contactG2 = gradientcalculator(imgchunk)
    [gx,gy] = gradient(imgchunk);
    g2 = (gx.^2 + gy.^2);
    contactG2 = sum(sum(g2));
% contactIntensities = nonzeros(contactIntensities);
% figure;
% histogram(contactIntensities,100);
% figure;
% contactnumbers = [particle.z];
% histogram(contactnumbers)
