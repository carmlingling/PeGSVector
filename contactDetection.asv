%function contactDetection(directory, fileNames, boundaryType,frameidind,verbose)
% %UNTITLED4 Summary of this function goes here
% %   Detailed explanation goes here
% 
% 
directory = './testdata/';
%topDirectory = '/Users/carmenlee/Desktop/20150731reprocesseduniaxial/'
% %topDirectory = './DATA/test/Step09/'
fileNames = '100Hz*.tif'; %image format and regex
frameidind = 15;
%

boundaryType = "annulus"; %if airtable use "airtable" if annulus use "annulus"
radiusRange = [40, 57];
%radiusRange = [45, 78]; %airtable

verbose = true;
% % % % % % % % % 

%% calibrating values
calibrate = false;
global particleNumber1 particleNumber2
particleNumber1 = 269;
particleNumber2 = 1247;

warning('off','signal:findpeaks:largeMinPeakHeight')

%% thresholding limits
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
    imadjust_limits = [0, 1];
    fineimadjust_limits = [0, 100/255];
    rednormal = 2;
    padding =1;
    sigma = 80; % ephraim
elseif boundaryType == "annulus"
    directory = [directory, 'warpedimg/'];
    minpeakheight = 0.10;
    minpeakprominence = 0.02;
    minpeakprom_main = 0.025;
    pxPerMeter = 0.015/939;
    fsigma = 141; %photoelastic stress coefficient
    g2cal = 145; %Calibration Value for the g^2 method, can be computed by joG2cal.m
    dtol = 30; % How far away can the outlines of 2 particles be to still be considered Neighbours
    contactG2Threshold = 0.5; %sum of g2 in a contact area larger than this determines a valid contact
    CR = 15; 
    imadjust_limits = [0, 0.6];
    fineimadjust_limits = [0/255, 30/255];%[13/255, 39/255]
    rednormal = 8;
    padding = 1;
    sigma  = 50; %for blurring large scale features
    polarizerstrip = [[2731,2719,3643,3666];[212,212,6099,6100]];
end

%% importing files

[directory, fileNames]
files = dir([directory, fileNames(1:end-4), '.tif']);
centersfile = dir([directory, 'centers_tracked.txt']);
pData = readmatrix([directory,centersfile.name],"NumHeaderLines", 1); %Read Position data from centers file


%% setting up mask



mask = abs(-CR:CR);
mask = mask.^2 + mask.^2';
maskCR = double(sqrt(mask) <= CR-1);





%% image manipulation
for imgnumb = 1:1
    
    clear particle %reinitialize the particle structure

    %read in image
    Img = imread([directory,files(imgnumb).name]);
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
        if verbose
            figure;
            imshow(Gimgp)
            viscircles([pData(data,3), pData(data,4)], round(pData(data,5)))
            title('Gimgp')
           
        end

    elseif boundaryType == "annulus"
        
        bckgnd = poly2mask(polarizerstrip(1,:),polarizerstrip(2,:), length(Gimg), length(Gimg));
        Gimg = inpaintCoherent(Gimg,bckgnd,'SmoothingFactor',5,'Radius',15);
        
        if calibrate == true
           figure;
           imshow(Gimg);
           hold on
           plot([2731,2719,3643,3666, 2731],[212,212,6099,6100, 212],'b','LineWidth',1)
           drawnow;
        end
        
        Gimg=imsubtract(Gimg,Rimg./rednormal);
   
        G = fspecial('gaussian', 3*sigma+1, sigma);
        yb = imfilter(imcomplement(Rimg), G, 'replicate');
        Gimg = bsxfun(@minus, Gimg,yb*.09);

        Gimg= im2double(Gimg);
        Gimg = Gimg.*(Gimg > 0);
    
    
        Gimgd = imadjust(Gimg,imadjust_limits); %regular contrast
    
        Gimgfine = imadjust(Gimg, fineimadjust_limits); %super boosted contrast
    

    end
   
%% initialize data structure
    
    frame = str2double(files(imgnumb).name(frameidind:frameidind+3))
    
    data = find(pData(:,1) == frame);%from the centers information, find the particles that are in a the current frame
    
    if ~isempty(data)

        N = size(data,1);

        particle(1:N) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[], 'nonContact', [], 'contactPos', [], 'contactInt', [], 'edge', 0);
        for n = 1:N %Bookkeeping from centers-tracked
            particle(n).id= pData(data(n),2);
            particle(n).x = pData(data(n),3); 
            particle(n).y = pData(data(n),4); 
            particle(n).r = round(pData(data(n),5));
            particle(n).edge = pData(data(n), 6);
            particle(n).rm = particle(n).r*pxPerMeter;
            particle(n).fsigma = fsigma;
        end   
        
        
        for n=1:N %loop over particles
    
    %create a circular mask
  
            r = particle(n).r;
            if round(particle(n).y+r)<size(Gimg, 1)&&round(particle(n).x+r)<size(Gimg,2)&&round(particle(n).y-r)>1&&round(particle(n).x-r)>1 %double check to make sure the bounds are within the image

                mask = abs(-r:r);
                mask = mask.^2 + mask.^2';
                mask1 = (sqrt(mask) <= r);

                %This crops out a particle
                cropXstart = round(particle(n).x-r);
                cropXstop = round(particle(n).x-r)+ size(mask1,1)-1;
                cropYstart = round(particle(n).y-r);
                cropYstop = round(particle(n).y-r)+ size(mask1,2)-1;


                particleImg= Gimgd(cropYstart:cropYstop, cropXstart:cropXstop).*mask1;
                particle(n).forceImage=particleImg; %save this so we can fit to this image later in diskSolve

                %create a circular mask with a radius that is one pixel smaller
                %for cropping out the relevant gradient

                mask2 = double(sqrt(mask) <= r-1);

                %Compute G^2 for each particle
                [gx,gy] = gradient(particleImg);
                g2 = (gx.^2 + gy.^2).*mask2;
                particle(n).g2 = sum(sum(g2));
                particle(n).f = particle(n).g2/g2cal; %saving some particle scale features
            else
                error('badimage!!')

            end
        end
%% look at neighbours

        xmat = pData(data,3);
        ymat = pData(data,4);
        rmat = pData(data,5);
        
        rmats = rmat; %Saves our radius matrix for later
        
        dmat = pdist2([xmat,ymat],[xmat,ymat]); %Creates a distance matrix for particle center locations
        rmat = rmat + rmat'; %Makes a combination of radii for each particle
        
        friendmat = dmat < (rmat + dtol) & dmat~=0; %Logical "friend" matrix
        
        friendmat = triu(friendmat); %Only examine the upper triangle portion (no repeats)
        [f1, f2] = find(friendmat == 1); %Creates an index of particles that are considered touching
%% 
        xpairs = [xmat(f1),xmat(f2)];
        ypairs = [ymat(f1),ymat(f2)];
        rpairs = [rmats(f1),rmats(f2)];
      
%% loop over friends
    
        for l = 1:length(f1)
            
            x = xpairs(l,:);
            y = ypairs(l,:);
            r = rpairs(l,:);
            [contactG2p, contactIp] = contactspot(x,y,r, CR, Gimgd, maskCR);
            
            if(contactG2p(1) > contactG2Threshold && contactG2p(2) > contactG2Threshold)
        
                %this is a valid contact, remember it
                particle(f1(l)).z= particle(f1(l)).z +1; %increase coordination number
                particle(f1(l)).contactG2s(particle(f1(l)).z)=contactG2p(1); %remember the g2 value of the current contact area
                particle(f1(l)).contactIs(particle(f1(l)).z)=contactIp(1); %changes to color
                particle(f1(l)).color(particle(f1(l)).z)='r'; %changes to color
                particle(f1(l)).neighbours(particle(f1(l)).z) = particle(f2(l)).id; %particle m is now noted as a neigbour in the particle l datastructure
                particle(f1(l)).betas(particle(f1(l)).z) = atan2(y(2)-y(1),x(2)-x(1)); %the contact angle to particle m is now noted in the particle l datastructure
                particle(f2(l)).z= particle(f2(l)).z+1; %increase coordination number
                particle(f2(l)).contactG2s(particle(f2(l)).z)=contactG2p(2); %remember the g2 value of the current contact area
                particle(f2(l)).contactIs(particle(f2(l)).z)=contactIp(2);
                particle(f2(l)).color(particle(f2(l)).z)='r'; %changes to color
                particle(f2(l)).neighbours(particle(f2(l)).z) = particle(f1(l)).id; %particle m is now noted as a neigbour in the particle l datastructure
                particle(f2(l)).betas(particle(f2(l)).z) = atan2(y(1)-y(2),x(1)-x(2));
  
            else %we try the more refined method of contact detection
       
            %find peaks in intensity for each particle, record the value of
            %the peak and the angular location relative to the x axis of
            %each particle
            [pks, locs] = peakfinder(x(1), y(1), r(1), f1(l), Gimgfine, minpeakheight, minpeakprominence, minpeakprom_main, padding);
            [pks2, locs2] = peakfinder(x(2), y(2), r(2), f2(l), Gimgfine, minpeakheight, minpeakprominence, minpeakprom_main, padding);
    
            %compare the locations to whatever the nominal angle is (center
            %to center)
            
            nominalAngle = atan2(y(2)-y(1), x(2)-x(1));
            nominalAngleEXP = [nominalAngle-2*pi, nominalAngle, nominalAngle+2*pi];

            [angle] = ismembertol(nominalAngleEXP,locs',  pi/6,'DataScale', 1);

            nominalAngle2 = atan2(y(1)-y(2), x(1)-x(2));
            nominalAngleEXP = [nominalAngle2-2*pi, nominalAngle2, nominalAngle2+2*pi];

            [angle2] = ismembertol(nominalAngleEXP,locs2',  pi/6,'DataScale', 1);
        
%             if f1(l) == particleNumber1 && f2(l) == particleNumber2
%                 error()
%             end
            if ~isempty(angle(angle==1)) && ~isempty(angle2(angle2==1))

                particle(f1(l)).z= particle(f1(l)).z +1; %increase coordination number
                particle(f1(l)).contactG2s(particle(f1(l)).z)=contactG2p(1); %remember the g2 value of the current contact area
                particle(f1(l)).contactIs(particle(f1(l)).z)=contactIp(1);
                particle(f1(l)).color(particle(f1(l)).z)='y';
                particle(f1(l)).neighbours(particle(f1(l)).z) = particle(f2(l)).id; %particle m is now noted as a neigbour in the particle l datastructure
                particle(f1(l)).betas(particle(f1(l)).z) = nominalAngle; %the contact angle to particle m is now noted in the particle l datastructure
                particle(f2(l)).z= particle(f2(l)).z +1; %increase coordination number
                particle(f2(l)).contactG2s(particle(f2(l)).z)=contactG2p(2); %remember the g2 value of the current contact area
                particle(f2(l)).contactIs(particle(f2(l)).z)=contactIp(2);
                particle(f2(l)).color(particle(f2(l)).z)='y';
                particle(f2(l)).neighbours(particle(f2(l)).z) = particle(f1(l)).id; %particle m is now noted as a neigbour in the particle l datastructure
                particle(f2(l)).betas(particle(f2(l)).z) = nominalAngle2;
            
            else
                s = length(particle(f1(l)).nonContact);
                w = length(particle(f2(l)).nonContact);
                particle(f1(l)).nonContact(s+1) = f2(l);
                particle(f2(l)).nonContact(w+1) = f1(l);
            end
        
        
        
        
            if isempty(particle(f1(l)).contactPos)
                particle(f1(l)).contactPos = locs;
                particle(f1(l)).contactInt = pks;
            end
            if isempty(particle(f2(l)).contactPos)
                particle(f2(l)).contactPos = locs2;
                particle(f2(l)).contactInt = pks2;
            end
        
            end
    

    
    
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
        [pks, locs] = peakfinder(x, y, r, Gimgfine, minpeakheight, minpeakprominence,minpeakprom_main, padding )
        

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
    
                     contactImg = contactImg.*maskCR;
%                      figure;
%                      imshow(contactImg);
                     contactIp2 = sum(sum(contactImg));
                     [gx,gy] = gradient(contactImg);
                     g2 = (gx.^2 + gy.^2);
                       contactG2p2 = sum(sum(g2));

                     contactXp1 = x1 + (r1 - CR) * cos(atan2(y-y1,x-x1));
                    contactYp1 = y1 + (r1 - CR) * sin(atan2(y-y1,x-x1));
                    contactImg1 = im2double(imcrop(Gimgfine,[contactXp1-CR contactYp1-CR CR*2 CR*2]));
    
                     contactImg1 = contactImg1.*maskCR;
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
        particle(n).id;
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
%end
function contactG2 = gradientcalculator(imgchunk)
    [gx,gy] = gradient(imgchunk);
    g2 = (gx.^2 + gy.^2);
    contactG2 = sum(sum(g2));

end


function [contactG2p, contactIp]=contactspot(x, y, r, CR, Gimgd, maskCR)
    contactangle = [atan2(y(2)-y(1),x(2)-x(1)), atan2(y(1)-y(2), x(1)-x(2))];
    contactXp = round(x + (r -  1 - CR).* cos(contactangle));
    contactYp = round(y + (r -1- CR).* sin(contactangle));
    
    contactImg = im2double(imcrop(Gimgd,[contactXp(1)-CR contactYp(1)-CR CR*2 CR*2]));
    contactImg = contactImg.*maskCR;
    
    contactG2p = [gradientcalculator(contactImg)];
    contactIp = [sum(sum(contactImg))];
    
    contactImg = im2double(imcrop(Gimgd,[contactXp(2)-CR contactYp(2)-CR CR*2 CR*2]));
    contactImg = contactImg.*maskCR;
    contactG2p(2,:)= gradientcalculator(contactImg);
    contactIp(2,:) = sum(sum(contactImg));
    %contactG2p = [G1 G2]
    
    
end

function [profile] = contactfind(croppedImg, r)
    
    [a,b] = size(croppedImg);
    [X, Y] = meshgrid( (1:b)-r, (1:a)-r);
    R = sqrt(X.^2 + Y.^2);
    m2 = (2*r/3<R&R<r-1);

    %maskedImg = double(croppedImg).*m2; %if you want to see what the
    %masked image looks like
    values = double(croppedImg(m2));


    theta = atan2(Y, X);
    angles = theta(m2);
    combo = [angles, values];
    profile = sortrows(combo);

    profile = smoothdata(profile,1,'sgolay', length(profile)/25);

    profX = [profile(1:uint16(length(profile)/4),1)+(2*pi), profile(1:uint16(length(profile)/4), 2)];
    profile = [profile; profX];
    profile = sortrows(profile);
end

function [finalpks, finallocs] = peakfinder(x, y, r, ind, Gimgfine, minpeakheight, minpeakprominence,minpeakprom_main, padding )
        global particleNumber1 particleNumber2
        croppedImg = (Gimgfine(round(y-r-padding):round(y+r+padding),round(x-r-padding):round(x+r+padding)));
        
        [profile] = contactfind(croppedImg, r-1);
        
        if any(profile(:,1)>minpeakheight)
            [pkints, locints] = findpeaks(profile(:,2),profile(:,1), "MinPeakHeight", minpeakheight,'MinPeakProminence', minpeakprominence, 'MaxPeakWidth', pi);
        else
            pkints = [];
            locints = [];
        end
        if any(profile(:,1)>minpeakprom_main)
            [pks,locs] =findpeaks(profile(:,2), profile(:,1),'MinPeakProminence',minpeakprom_main,'MaxPeakWidth', pi);
        else
            pks = [];
            locs = [];
        end
        pks = [pks;pkints];
        locs = [locs;locints];
        inds = find(locs>(pi));
        locs(inds) = locs(inds)-2*pi;

        [a, in] = uniquetol(locs, pi/12);
        [b] = uniquetol(locs, pi/12,'highest');
        finallocs = mean([a, b], 2);
        finalpks = pks(in);

        if ind == particleNumber1 || ind == particleNumber2
            figure1 = figure;
            subplot(2,1,1, 'Parent', figure1)

       
        
            imshow(croppedImg);
            title(num2str(ind));
        
            hold on;
            contactLoc2 = [r+r.*cos(locs), r+r.*sin(locs)];
            contactLoc = [r+r.*cos(finallocs), r+r.*sin(finallocs)];
            %plot(r+r.*cos(nominalAngle), r+r.*sin(nominalAngle-pi),'yo')
            if ~isempty(contactLoc2)
            plot(contactLoc2(:,1), contactLoc2(:,2), 'ro');
            plot(contactLoc(:,1), contactLoc(:,2), 'bo');
            plot(r, r, 'gx')
            end
            subplot(2,1,2, 'Parent', figure1);
            plot(profile(:,1), profile(:,2));
            title(num2str(ind))
            hold on;
            plot(locs, pks, 'o');
            plot(finallocs, pks(1:length(finallocs)), 'o');
        
        
        end
        
%         subplot(2,3,6, 'Parent', figure1)
%         plot(c);
%         hold on;
%         plot(loc, peaks, 'o');
        
end
% function calibrateparameters()

