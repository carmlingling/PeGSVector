%rebalancing forces for particles that have a fit error greater than a given
%threshold
function newtonize(directory, fileNames, boundaryType, verbose)

if boundaryType == "annulus"

%directory = '~/Desktop/230502_2/warpedimg/'
filename = '*solved.mat';
imname = '*Synth.jpg';
imagefiles = dir([directory,'/synthImg/', imname]);
%directory = [directory,'warpedimg/'];
forcefiles = dir([directory,'solved/', filename]);


fiterrorcutoff = 1850;
cutoff = 1850;
%troubleid = [1700,1760,1765, 2066, 2039];
maskradius = 0.96/2;
troubleid = []
elseif boundaryType == "airtable"
forcefiles = dir([directory,'solved/', fileNames(1:end-4),'_solved.mat']);
imagefiles = dir([directory,'synthImg/', fileNames(1:end-4), '-Synth.jpg']);
cutoff = 250;
maskradius = 0.96/2;
troubleid = [];
end

fig1 =figure(1);
hAx1 = subplot(1,1,1,'Parent', fig1);
fig2 =figure(2);
hAx2 = subplot(1,1,1,'Parent', fig2);
for frame=1:length(forcefiles)
    %frame=frame+105
    pres = load([forcefiles(frame).folder,'/',forcefiles(frame).name]);
    pres = pres.pres;
    id2ind = [pres.id];
    for k=1:length(troubleid)
        new = find(id2ind == troubleid(k));
    
        if ~isempty(new)
            pres(new).edge = -1;
        end
    
    end
    N = length(pres);
    errors = zeros(N,1);
    for n=1:N
        if ~isempty(pres(n).fitError) 
        errors(n) = pres(n).fitError;
        end
    end
    if boundaryType == "annulus"
    othererrors = zeros(N,1);
    
    for n=1:N
       
        if ~isempty(pres(n).fitError)
            
        im1 = pres(n).synthImg;
        im2 = pres(n).forceImage;
        if size(im1, 1) ~= size(im2, 1)
            im1 = im1(1:end-1, 1:end-1);
        end

        residual = imsubtract(im1, im2);
        residual = abs(residual);
        int = sum(sum(residual));
        othererrors(n) = int;
        pres(n).fitError=int;
        end
    end
    bigerr = find(othererrors>fiterrorcutoff); %gives index
    bigothererr = find(othererrors>cutoff);
    elseif boundaryType =="airtable"
        bigerr = find(errors>cutoff);
    end
    if verbose
    figure;
    histogram(errors);
    title('errors')
    end
%     figure;
%     histogram(othererrors)
%     title('othererrors')
    
    
    edges = find([pres.edge]~=0);



    I = imread([imagefiles(frame).folder,'/', imagefiles(frame).name]);
    
    imshow(I, 'Parent', hAx1);
    x = [pres.x];
    y = [pres.y];
    r = [pres.r];
    x2 = x(bigerr)';
    y2 = y(bigerr)';
    r2 = r(bigerr)';
%     x3 = x(bigothererr)';
%     y3 = y(bigothererr)';
%     r3 = r(bigothererr)';
    x4 = x(edges)';
    y4 = y(edges)';
    r4 = r(edges)';
    viscircles(hAx1, [x2, y2], r2, 'EdgeColor', 'y', 'Linewidth', 2);
    %viscircles(ax, [x3, y3], r3, 'EdgeColor', 'b', 'Linewidth', 1);
    viscircles(hAx1, [x4, y4], r4, 'EdgeColor', 'r', 'Linewidth', 1);
    for n=1:N
        text(hAx1, x(n), y(n), num2str(pres(n).id), 'Color', 'b')
    end
    axis on;
    drawnow;

    


%%
%sort the bad particles by the number of bad neighbours/edges they have,
%fewest first
id2ind = [pres.id];
rankerror = zeros(length(bigerr),1);
for badparticle = 1:length(bigerr)
    %bigerr(badparticle) is an index
    neighbours = pres(bigerr(badparticle)).neighbours; %id
    nonedge = neighbours(neighbours>0);
    edge = neighbours(neighbours<0);
    %id2ind(nonedge)
    badneighbours = ismember(nonedge, id2ind(bigerr));
    rankerror(badparticle) = sum(badneighbours)+length(edge);
end
[rankerror,sortIdx] = sort(rankerror);
bigerr = bigerr(sortIdx);
        


% if boundaryType == "annulus"
%            iterate = 0;
%         
%         maxiter = length(bigerr);
%         while iterate <maxiter & ~isempty(bigerr)
%     for n=1:N
%         if pres(n).edge >= 0 & errors(n)>cutoff %check to only work with non-inner edge particles and those above the cutoff
%                 IDN = pres(n).id; %get the id of the particle in question
%                 Nneighbours = [pres(n).neighbours]; %find the neighbours (ids)
% 
%                 z = pres(n).z; %number of neighbours
%                 badneighbour = zeros(z, 1); %make an array for potential bad neighbours that are also above the cutoff
% 
% 
%                 for m=1:z %loop over the neighbours
% 
%                     indN = find(id2ind==Nneighbours(m)); %find the index of the neighbour
%             		if length(indN)>1
%                 		indN = indN(1); %sometimes the contact finding algorithm finds the same neighbour multiple times this is a needed bug to fix but this should work for now
%             		end
%                     if ~ismember(indN, bigerr) && pres(indN).edge >=0 %make sure it's not bad either
%                         neighbours = [pres(indN).neighbours];
%                         if length(neighbours)>1
%                             positionNeighbour = find(neighbours == IDN);
%                             if length(positionNeighbour)>1
%                 				positionNeighbour=positionNeighbour(1);
%             				end
%                         elseif neighbours == IDN;
%                             positionNeighbour = 1;
%                         end
%                         
%                         pres(n).forces(m) = pres(indN).forces(positionNeighbour);
%                         pres(n).alphas(m) = -pres(indN).alphas(positionNeighbour);
%                         pres(indN).alphas(positionNeighbour)
%                     else
%                         badneighbour(m) = indN;
%                     end
%                 end
% 
%                 force = pres(n).forces;
%                 alpha =pres(n).alphas;
%                 betas = pres(n).betas;
%                 fsigma = pres(n).fsigma;
%                 rm = pres(n).rm;
%                 scaling = 0.5;
%                 verbose = false;
%                 template= pres(n).forceImage;
%                 template = imadjust(particle(n).forceImage);
%                 template = imresize(template,scaling);
%                 if verbose
%                     subplot(1,2,1)
%                     imshow(template)
%                     hold on;
%                     plot(pres(n).r*(1+cos(betas)), pres(n).r*(1+sin(betas)), 'o');
%                     for m=1:z
%                         plot([pres(n).r*(1+cos(betas(m))),pres(n).r*(1+cos(pres(n).betas(m)+alphas(m)))] , [pres(n).r*(1+sin(pres(n).betas(m))),pres(n).r*(1+sin(pres(n).betas(m)+sin(pres(n).alphas(m))))], '-')
%                     end
%                 end
%                 
% 
%                 size of the force image
%                 px = size(template,1);
%                 badneighbour=nonzeros(badneighbour);
%                 beta = -betas+pi/2;
% 
%                 if length(badneighbour)==1
% 
%                     sum1 = 0;
%                     sum2 = 0;
%                     for k = 1:z
%                         if(Nneighbours(k)~=badneighbour)
% 
%                             sum1 = sum1 + force(k)*sin(alpha(k)+beta(k)); %xforces
%                             sum2 = sum2 + force(k)*cos(alpha(k)+beta(k)); %yforces
% 
%                         end
%                     end
%                     f = sqrt(sum1^2+sum2^2);
% 
%                     a = asin(-sum1/f);
%                     loc = Nneighbours == badneighbour;
%                     pres(n).forces(loc) = f;
%                     pres(n).alphas(loc) = -a;
% 
%                     img = joForceImg (z, pres(n).forces, pres(n).alphas, betas, fsigma, rm, px, verbose);
%                     bigerr(bigerr ==n)=[];
%                     pres(n).synthImg = img;
%                     pres(n).fitError = sum(sum(abs(imsubtract(img, template))));
%                 elseif isempty(badneighbour)
%                     img = joForceImg (z, pres(n).forces, pres(n).alphas, betas, fsigma, rm, px, verbose);
%                     bigerr(bigerr ==n) =[];
%                   
%                     pres(n).synthImg = img;
%                     
%                     err = abs(sum(sum( ( c_mask.*(template-img).^2) )))
%                     pres(n).fitError = sum(sum(abs(imsubtract(img, template))));
%                 elseif length(badneighbour)>1
%                     loc = Nneighbours == badneighbour;
%                     pres(n).forces(loc) = f;
%                     pres(n).alphas(loc) = -a;
%                 end
%                 if maxbad < length(badneighbour)
%                     maxbad = length(badneighbour);
%                 end
%             end
%     end
%         end

        
%   if boundaryType =="airtable" | "annulus"
        %sort bad particles by the number of bad neighbours
        iterate = 0;
        
        maxiter = length(bigerr);
        while iterate <maxiter & ~isempty(bigerr)
            
            badparticle = 1; %check only particles above the cutoff
            n = bigerr(badparticle); %index of the particle in question in the particle structure
            IDN = pres(n).id;        %ID of the particle in question


            %load some stuff from the particle structure
            fsigma = pres(n).fsigma;
            rm = pres(n).rm;

            verbose = false;
            template= pres(n).forceImage;
            px = size(template,1);

            rank = rankerror(badparticle);
            
            neighbours = pres(n).neighbours;
                if rank>0
                nonedge=neighbours((neighbours>0));
                badneighbours = find(ismember(nonedge, id2ind(bigerr)));
                
                bad2 = find(neighbours < 0);
                if ~isempty(badneighbours) & ~isempty(bad2)
                badneighbours = [badneighbours,bad2];
                elseif isempty(badneighbours) & ~isempty(bad2)
                    badneighbours = bad2;
                elseif ~isempty(badneighbours) & isempty(bad2)
                    badneighbours = badneighbours;
                end
                good = 1:length(neighbours);
                good =good(setdiff(1:end,badneighbours));
                %good = find(neighbours~=neighbours(badneighbours))
                else
                    good = 1:length(neighbours);
                end
                for m=1:length(good) %loop over the neighbours, 
                    %now find the force from the neighbour particles' list
                    %of forces
                    k = good(m);
                    indexN = find(id2ind==neighbours(k));
                    positionNeighbour = find(pres(indexN).neighbours==IDN);
                    if length(positionNeighbour)>1
                        positionNeighbour = positionNeighbour(1)
                    end
                    pres(n).forces(k) = pres(indexN).forces(positionNeighbour);
                    pres(n).alphas(k) = pres(indexN).alphas(positionNeighbour);
                    
                end
                
                force = pres(n).forces;
                alpha =pres(n).alphas;
                betas = pres(n).betas;
                
            if rankerror(badparticle) == 0
                bigerr(badparticle) =[];
                rankerror(badparticle) = [];
                    
            elseif rankerror(badparticle)==1
                %only one or zero neighbour/edges we cannot match
                'small route'


           
                    sum1 = 0;
                    sum2 = 0;
                    sum3 = 0;
                    
                    for m = 1:length(good)
                        k = good(m);
                            
                            sum1 = sum1 + force(k)*sin(alpha(k)-betas(k)+pi/2); %xforces
                            sum2 = sum2 + force(k)*cos(alpha(k)-betas(k)+pi/2); %yforces
                            sum3 = sum3 + force(k)*sin(alpha(k));
                            
                        
                    end
                    
                    f = sqrt(sum1^2+sum2^2);
                    
                    a = asin(-sum3/f);
                    
                    pres(n).forces(badneighbours) = f;
                    pres(n).alphas(badneighbours) = a;
                    bigerr(badparticle) =[];
                    rankerror(badparticle) = [];
                    
                   

            else
                %now we have more than 1 unconstrained force so we're going
                %to try to fix the forces we know from the neighbours and
                %fit the others
                'big route'
                
                fixedforce =[];
                fixedbeta =[];
                fixedalpha =[];
                for m = 1:length(good)
                    fixedforce = [fixedforce,pres(n).forces(good(m))];
                    fixedalpha = [fixedalpha,pres(n).alphas(good(m))];
                    fixedbeta = [fixedbeta, pres(n).betas(good(m))];
                end
                freeforce =[];
                freebeta =[];
                freealpha =[];
                for m = 1:length(badneighbours)
                    freeforce = [freeforce,pres(n).forces(badneighbours(m))];
                    freealpha = [freealpha,pres(n).alphas(badneighbours(m))];
                    freebeta = [freebeta, pres(n).betas(badneighbours(m))];
%                         pres(n).alphas(m) = -pres(indN).alphas(positionNeighbour);
%                         pres(indN).alphas(positionNeighbour);
                end
                    cx=px/2;cy=px/2;ix=px;iy=px;r=maskradius*px;
                    [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
                    c_mask=((x.^2+y.^2)<=r^2);   
                    z = pres(n).z;
                    q = length(freeforce);
                    length(fixedforce);
                    if q>length(fixedforce)
                        func = @(par) joForceImgFixed(z, par(1:q),par(q+1:q+q), freebeta(1:q), fsigma, rm, px, verbose, fixedforce', fixedalpha', fixedbeta); %+par(2*z+1); %this is the function I want to fit (i.e. synthetic stres image), the fitting paramters are in vector par
                    else
                        func = @(par) joForceImgFixed(z, par(1:q),par(q+1:q+q), freebeta(1:q), fsigma, rm, px, verbose, fixedforce', fixedalpha', fixedbeta, false);
                    end
                        %This is the error function we are actually fitting,
                    %that is, the distance between our fit function and the
                    %real particle image in terms of the sum of squares of the pixelwise differnce.
                    %Also a mask is applied to crop
                    %out only the circular part of the particle. 
                    
                    err = @(par) abs(sum(sum( ( c_mask.*(template-func(par)).^2) ))); %BUG: for some reason I sometimes get imaginary results, this should not happen

                    %Set up initial guesses
                    p0 = zeros(2*q, 1);
                    p0(1:q) = freeforce;
                    p0(q+1:2*q) = freealpha;

                    %Do the fit, will also work with other solvers
                    %TODO: make a user defined option to select between
                    %different solvers
                    fitoptions = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxIter',200,'MaxFunEvals',400,'TolFun',0.01,'Display','final-detailed');
                    
                    p=lsqnonlin(err,p0,[],[],fitoptions);

                    %get back the result from fitting
                    forcesfit = p(1:q);
                    alphasfit = p(q+1:q+q);

                    %resudual
                    fitError = err(p);
                    length(forcesfit);
                    for m = 1:length(forcesfit)
                        pres(n).forces(badneighbours(m))= forcesfit(m);
                        pres(n).alphas(badneighbours(m))= alphasfit(m);
                    end
                    bigerr(badparticle) = [];
                    rankerror(badparticle) = [];
        end
            z = pres(n).z;
            img = joForceImg(z, pres(n).forces, pres(n).alphas, betas, fsigma, rm, px, verbose, false);
             pres(n).synthImg = img;     
                
                rankerror = zeros(length(bigerr),1);
            for badparticle = 1:length(bigerr)
                    neighbours = pres(bigerr(badparticle)).neighbours;
                    nonedge = neighbours(neighbours>0);
                    edge = neighbours(neighbours<0);
                    badneighbours = ismember(nonedge, id2ind(bigerr));
                    rankerror(badparticle) = sum(badneighbours)+length(edge);
            end
            [rankerror,sortIdx] = sort(rankerror);
            bigerr = bigerr(sortIdx);
            iterate= iterate+1;
        end
   %     end
    
%                 for m=1:z %loop over the neighbours
% 
%                     indN = find(id2ind==Nneighbours(m)) %find the index of the neighbour
%             		if length(indN)>1
%                 		indN = indN(1); %sometimes the contact finding algorithm finds the same neighbour multiple times this is a needed bug to fix but this should work for now
%             		end
%                     if ~ismember(indN, bigerr) && pres(indN).edge >=0 %make sure it's not bad either
%                         neighbours = [pres(indN).neighbours];
%                         if length(neighbours)>1
%                             positionNeighbour = find(neighbours == IDN);
%                             if length(positionNeighbour)>1
%                 				positionNeighbour=positionNeighbour(1);
%             				end
%                         elseif neighbours == IDN;
%                             positionNeighbour = 1;
%                         end
%                         
%                         pres(n).forces(m) = pres(indN).forces(positionNeighbour);
%                         pres(n).alphas(m) = -pres(indN).alphas(positionNeighbour);
%                         pres(indN).alphas(positionNeighbour);
%                     else
%                         badneighbour(m) = indN;
%                     end
%                 end
% 
%                 force = pres(n).forces;
%                 alpha =pres(n).alphas;
%                 betas = pres(n).betas;
%                 fsigma = pres(n).fsigma;
%                 rm = pres(n).rm;
%                 scaling = 0.5;
%                 verbose = false;
%                 template= pres(n).forceImage;
%                 %template = imadjust(particle(n).forceImage);
%                 %template = imresize(template,scaling);
%                 if verbose
%                     subplot(1,2,1)
%                     imshow(template)
%                     hold on;
%                     plot(pres(n).r*(1+cos(betas)), pres(n).r*(1+sin(betas)), 'o');
%                     for m=1:z
%                         plot([pres(n).r*(1+cos(betas(m))),pres(n).r*(1+cos(pres(n).betas(m)+alphas(m)))] , [pres(n).r*(1+sin(pres(n).betas(m))),pres(n).r*(1+sin(pres(n).betas(m)+sin(pres(n).alphas(m))))], '-')
%                     end
%                 end
%                 %
% 
%                 %size of the force image
%                 px = size(template,1);
%                 badneighbour=nonzeros(badneighbour);
%                 beta = -betas+pi/2;
% 
%                 if length(badneighbour)==1
% 
%                     sum1 = 0;
%                     sum2 = 0;
%                     for k = 1:z
%                         if(Nneighbours(k)~=badneighbour)
% 
%                             sum1 = sum1 + force(k)*sin(alpha(k)+beta(k)); %xforces
%                             sum2 = sum2 + force(k)*cos(alpha(k)+beta(k)); %yforces
% 
%                         end
%                     end
%                     f = sqrt(sum1^2+sum2^2);
% 
%                     a = asin(-sum1/f);
%                     loc = Nneighbours == badneighbour;
%                     pres(n).forces(loc) = f;
%                     pres(n).alphas(loc) = -a;
% 
%                     img = joForceImg (z, pres(n).forces, pres(n).alphas, betas, fsigma, rm, px, verbose);
%                     bigerr(bigerr ==n)=[];
%                     pres(n).synthImg = img;
%                     pres(n).fitError = sum(sum(abs(imsubtract(img, template))));
%                 elseif isempty(badneighbour)
%                     img = joForceImg (z, pres(n).forces, pres(n).alphas, betas, fsigma, rm, px, verbose);
%                     bigerr(bigerr ==n) =[];
%                     %             figure
%                     %             imshow(template)
%                     %             hold on;
%                     %             for m=1:length(force)
%                     %             lineX(1)=px/2;
%                     %                 lineY(1)=px/2;
%                     %                 lineX(2) = lineX(1) + 100*force(m)*cos(alpha(m));
%                     %                 lineY(2) = lineY(1) + 100*force(m)*sin(alpha(m));
%                     %             plot(lineX, lineY)
%                     %             end
%                     pres(n).synthImg = img;
%                     
%                     %err = abs(sum(sum( ( c_mask.*(template-img).^2) )))
%                     pres(n).fitError = sum(sum(abs(imsubtract(img, template))));
% %                 elseif length(badneighbour)>1
% %                     loc = Nneighbours == badneighbour;
% %                     pres(n).forces(loc) = f;
% %                     pres(n).alphas(loc) = -a;
%                 end
%                 if maxbad < length(badneighbour)
%                     maxbad = length(badneighbour);
%                 end
%             end


%         end
%     end
%     
%  iteration = iteration+1;
% end
if boundaryType == "annulus"
for l=1:length(edges)
    n = edges(l);
    edgeneighbours = pres(n).neighbours;
    z = pres(n).z;
    IDN = pres(n).id;
    for m=1:z
        Nindex = find(id2ind == edgeneighbours(m));
        if pres(Nindex).edge ==0
            steal = find(pres(Nindex).neighbours == IDN);
            if length(steal) ==1
            pres(n).forces(m) = pres(Nindex).forces(steal);
            pres(n).alphas(m) = -pres(Nindex).alphas(steal);
            elseif length(steal) >1
            pres(n).forces(m) = pres(Nindex).forces(steal(1));
            pres(n).alphas(m) = -pres(Nindex).alphas(steal(1));
            end
        else
            pres(n).forces(m) = 0;
            pres(n).alphas(m) = 0;
        
        end
    end
    template= pres(n).forceImage;
    px = size(template,1);
    
    img = joForceImg (z, pres(n).forces, pres(n).alphas, pres(n).betas, pres(n).fsigma, pres(n).rm, px, verbose);
    pres(n).synthImg = img;
end
end
 %%


% %forcebalance
% for n = 1:N
%     error = pres(n).fitError;
%     Nneighbours = pres(n).neighbours;
%     if length(Nneighbours)~=length(unique(Nneighbours))
%         %find non-unique contacts
%         'oop'
%     end
%     z=pres(n).z;
%     IDN = pres(n).id;
%     for m=1:z
%         
%         indi = find(id2ind==Nneighbours(m));
%         if n < indi
%         f1 = pres(n).forces(m);
%         if pres(indi).edge ~=-1
%             steal = find(pres(indi).neighbours == IDN);
%             if length(steal) == 1
%             errorM = pres(indi).fitError;
%             f2 = pres(indi).forces(steal);
%             elseif length(steal) >1
% 		steal = steal(1);
% 		errorM = pres(indi).fitError;
%             	f2 = pres(indi).forces(steal);
%             	end
%             favg = (f1+f2)/2;
%             pres(indi).forces(steal) = favg;
%             pres(n).forces(m) = favg;
%             %contactforces = [contactforces, favg];
%         end
%         end
%     end
%     
%     template = pres(n).forceImage;
%     px = size(template, 1);
%     img = joForceImg (z, pres(n).forces, pres(n).alphas, pres(n).betas, pres(n).fsigma, pres(n).rm, px, verbose);
%     pres(n).synthImg = img;
% end
if boundaryType == "annulus"
    save([directory, 'solved/',forcefiles(frame).name(1:end-4),'_update.mat'],'pres');
else
    save([directory,forcefiles(frame).name(1:end-4),'_update.mat'],'pres');
end
NN = length(pres);

bigSynthImg = zeros(size(I,1),size(I,2)); %make an empty image with the same size as the camera image
    for n=1:NN %for all particles
        %display(['fitting force(s) to particle ',num2str(n)]); %status indicator
        if (pres(n).z > 0 )
            %Add the syntetic peImage for the particle to the
            %synthetic image of our whole packing 
            x = floor(pres(n).x); %interger rounded x coordinate of the current particle
            y = floor(pres(n).y); %interger rounded y coordinate of the current particle
            r = pres(n).r; %radius (in pixels) of the current particle
            sImg = pres(round(n)).synthImg; %synthetic force image for the current particle
            sx = size(sImg,1)/2; %width of the synthetic force image of the current particle
            sy = size(sImg,2)/2; %heights of the synthetic force image of the current partice
            bigSynthImg(round(y-sy+1):round(y+sy),round(x-sx+1):round(x+sx)) = bigSynthImg(round(y-sy+1):round(y+sy),round(x-sx+1):round(x+sx))+sImg; %Add the syntetic Force Image of the current particle to the appropriate location
        end
    end

imshow(bigSynthImg);
drawnow;
imwrite(bigSynthImg,[directory,'synthImg/',imagefiles(frame).name(1:end-4),'update.jpg'])
%%

% figure;
% histogram([pres.forces])
% if verbose
% figure;
% histogram(rmoutliers(contactforces, "mean" ))
% end
%clear pres
end
end