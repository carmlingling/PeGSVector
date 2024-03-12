%function particle_detect(directory)
% A script to find particle locations
function particleDetect(directory, imname, radiusRange, boundaryType, verbose)


    if boundaryType == "annulus"
        
        images=dir([directory, '/warpedimg/',imname(1:end-4),'warped.tif']);
        nFrames = length(images);
         cen = [71+5313/2, 110+5313/2];
         rad = [2783/2, 5313/2]; %measured in imageJ, pixels
         R_overlap = 15; %how much do the particle have to overlap by to remove double detected particles

  

        for frame = 1:nFrames
            
            im = imread([images(frame).folder,'/', images(frame).name]);
            red = im(:,:,1);
            green = im(:,:,2);
            red = imsubtract(red, green*0.2);
            

            red = imadjust(red, [0.20,0.65]); %this works for annulus, might need to tweak
%             figure;
%             imshow(red)
%             drawnow
            sigma = 50; % choosen by visual inspection
            G = fspecial('gaussian', 3*sigma+1, sigma);
            yb = imfilter(red, G, 'replicate');
            %imshow(yb)
            red = bsxfun(@minus, red,yb);
%             figure;
%             imshow(red)
%             title('red2')
%             drawnow
            green =imadjust(green);
            %
            if verbose == true
                h1 = figure(1);
                ax1 = axes('Parent', h1);

                subplot(1, 2, 1)
                imshow(red)

                subplot(1, 2, 2)
                imshow(green);
                hold on;
                axis on
            end


            [centers,radii,metrics]=imfindcircles(red,radiusRange,'objectpolarity','dark','sensitivity',0.945,'method','twostage','EdgeThreshold',0.02);
%% 
            xt = centers(:,1);
            yt = centers(:,2);
            rt = radii;

            for r=1:length(rt) %binarize radius
                if rt(r)>49
                    rt(r) = 55;
                else
                    rt(r) = 44;
                end
            end


            %convert back to real space
            [midx,~] = size(red);
            [theta,r] = cart2pol(xt-midx/2,yt-midx/2);

            d = -6.5*r.^2/(200*(925+6.5)); %6.5 is the thickness of the particles in mm, 1000 is distance between particles and camera lens in mm
            rmax = max(r(:));
            s1 = d+r;
            [ut,vt] = pol2cart(theta,s1);
            ui = ut + midx/2;
            vi = vt + midx/2;

            ifcn = @(c) [ui(:) vi(:)];
            tform = geometricTransform2d(ifcn);
            [uv] = transformPointsInverse(tform, [0,0]); %particle original coordinates
            u = uv(:,1)-400;
            v = uv(:,2)-400;
%             figure;
%             imold = imread([directory, '/', images(frame).name]);
%             imshow(imold)
%             if verbose
%                 viscircles([u, v], rt)
%                 N = length(uv);
%                 for n=1:N
%             	    text(u(n),v(n),num2str(n),'Color','w');
%                 end
%                 hold on;
%             end

            %remove some of the misfound particles too close to the center
            radialPos = sqrt((u-cen(1)).^2+(v-cen(2)).^2);
            closeind = find(radialPos <= rad(1)+15 );
            closeind = sortrows(closeind, 'descend');


            xt(closeind) = [];
            yt(closeind) = [];
            radii(closeind) = [];
            rt(closeind) = [];
            metrics(closeind)= [];
            u(closeind) = [];
            v(closeind) = [];

            if verbose
                viscircles([u, v], rt,'EdgeColor', 'b')
            end
            %%
            %now we look for particles with a dramatic overlap
            dmat = pdist2([u,v],[u,v]); %Creates a distance matrix for particle center locations
            rmat = rt + rt'; %Makes a combination of radii for each particle

            friendmat = dmat < (rmat - 25) & dmat~=0; %Logical "friend" matrix
            [f1, f2] = find(friendmat == 1);


            badind = zeros(length(f1),1);

            M = length(f1);
            for n=1:M
                if metrics(f1(n)) > metrics(f2(n))
                    badind(n) = f2(n);
                else
                    badind(n) = f1(n);
                end
            end
            badind = badind(badind~=0);
            badind = unique(badind);
            badind = sortrows(badind, 'descend');


            xt(badind) = [];
            yt(badind) = [];
            radii(badind) = [];
            rt(badind) = [];
            metrics(badind)= [];
            u(badind) = [];
            v(badind) = [];

            if verbose
                viscircles([u, v], rt,'EdgeColor','g'); %draw particle outline
                hold on;
            end
            %%
            dmat = pdist2([u,v],[u,v]);
            rmat = rt + rt';
            friendmat = dmat < (rmat -8 ) & dmat~=0; %Logical "friend" matrix

            % %friendmat = triu(friendmat); %Only examine the upper triangle portion (no repeats)
            [f1, f2] = find(friendmat == 1);
            % [f3, f4] = find(friendmat2 == 1);
            badind = unique(f2);
            M = length(badind);
            toobig = zeros(M, 1);
            %
            badin2 = unique(f2);

            for n=1:M
                sum(f2 == badin2(n))
                if  sum(f2 == badin2(n))>2
                    if rt(badind(n)) > 49
                        rt(badind(n)) = 44;
                        toobig(n) = badind(n);

                    end
                end
            end
            toobig = toobig(toobig~=0);

            %viscircles([u(toobig),v(toobig)], rt(toobig), edgecolor = 'y');



            %%
            dtol = 25;
            radialPos = sqrt((u-cen(1)).^2+(v-cen(2)).^2);
            owi= find(radialPos <= rad(2)+2.5*dtol &radialPos >=rad(2)-2.5*dtol);
            iwi = find(radialPos <= rad(1)+3.5*dtol &radialPos >=rad(1)-2.5*dtol);
            edges = zeros(length(u), 1);
            edges(owi) = 1;
            edges(iwi) = -1;
            particle = [xt, yt, rt, edges];
            %viscircles([u(owi), v(owi)], rt(owi), 'Edgecolor', 'c');
            %viscircles([u(iwi), v(iwi)], rt(iwi), 'Edgecolor', 'm');
%             imshow(red)
%             viscircles([xt, yt], rt, 'Edgecolor', 'm');
            dlmwrite([directory,'warpedimg/', images(frame).name(1:end-4),'_centers.txt'], particle)
            %save([directory, images(frame).name(1:end-4),'centers_Improved.mat'], 'particle')
            % save('./DATA/test5centers_Improved.mat', 'particle')
        end

    elseif boundaryType == "airtable"
        dtol = 10;

        images=dir([directory,imname]);
        nFrames = length(images);
        for frame = 1:nFrames
            %for frame = 1:1
            %frame = 663
            im = imread([directory,images(frame).name]);
            red = im(:,:,1);
            green = im(:,:,2);
            red = imsubtract(red, green*0.05);
            %imfindcircles(red,radius,'objectpolarity','dark','sensitivity',0.945,'method','twostage','EdgeThreshold',0.02)
            [centers, radii, metrics] = imfindcircles(red,radiusRange,'ObjectPolarity','bright','Method','TwoStage','Sensitivity',0.945);
%             for r=1:length(radii) %binarize radius
%                 if radii(r)>65
%                     radii(r) = 70;
%                 else
%                     radii(r) = 49;
%                 end
%             end

            
            
            x = centers(:,1);
            y = centers(:,2); %for brevity in coding the next bit
            dmat = pdist2([x,y],[x,y]); %Creates a distance matrix for particle center locations
            rmat = radii + radii'; %Makes a combination of radii for each particle

            friendmat = dmat < (rmat - 25) & dmat~=0; %Logical "friend" matrix
            [f1, f2] = find(friendmat == 1);


            badind = zeros(length(f1),1);

            M = length(f1);
            for n=1:M
                if metrics(f1(n)) > metrics(f2(n))
                    badind(n) = f2(n);
                else
                    badind(n) = f1(n);
                end
            end
            badind = badind(badind~=0);
            badind = unique(badind);
            badind = sortrows(badind, 'descend');


            centers(badind,:) = [];
            
            radii(badind) = [];

            imshow(red);
            viscircles(centers, radii);
            if length(radii)>29
                for m = 1:length(radii)-29
                 text(centers(m+29,1), centers(m+29,2), num2str(m+29) )
                end
            end
            drawnow


            if length(radii)>29
                for m = length(radii)-29:-1:1
                    centers(m+29,:)=[];
                    radii(m+29)=[];
                end
            end
            
            lpos = min(centers(:,1)-radii);
            rpos = max(centers(:,1)+radii);
            upos = max(centers(:,2)+radii);
            bpos = min(centers(:,2)-radii);
            lwi = find(centers(:,1)-radii <= lpos+dtol);
            rwi = find(centers(:,1)+radii >= rpos-dtol);
            uwi = find(centers(:,2)+radii >= upos-dtol);
            bwi = find(centers(:,2)-radii <= bpos+dtol); %need to add edge case of corner particle
            viscircles(centers(lwi,:), radii(lwi), 'EdgeColor', 'b');
            viscircles(centers(rwi,:), radii(rwi),'EdgeColor', 'k');
            viscircles(centers(bwi,:), radii(bwi),'EdgeColor', 'w');
            viscircles(centers(uwi,:), radii(uwi),'EdgeColor', 'g');
            edges = zeros(length(radii), 1);
            edges(rwi) = 1;
            edges(lwi) = -1;
            edges(uwi) = 2;
            edges(bwi) = -2;
            for q = 1:length(uwi)
                d = centers(uwi(q),2)+radii(uwi(q))-799;
                if d>=0
                    
                    centers(uwi(q),2)=centers(uwi(q),2)-d-0.5;
                end
            end
            particle = [centers(:,1), centers(:,2), radii, edges];
            [directory, images(frame).name(1:end-4),'_centers.txt']
            dlmwrite([directory, images(frame).name(1:end-4),'_centers.txt'], particle);

        end
    end
    %dlmwrite([directory,images(frame).name(1:end-4),'centers_Improved.txt'],particle)