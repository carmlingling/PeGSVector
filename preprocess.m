%function particle_detect(directory)
% A script to find particle locations
function preprocess(directory, imname, boundaryType, verbose)


images=dir([directory,imname]);
nFrames = length(images);
%nFrames = 15

cen = [71+5313/2, 110+5313/2];
rad = [2783/2, 5313/2];
%% for doing camera transforms for both camera lens distortion and particle distance correction on annulus


if boundaryType == "annulus"
    if not(isfolder(append(directory,'warpedimg'))) %make a new folder with warped images
    mkdir(append(directory,'warpedimg'))
    end
    im1=imread([directory,images(1).name]);
       
        %trim the image to get rid of edge and center points 
        ncols = size(im1,2);
        nrows = size(im1,1);
        for i =1:ncols
          for j=1:nrows
             CC(j,i)=sqrt((j-cen(2))^2+(i-cen(1))^2);
          end
        end
        Q=((CC>=rad(1)) & (CC<=rad(2)));
    parfor (frame = 1:nFrames)
    %frame = frame+130
   
    
    if ~isfile(append(directory,'warpedimg/',images(frame).name(1:end-4),'warped.tif')) %check to see if images are already warped
        frame
        im1=imread([directory,images(frame).name]);
        ncols = size(im1,2);
        nrows = size(im1,1);
        for i =1:nrows
            for j =1:ncols
               if(Q(i,j)~=1)
               im1(i,j,:)=0;
               end
            end
        end
        'trimmed'
        if verbose
          imshow(im1)
          display now;
        end
        im1 = padarray(im1,[400 400], 'both'); % have to pad image with some amount of extra pixels to not lose any edge material when warping
        [nrows, ncols,~] = size(im1);

        [xi,yi] = meshgrid(1:ncols,1:nrows); %transform the points to polar coordinates
        xt = xi - ncols/2;
        yt = yi - nrows/2;
        [theta,r] = cart2pol(xt,yt);

        d = -6.5*r.^2/(200*(925+6.5)); %6.5 is the thickness of the particles in mm, 1000 is distance between particles and camera lens in mm
        rmax = max(r(:));
        s1 = d+r;
        [ut,vt] = pol2cart(theta,s1);
        ui = ut + ncols/2;
        vi = vt + nrows/2;
        %fill = [187, 22, 54];
        ifcn = @(c) [ui(:) vi(:)];
        tform = geometricTransform2d(ifcn);
        'done transform' %finishes creating the geometric transform
    
        [h,w,~] = size(im1);
        [xLimitsOut,yLimitsOut] = outputLimits(tform,[1,w],[1,h]);
        Rout = imref2d(size(im1));
        Rout.XWorldLimits = xLimitsOut;      
        Rout.YWorldLimits = yLimitsOut; 
    
        im =imwarp(im1,tform, "OutputView",Rout);

    
        imwrite(im, append(directory,'warpedimg/',images(frame).name(1:end-4),'warped.tif'))

        'done image warp'
    
%     else
%     im = imread(append(directory,'warpedimg/',images(frame).name(1:end-4),'warped.tif'));

    end
    end
end

    %dlmwrite([directory,images(frame).name(1:end-4),'centers_Improved.txt'],particle)