% Jonathan's Photoelastic Disk Solver
%
% Takes input from preprocessing script joDiskPrep.m as of 2016/09/27
% inspired from peDiskSolve by James Puckett (Phd-Thesis 2012)
% http://nile.physics.ncsu.edu
%
% If you use this sovler please cite the follwoing paper
% K.E. Daniels, J. E. Kollmer & J. G. Puckett, "Photoelastic force measurements in granular materials", Rev. Sci. Inst. (201X)
% DOI: XXXXXX
%
% The generation of the synthetic force images is parallelized 
% in a way that each row of the ouput image can be its own worker/thread
% it is usually limited by your number of processing cores though. 
%
% Running this script may take a long time (up to one minute per particle
% on a 2016 desktop computer) so it makes sense to send this off to a high
% performance computer. 
%
% Depending on your experimental data you can play around with the fit
% options, which might speed up processing considerabely.
%
% Last edit on 2016/09/28 by Jonathan Kollmer (jekollme@ncsu.edu)

%close all %housekeeping
%clear all %housekeeping

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%User defined values:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%which files are we processing ?
%directory = '/eno/cllee3/DATA/230502_2/warpedimg/';
%directory = './DATA/test/'
%files = dir([directory, '*preprocessing.mat']); 
%files = dir([directory, '100Hz_1fps_Img002*preprocessing.mat']); 

%how much of the particle diameter is used to fit the synthetic image 
%(1 = use everything). Change this parameter only if the fit doesn't work 
%correctly because of imaging artifacts at the particle boundary.
%maskradius = 0.96;% 
%scaling = 1; %scale the image by this factor before doing the fit

%do we want to see each particle on screen while it is fitted ?
%verbose = false; 

%fit options: play around with these, depending on the quality of your
%images and the accuracy you want your result to have. 
%Setting a good TolFun value can considerabely speed up your processing.

%fo_c = parallel.pool.Constant(@() fitoptions('lsqnonlin','Algorithm','levenberg-marquardt','MaxIter',200,'MaxFunEvals',400,'TolFun',0.01,'Display','final-detailed'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%There should be no need for user input below this line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%nFrames = length(files) %how many files are we processing ?
%nFrames = 3


%for frame = 1:nFrames %loop over these frames 
    %frame = 2
    %frame = frame +3

function particle = disksolve(particle, maskradius, scaling, directory, fileName);
    %fileName = [directory,files(frame).name]; %which file/frame are we processing now ?
    maskradius = maskradius / 2; %I did an unwise choice in naming this
    %load(fileName); %load the particle data file to process
    pres=particle; %initialize a copy of the particle vector to mess around with
    N = length(particle); %number of particles in this frame
    verbose = false;
    display(['processing file ',fileName, 'containing ' ,num2str(N), 'particles']); %status indicator
    %figure;
    parfor (n = 1:N)
        display(['fitting force(s) to particle ',num2str(n)]); %status indicator
        if (particle(n).z > 0 )
                    %bookkeeping
                    fsigma = particle(n).fsigma;
                    %fsigma = 80;
                    rm = particle(n).rm;

%                     %This is the Camera Image
%   %template = im2double(particle(n).forceImage);
                    template= particle(n).forceImage;
                     %template = imadjust(particle(n).forceImage);
                     template = imresize(template,scaling);
%                     template = template-0.1;
%                     %template = (template -0.2); %fine tunes the image, this should happen in preprocessing!
%                     template = template.*(template > 0); %fine tunes the image, this should happen in preprocessing!
%                     template = template*3; %fine tunes the image, this should happen in preprocessing!

                    %size of the force image
                    px = size(template,1); 
                    
                    %plot the experimental image that is going to be fitted
                    %onto
                    if verbose
                        subplot(1,2,1)
                        imshow(template);
                    end
                    
                    %Create initial guesses for each contact force, based
                    %on the gradient squared in the contact region (also
                    %better guess lower than too high)
                    z = particle(n).z;
                    forces = zeros(z,1);
                    cg2s = sum(particle(n).contactG2s);

                    for i=1:z
                            forces(i)=particle(n).f*particle(n).contactG2s(i)/cg2s;%initial guess for force
                    end

                    %Initial guesses for the angles of attack of the forces
                    %for each contact
                    alphas = zeros(z,1);
                    beta = particle(n).betas;

                    %apply force balance to the initial guesses
                    %[alphas,forces] = forceBalance(forces,alphas,beta);
                    %imshow(joForceImg(z, forces,alphas, beta(1:z), fsigma, rm, px, verbose))
                    %create a circular mask
                    cx=px/2;cy=px/2;ix=px;iy=px;r=maskradius*px;
                    [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
                    c_mask=((x.^2+y.^2)<=r^2);   

                    %This is the function we want to fit, i.e. a synthetic
                    %version of the force image with free parameters f (here called par(1:z)) and
                    %alpha (here called par(z+1:z*z) since we have to stuff
                    %them into one vector.
                    func = @(par) joForceImg(z, par(1:z),par(z+1:z+z), beta(1:z), fsigma, rm, px, verbose); %+par(2*z+1); %this is the function I want to fit (i.e. synthetic stres image), the fitting paramters are in vector par
                    %This is the error function we are actually fitting,
                    %that is, the distance between our fit function and the
                    %real particle image in terms of the sum of squares of the pixelwise differnce.
                    %Also a mask is applied to crop
                    %out only the circular part of the particle. 
                    
                    err = @(par) abs(sum(sum( ( c_mask.*(template-func(par)).^2) ))); %BUG: for some reason I sometimes get imaginary results, this should not happen

                    %Set up initial guesses
                    p0 = zeros(2*z, 1);
                    p0(1:z) = forces;
                    p0(z+1:2*z) = alphas;

                    %Do the fit, will also work with other solvers
                    %TODO: make a user defined option to select between
                    %different solvers
                    fitoptions = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxIter',200,'MaxFunEvals',400,'TolFun',0.01,'Display','final-detailed');
                    
                    p=lsqnonlin(err,p0,[],[],fitoptions);

                    %get back the result from fitting
                    forces = p(1:z);
                    alphas = p(z+1:z+z);

                    %resudual
                    fitError = err(p);
                    
                    %generate an image with the fitted parameters
                    imgFit = joForceImg(z, forces, alphas, beta, fsigma, rm, px*(1/scaling), verbose);
                    
                    %since the image gnerator also forces force balance we
                    %have to explicitly do it once more to the values we
                    %are gonna save 
                    [alphas,forces] = forceBalance(forces,alphas,beta);
                    
                    %store the new information into a new particle vector 
                    pres(n).fitError = fitError;
                    pres(n).forces = forces;
                    pres(n).alphas = alphas;
                    pres(n).synthImg = imgFit;
        end
    end
 
    %save the result
    save([directory,'solved/',fileName(1:end-17),'solved.mat'],'pres');
    particle = pres;
    %clearvars('-except', initialVars{:})


%     imread()
%     bigSynthImg = zeros(size(img,1),size(img,2)); %make an empty image with the same size as the camera image
%     for n=1:NN %for all particles
%         %display(['fitting force(s) to particle ',num2str(n)]); %status indicator
%         if (particle(n).z > 0 )
%             %Add the syntetic peImage for the particle to the
%             %synthetic image of our whole packing 
%             x = floor(pres(n).x)-300; %interger rounded x coordinate of the current particle
%             y = floor(pres(n).y)-300; %interger rounded y coordinate of the current particle
%             r = pres(n).r; %radius (in pixels) of the current particle
%             sImg = pres(round(n)).synthImg; %synthetic force image for the current particle
%             sx = size(sImg,1)/2; %width of the synthetic force image of the current particle
%             sy = size(sImg,2)/2; %heights of the synthetic force image of the current partice
%             bigSynthImg(round(y-sy+1):round(y+sy),round(x-sx+1):round(x+sx)) = bigSynthImg(round(y-sy+1):round(y+sy),round(x-sx+1):round(x+sx))+sImg; %Add the syntetic Force Image of the current particle to the appropriate location
%         end
%     end
%     imwrite(bigSynthImg, [fileName(1:end-17),'_Synth-Img.jpg'])
end
