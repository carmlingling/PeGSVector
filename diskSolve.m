%function particleSolve(directory, fileNames, boundaryType, verbose)

directory = '/eno/cllee3/DATA/240506/test/';
%topDirectory = '/Users/carmenlee/Desktop/20150731reprocesseduniaxial/'
% %topDirectory = './DATA/test/Step09/'
fileNames = '200Hz*.tif'; %image format and regex
frameidind = 16;
%

boundaryType = "annulus"; %if airtable use "airtable" if annulus use "annulus"
radiusRange = [40, 57];
directoryini = directory
if boundaryType == "annulus"
    directory = [directory, 'warpedimg/']
    images = dir([directory, fileNames(1:end-4), '.tif'])
else
    images = dir([directoryini, fileNames])
end
[directoryini,'/particles' fileNames(1:end-4),'_preprocessings.mat']
files = dir([directoryini, 'particles/', fileNames(1:end-4),'_preprocessings.mat'])

if not(isfolder(append(directoryini,'synthImg'))) %make a new folder with warped images
    mkdir(append(directoryini,'synthImg'));
end

if not(isfolder(append(directoryini,'solved'))) %make a new folder with warped images
    mkdir(append(directoryini,'solved'));
end
%how much of the particle diameter is used to fit the synthetic image 
%(1 = use everything). Change this parameter only if the fit doesn't work 
%correctly because of imaging artifacts at the particle boundary.
maskradius = 0.96;% 
scaling = 0.5; %scale the image by this factor before doing the fit

%do we want to see each particle on screen while it is fitted ?
%verbose = false; 

%fit options: play around with these, depending on the quality of your
%images and the accuracy you want your result to have. 
%Setting a good TolFun value can considerabely speed up your processing.



%% 

%fo_c = parallel.pool.Constant(@() fitoptions('lsqnonlin','Algorithm','levenberg-marquardt','MaxIter',200,'MaxFunEvals',400,'TolFun',0.01,'Display','final-detailed'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%There should be no need for user input below this line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
nFrames = length(files) %how many files are we processing ?
%nFrames = 145
h3 = figure();
hAx1 = subplot(1,1,1,'Parent', h3);

for frame = 1:nFrames %loop over these frames 
%for frame =1:9
    frame = frame
    fileName = [directoryini,'particles/',files(frame).name]
    data = load(fileName);
    particle = data.particle;
    particle = PeGSDiskSolve(particle, maskradius, scaling, directoryini, files(frame).name);

    img = imread([directory, images(frame).name(1:end-4), '.tif']);
    bigSynthImg = zeros(size(img,1),size(img,2)); %make an empty image with the same size as the camera image
    NN = length(particle);
    
    for n=1:NN %for all particles
        %display(['fitting force(s) to particle ',num2str(n)]); %status indicator
        if (particle(n).z > 0 )
            %Add the syntetic peImage for the particle to the
            %synthetic image of our whole packing 
            x = floor(particle(n).x); %interger rounded x coordinate of the current particle
            y = floor(particle(n).y); %interger rounded y coordinate of the current particle
            r = particle(n).r; %radius (in pixels) of the current particle
            sImg = particle(round(n)).synthImg; %synthetic force image for the current particle
            sx = size(sImg,1)/2; %width of the synthetic force image of the current particle
            sy = size(sImg,2)/2; %heights of the synthetic force image of the current partice
            bigSynthImg(round(y-sy+1):round(y+sy),round(x-sx+1):round(x+sx)) = bigSynthImg(round(y-sy+1):round(y+sy),round(x-sx+1):round(x+sx))+sImg; %Add the syntetic Force Image of the current particle to the appropriate location
            
        end
    
    end
   
    imshow(bigSynthImg, 'Parent', hAx1);
%     hold (hAx1, 'on');
%     for n=1:NN
%         viscircles([particle(n).x, particle(n).y], particle(n).r)
%         text(particle(n).x, particle(n).y, num2str(particle(n).fitError))
%     end     
drawnow;
    imwrite(bigSynthImg, [directoryini,'synthImg/', images(frame).name(1:end-4), '-Synth.jpg']);
end
