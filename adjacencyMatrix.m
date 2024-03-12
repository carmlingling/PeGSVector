function adjacencyMatrix(directory, fileNames, boundaryType, frameidind,verbose)
%  directory = '/eno/cllee3/DATA/esbilili/biaxial/img2/' % location of image files
% % % %topDirectory = '/eno/cllee3/DATA/230502_2/'
% % % %topDirectory = './DATA/test/Step09/'
% fileNames = 'img*' %image format and regex
% %
if boundaryType == "annulus"
    directory = [directory, 'warpedimg/']
end
files = dir([directory,fileNames(1:end-4),'_solved_update.mat']) %which files are we processing ?
nFrames = length(files) %how many files are we processing ?
%nFrames = 1
%camImageFileName = '/eno/cllee3/DATA/230420/150Hz_2fps_run1/warpedimg/150Hz_*warped.tif'
%camImageFileName = [directory, '100Hz_1fps_Img0021warped.tif']
%camImageFileName = '/eno/cllee3/DATA/230428/run1/warpedimg/150Hz_*warped.tif'
%PARAMETERS NEEDED TO RUN THIS SCRIPT ARE SET HERE
go = true;
fmin = 0.000001; %minimum force (in Newton) to consider a contact a valid contact
fmax = 1000; %maximum force (in Newton) to consider a contact a valid contact
emax = 2800; %maximum fit error/residual to consider a contact a valid contact
fs=16; %plot font size
verbose = false; %make lots of plots as we go

% jut in case your data coordinate system is offset from the image
% coordinate system (i.e. you only processed a small part of a larger
% image)
% xoffset = 1000;
% yoffset = 700;
% xsize = 2470-xoffset;
% ysize = 2260-yoffset;

%%
%Global Metrics will be stored in these structures
%allContacts = struct('fAbs',0,'fNorm',0,'fTan',0); %data structure to store information about contacts
%aID = 1; %Global contact counter over all contacts in all cycles.

if go==true
for cycle = 1245:nFrames %loop over these cycles 
    %cycle = cycle+820
    clearvars particle;
    clearvars contact;
    
    %input filnames
    peOutfilename = [directory,files(cycle).name] %input filename 
    if boundaryType =="airtable"
        camImageFileName = [peOutfilename(1:end-18),'.jpg'];  %adjusted force image filename
    elseif boundaryType == "annulus"
        camImageFileName = [peOutfilename(1:end-18),'.tif'];  %adjusted force image filename
    end

    %camImageFileName = [directory,'100Hz_1fps_Img0021warped.tif'];  %adjusted force image filename 
   
    %output filenames
    %workspacefilename =  [peOutfilename(1:end-9),'-postProcessingWorkspace.mat']; 
    workspacefilename =  [peOutfilename(1:end-10),'-postProcessingWorkspace.mat']; 
    %synthImgFilename = [peOutfilename(1:end-10),'update-Synth.jpg'];  %output filename 
    adjMatAbsFilename = [peOutfilename(1:end-11),'-joAdjacencyAbs.dlm'];  %output filename 
    adjMatNorFilename = [peOutfilename(1:end-11),'-joAdjacency_N.dlm'];
     adjMatTanFilename = [peOutfilename(1:end-11),'-joAdjacency_T.dlm'] ;
     particlelocFilename = [peOutfilename(1:end-10),'-particle.dlm'];
    % NO PARAMETERS SHOULD BE SET BY HAND BELOW THIS LINE

    %check if the data we want to read exists
    %if it does, load it, else abort
    if ~(exist(peOutfilename, 'file') == 2) %if the file we try to open does not exist
        display(['File not Found:', peOutfilename]); %complain about it
        return %and end the execution of this script
    else
        load(peOutfilename); %read peDiscSolve ouput
        particle = pres;
        NN = length(particle);
        IDN = max([particle.id]);
    end

    % Gernate Combined Syntehtic Force Image
    %img = imread(camImageFileName); %camera force image
    %img = imcrop(img,[xoffset, yoffset, xsize, ysize]); %particle image
    
    
    %DATA EVALUATION AND ANALYSIS STARTS HERE

    %particle(1:size(data,1)) = struct('id',0,'x',0,'y',0,'z',0,'fx',0,'fy',0); %data structure to store particle information
    contact = struct('id1',0,'id2',0,'x',0,'y',0,'fAbs',0,'fNorm',0,'fTan',0,'alpha',0,'beta',0,'contactX',0,'contactY',0,'error',0); %data structure to store information about contacts
    cID = 1; %contact counter
    %A = zeros(NN); %empty binary adjacency matrix
    W = NaN(IDN); %empty force weighted adjacency matrix
    N = NaN(IDN); %empty normal force weighted adjacency matrix
    T = NaN(IDN); %empty tangential force weighted adjacency matrix
    %P = NaN(NN);
    for n = 1:NN %for each particle
        err = particle(n).fitError; %get fit error 
        z = particle(n).z; %get coordination number
        r = particle(n).r; %get particle radius in pixel
        rm = particle(n).rm; %get particle radius in meters
        if(length(particle(n).neighbours) > 0) % particle is in contact
            contacts = particle(n).neighbours; %get IDs of all contacting particles
            betas = particle(n).betas+pi; %get the beta angle (position of contact point) associated with each contact
            forces = particle(n).forces; %get the force associated with each contact
            alphas = particle(n).alphas; %get the alpha angle (direction of force) associated with each contact
            
            for m=1:length(forces) %for each contact

                %if(forces(m) > fmin && err < emax && forces(m) < fmax) %is this a valid contact ?
                if (forces(m) < fmax && forces(m)> 0)% && err <emax)
                    %put information about the first particle involved in this
                    %contact in the corresponding particle structure vector

                    %ideally the accumulated fx and fy should be zero, that is
                    %the particle is in force balance
                    %particle(n).fx = particle(n).fx + forces(m) * cos(betas(m)-pi); %x component of total force vector %CHECK AGAIN IF THIS IS GEOMETRICALLY CORRECT
                    %particle(n).fy = particle(n).fy + forces(m) * sin(betas(m)-pi); %y component of total force vector %CHECK AGAIN IF THIS IS GEOMETRICALLY CORRECT
                    %particle(n).z = particle(n).z+1; %increment the real contact number for the current particle

                    %put all the information about this contact
                    %into the contact struct vector
%                     contact(cID).id1 = particle(n).id; %first particle involved in this contact
                    targetid = contacts(m);
                    ids = [particle.id];
                    tind1m = ids == targetid;
                    tind1 = find(tind1m);
%                     contact(cID).id2 = targetid; %second particle involved in this contact 
%                     contact(cID).x = particle(n).x;
%                     contact(cID).y = particle(n).y;
%                     contact(cID).fAbs = forces(m); %absolute force
%                     contact(cID).fNorm = forces(m)*cos(alphas(m)); %normal force (see Eq. 4.16)
%                     contact(cID).fTan = forces(m)*sin(alphas(m)); %tangential force
%                     contact(cID).alpha = alphas(m); %the alpha angle (direction of force) associated with this contact
%                     contact(cID).beta = betas(m); %the beta angle (position of contact point) associated with this contact  
%                     contact(cID).contactX =  r * cos(betas(m)-pi); %x component of vector to contact point
%                     contact(cID).contactY =  r * sin(betas(m)-pi); %y component of vector to contact point
%                     contact(cID).error = err; %fit error for this particle (the first particle in the contact)
% % 
%                     cID = cID + 1; %increment contact counter
%                     
%                     allContacts(aID).fAbs = forces(m); %absolute force
%                     allContacts(aID).fNorm = forces(m)*cos(alphas(m)); %normal force (see Eq. 4.16)
%                     allContacts(aID).fTan = forces(m)*sin(alphas(m)); %tangential force
%                     
%                     aID = aID+1;

                    %build some adjacency matrices
                    %if (contacts(m)>0) %correct for negative contact IDs in peDiscsolve, i.e.non-wall contacts only
                         %A(n,contacts(m)) = 1; %mark contact in the binary adjacency matrix
                         W(particle(n).id,particle(tind1).id) = real(forces(m)); %write the corrsponding force as a weight into an adjacency matrix
                         N(particle(n).id,particle(tind1).id) = real(forces(m))*cos(alphas(m)); %write the corrsponding normal force as a weight into an adjacency matrix
                         T(particle(n).id,particle(tind1).id) = real(forces(m))*sin(alphas(m)); %write the corrsponding tangential force as a weight into an adjacency matrix
                         %P(n, tind1) = [x, y, rm]
                    %end
                end
            end
        end
       
    end
    %%
% length(nonzeros(W))
% edges = 10.^(-5:0.1:2);
% [N,edges] = histcounts(W,edges,'Normalization','countdensity');
% figure;
% g = histogram('BinEdges',edges,'BinCounts',N);
% set(gca, "Xscale", "log")
%     figure;
%     plot(nonzeros(W), '.')
%      drawnow;
    %make sure our adjacency matrix is nice
    %discard singular contacts
%         B = double((W.*W') > 0);
%         W = (W.*B);
% 
%         B = double((N.*N') > 0);
%         N = (N.*B);
% 
%         B = double((T.*T') > 0);
%         T = (T.*B);

    %average reciprocal contacts so the matrix is fully symmetric
%         W = (W+W')/2;
%         N = (N+N')/2;
%         T = (T+T')/2;
    %DATA OUTPUT STARTS HERE

    
    %write out the weighte adjacency matrix for the current frame
    dlmwrite(adjMatAbsFilename,W); 
    dlmwrite(adjMatNorFilename,N);
    dlmwrite(adjMatTanFilename,T);
    %write out the synthetic force image for the current frame
    %imwrite(bigSynthImg,synthImgFilename);
   

    %Plot what we have learned so far
    
        %set up a new figure environment and scale it appropriately
        %close all;
        %pFig = figure ;
        %set(pFig,'Position',[0 0 2000 1200]);
        %set(pFig,'paperorientation','landscape');

        %this suplot shows the original image, overlayed with the
        %detected contacts

        if verbose
        %figure(1)
            %read and display the original image used as input to peDisc
            %img = imcrop(imread(camImageFileName),[xoffset, yoffset, xsize, ysize]); %force image
            img = imread(camImageFileName); %force image
            imshow(img); hold on;
            colormap(gray);
            %plot the centers of all particles associated with a contact
            plot([contact.x],[contact.y],'or')
            %plot arrows from the centers of all particles associated with a contact to
            %the contact point
            f = [contact.fAbs];
            norm = max(max(f));
            shift = min(min(f));
            linewidths = 10*(f-shift+0.001)/norm;
            for m= 1:length(f)
                quiver([contact(m).x],[contact(m).y],[contact(m).contactX],[contact(m).contactY],0,'LineWidth',linewidths(m), Color='b')
            end
                
            %set font sizes and labels
            
            set(gca,'FontSize',fs);
            title('camera image','FontSize',fs);
            drawnow;
%          figure(2)
%             %read and display the original image used as input to peDisc
%             %img = imcrop(imread(camImageFileName),[xoffset, yoffset, xsize, ysize]); %force image
%             img = bigSynthImg; %force image
%             imshow(img); hold on;
%             colormap(gray);
%             %plot the centers of all particles associated with a contact
%             %plot([contact.x],[contact.y],'or')
%             %plot arrows from the centers of all particles associated with a contact to
%             %the contact point
%             %quiver([contact.x],[contact.y],[contact.contactX],[contact.contactY],0,'LineWidth',1.5)
%             %set font sizes and labels
%             set(gca,'FontSize',fs);
%             title('camera image','FontSize',fs);
%             hold off;
%          figure(3)
%             imagesc(W); 
%             colormap(jet);
        end
    
 
end
end
%%
[directory,fileNames(1:end-4),'_solved-joAdjacencyAbs.dlm']
Adj_ABS = dir([directory,fileNames(1:end-4),'_solved-joAdjacencyAbs.dlm']);
Adj_T = dir([directory,fileNames(1:end-4),'_solved-joAdjacency_T.dlm']);
Adj_N = dir([directory,fileNames(1:end-4),'_solved-joAdjacency_N.dlm']);
if ~isfile([directory, 'centers_tracked.txt']);
    centers = dlmread([directory, 'centers_tracked.txt'],',', 1,0);
    %Adj_T_Newton = dir([directory,'*newtonized2-joAdjacency_NewtonizedT.dlm']);
    %Adj_N_Newton = dir([directory,'*newtonized2-joAdjacency_NewtonizedN.dlm']);
    % Adj_theta = dir([directory,'*solved-joAdjacencyN.dlm']);
    % add as many as you want...

    %% Make a list
    Adj_list = [];
    %numbers = [10, 18]

    unique(centers(:,1));
    for i = 1:nFrames
        %for i = 1:length(2)
        frameid = str2num(Adj_T(i).name(frameidind:frameidind+3))
        T = dlmread([directory,Adj_T(i).name]);
        
        N = dlmread([directory,Adj_N(i).name]);
        %T_Newtonized = dlmread([directory,Adj_T_Newton(i).name]);
        %N_Newtonized = dlmread([directory,Adj_N_Newton(i).name]);

        %     Theta = dlmread([directory,Adj_theta(i).name]);

        %list = [];
        con = find(centers(:,1) == i); %the data that goes with the frame
        con;
        pos2ind = centers(con, 2); %
        d = ~isnan(T); %where there is data that is actually a force
        [row , col] = find(d==1); %which indices correspond to that contact
        pos2ind;
        newrow = pos2ind(row);
        newcol = pos2ind(col);
        ind=sub2ind(size(T),row,col);
        ind;
	    length(row), length(newrow), length(newcol), real(T(ind)), real(N(ind));
        %     list_T = [row , col , T(row,col) , N(row,col) , Theta(row,col)];
        list = [ones(length(row),1).*frameid,newrow , newcol , T(ind) , N(ind)];%, T_Newtonized(ind),N_Newtonized(ind)];


        Adj_list=[Adj_list;list];
    end
else
    Adj_list = [];

    
    for i = 1:length(Adj_T)
    Adj_T(i).name;
    frameid = str2num(Adj_T(i).name(frameidind:frameidind+3))
    T = dlmread([directory,Adj_T(i).name]);
    N = dlmread([directory,Adj_N(i).name]);
    

%     Theta = dlmread([directory,Adj_theta(i).name]);

    list = [];
    d = ~isnan(T);
    [row , col] = find(d==1);
    ind=sub2ind(size(T),row,col);
    
%     list_T = [row , col , T(row,col) , N(row,col) , Theta(row,col)];
    list = [ones(length(row),1).*frameid,row , col , T(ind) , N(ind)];


    Adj_list=[Adj_list;list];

end

    


    end

%fram number, particle 1, particle id 2, tangential force, normal force
dlmwrite([directory,'Adjacency_list.txt'],Adj_list);




