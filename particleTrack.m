%this function links particles from their position file. exports a
%centers_tracked.txt file that contains frame, particle id, x, y, r and
%edge status. If it is an annulus, it also exports the unwarped positions

function particleTrack(directory,imname, boundaryType, frameidind, verbose)

%handling specific file structure stuff
if boundaryType == "annulus"
    directory = [directory, 'warpedimg/'];
    datafiles = dir([directory, imname(1:end-4),'warped_centers.txt']);
else
    %datafiles = dir([directory, imname(1:end-4),'_centers.txt'])
    datafiles = dir([directory, imname]);
end

dtol=10;

%find guess for array size
posData = load([directory, datafiles(1).name]);
nFrames = length(datafiles);    
skipamount = length(posData)+200;
centers = nan(nFrames*skipamount, 6);
if boundaryType == "airtable"
    for frame = 1:nFrames
        frame
        %posData = dlmread([directory, datafiles(n).name]);
    
        posData = load([directory, datafiles(frame).name]);
    
        frameid = str2double(datafiles(frame).name(frameidind:frameidind+3))
        %posData = new;
        frameNumber = ones(length(posData),1)*frameid;
        id = 1:length(posData);
        
        %posArray = [posArray; [posData(:,1), posData(:,2),  frameNumber]];
        %radArray = [radArray; posData(:,3)];
        %edgesArray =[edgesArray;posData(:,4)];
        
        x = extractfield(posData, 'x')'; %for using Ephraim's already detected particles
        y = extractfield(posData, 'y')';
        r = (extractfield(posData, 'r'))';
        lpos = min(x-r);
        rpos = max(x+r);
        upos = max(y+r);
        bpos = min(y-r);
        edges = zeros(length(r), 1);
        for k= 1:length(x)
            if x(k)-r(k) <=lpos+dtol
                edges(k) = -1;
            elseif x(k)+r(k) >=rpos-dtol
                edges(k) = 1;
            elseif y(k)+r(k) >= upos-dtol
                edges(k) = 2;
            elseif y(k)-r(k) <= bpos+dtol
                edges(k) = -2;
            end
        end
        centers((frame-1)*skipamount+1:(frame-1)*skipamount+length(x)+1) = [frameNumber, id', x, y, round(r), edges];
    end
else
    for frame = 1:nFrames
        frame
        posData = readmatrix([directory, datafiles(frame).name]);
    
        %posData = load([directory, datafiles(n).name]);
        
        frameid = str2double(datafiles(frame).name(frameidind:frameidind+3));
        
        frameNumber = ones(length(posData),1)*frameid;
        id = 1:length(posData);
        %posArray = [posArray; [posData(:,1), posData(:,2),  frameNumber]];
        
        centers((frame-1)*skipamount+1:(frame-1)*skipamount +length(posData),:) = [frameNumber, id', posData(:,1), posData(:,2), posData(:,3), posData(:,4)];
    end
end


%tidying up the centers dataframe
centers(any(isnan(centers),2),:)=[];
posArray=[centers(:,3), centers(:,4), centers(:,1)];

%%uncomment this if you don't care about particle order
% out = fopen([directory,'centers_tracked.txt'],'w');
% fprintf(out,['frame', ',', 'particleID', ',', 'x' ',' 'y', ',','r',',','edge''\n']);
% fclose(out);
% dlmwrite([directory,'centers_tracked.txt'], centers, 'delimiter',',','-append');
%   exit



%linking the particles, where we remember the particles over 50 frames, iin
%2 dimensions, limit is 30 pixels of motion per frame
param = struct('mem', 50, 'good', 0,'dim',2,'quiet',0);
linked = track(posArray, 30, param);



%carry radii information along too
radind = zeros(length(centers),1);
edgesind = zeros(length(centers),1);
for l = 1:length(centers)
    ind = find(posArray(l,1) == linked(:,1) & posArray(l,2) ==linked(:,2));
    radind(ind) = centers(l,5);
    edgesind(ind) = centers(l,6);
end

rearrange = [linked(:,3),linked(:,4), linked(:,1), linked(:,2), radind, edgesind];
rearrange = sortrows(rearrange);
% %% don't know why I wrote this but I am scared to delete
% P = length(rearrange);
% rearrange = [rearrange; zeros(1700, 6)];
% all_particleIDs = unique(rearrange(:,2));
% for l=1:length(all_particleIDs)
%     framevalue = rearrange(find(rearrange(:,2)==all_particleIDs(l)),1);
%     frameind = find(rearrange(:,2)==all_particleIDs(l));
%     for m=1:length(framevalue)-1
%         if framevalue(m+1)-framevalue(m)>1
%             counter = 1;
%             filler = framevalue(m)+1:framevalue(m+1)-1;
%             x = rearrange(frameind(m+1),3);
%             x2 = rearrange(frameind(m), 3);
%             y = rearrange(frameind(m+1), 4);
%             y2 = rearrange(frameind(m), 4);
%             for n=1:length(filler)
%                 xi = x+uint16(n*(x2-x)/length(filler));
%                 yi = y+uint16(n*(y2-y)/length(filler));
%                
%                 rearrange(P+counter,:) = [filler(n), all_particleIDs(l), xi, yi, rearrange(frameind(m),5), 0];
%                 counter = counter+1;
%             end
%        
%         end
%     end
% end
% rearrange(~any(rearrange,2),:)=[];
% rearrange = sortrows(rearrange);
% length(rearrange)
% 
% 
% 
% 
out = fopen([directory,'centers_tracked.txt'],'w');
fprintf(out,['frame', ',', 'particleID', ',', 'x' ',' 'y', ',','r',',','edge''\n']);
fclose(out);
writematrix(rearrange,[directory,'centers_tracked.txt'], 'WriteMode', 'append');


%exports original positions
if boundaryType == "annulus"
    rearrange2 = rearrange;
    x = rearrange(:,3);
    y = rearrange(:,4);

    midx = 6304; %size of image
    [theta,r] = cart2pol(x-midx/2,y-midx/2);

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

    rearrange2(:,3) = u;
    rearrange2(:,4) = v;
    out = fopen([directory,'centers_tracked_original.txt'],'w');
    fprintf(out,['frame', ',', 'particleID', ',', 'x' ',' 'y', ',','r',',','edge''\n']);
    fclose(out);
    writematrix(rearrange2,[directory,'centers_tracked_original.txt'], 'Writemode', 'append');
% %% 
end
% 

%if you want to see the particle traces
if verbose
    figure;
    image = dir([directory, imname(1:end-4),'warped.tif']);
    pic = imread([directory, image(1).name]);
    imshow(pic);
    hold on;
    N = unique(rearrange(:,2));
    cm = colormap(parula(size(N,1))); 
    for frame = 1:length(N)
        ind = find(rearrange(:,2) ==frame);
        plot(rearrange(ind,3), rearrange(ind,4),'Color',cm(frame,:));
        hold on;
    end
end
