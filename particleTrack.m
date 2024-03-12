% %try to link particle positions over sample dataset
% clear all;
% %directory='/eno/cllee3/DATA/230420/150Hz_2fps_run1/warpedimg/';
% %directory = './DATA/warpedimg/'
% directory = '/eno/cllee3/DATA/230502_2/warpedimg/'
% %directory = '~/Desktop/230502/warpedimg/'
% %directory = './Calibration/';
% %directory = '/eno/cllee3/DATA/Calibration/'
% %directory = '/mnt/ncsudrive/c/cllee3/Documents/'
% %directory = '~/Desktop/230410/'
% imname = '*centers.txt';
% %imname = 'background000.tif'
% %mname = '*jpg'
% %imname = '*tif'
% datafiles=dir([directory,imname]);
% nFrames = length(datafiles);
function trackParticles(directory,imname, boundaryType, frameidind, verbose)
if boundaryType == "annulus"
    directory = [directory, 'warpedimg/']
    [directory,imname(1:end-4),'warped_centers.txt'];
    datafiles = dir([directory, imname(1:end-4),'warped_centers.txt'])
else
    %datafiles = dir([directory, imname(1:end-4),'_centers.txt'])

     datafiles = dir([directory, imname])
end
dtol=10
nFrames = length(datafiles)
posArray = [];
radArray = [];
edgesArray = [];
centers = [];
if boundaryType == "airtable"
for n = 1:nFrames
    n
    %posData = dlmread([directory, datafiles(n).name]);
    
    posData = load([directory, datafiles(n).name]);
    
    frameid = str2num(datafiles(n).name(frameidind:frameidind+3))
    %posData = new;
    frameNumber = ones(length(posData),1)*frameid;
    id = [1:length(posData)];
    %posArray = [posArray; [posData(:,1), posData(:,2),  frameNumber]];
    %radArray = [radArray; posData(:,3)];
    %edgesArray =[edgesArray;posData(:,4)];
    x = extractfield(posData, 'x')';
    y = extractfield(posData, 'y')';
    r = (extractfield(posData, 'r'))';
    lpos = min(x-r)
    rpos = max(x+r)
    upos = max(y+r)
    bpos = min(y-r)
    edges = zeros(length(r), 1);
    for k= 1:length(x)
            if x(k)-r(k) <=lpos+dtol;
                edges(k) = -1;
            elseif x(k)+r(k) >=rpos-dtol
                edges(k) = 1;
            elseif y(k)+r(k) >= upos-dtol
                edges(k) = 2;
            elseif y(k)-r(k) <= bpos+dtol
                edges(k) = -2;
            end
    end
    centers = [centers; [frameNumber, id', x, y, round(r), edges]];
end
else
    for n = 1:nFrames
    n
    posData = dlmread([directory, datafiles(n).name]);
    
    %posData = load([directory, datafiles(n).name]);
    datafiles(n).name(frameidind:frameidind+3)
    frameid = str2num(datafiles(n).name(frameidind:frameidind+3))
    %posData = new;
    frameNumber = ones(length(posData),1)*frameid;
    id = [1:length(posData)];
    posArray = [posArray; [posData(:,1), posData(:,2),  frameNumber]];
    radArray = [radArray; posData(:,3)];
    edgesArray =[edgesArray;posData(:,4)];
    
    centers = [centers; [frameNumber, id', posData(:,1), posData(:,2), posData(:,3), posData(:,4)]];
    end
end
centers
% out = fopen([directory,'centers_tracked.txt'],'w');
% fprintf(out,['frame', ',', 'particleID', ',', 'x' ',' 'y', ',','r',',','edge''\n']);
% fclose(out);
% dlmwrite([directory,'centers_tracked.txt'], centers, 'delimiter',',','-append');
param = struct('mem', 50, 'good', 0,'dim',2,'quiet',0);
linked = track(posArray, 30,param);



%carry radii information along too
radind = zeros(length(radArray),1);
edgesind = zeros(length(edgesArray),1);
for l = 1:length(radArray)
    ind = find(posArray(l,1) == linked(:,1) & posArray(l,2) ==linked(:,2));
    radind(ind) = radArray(l);
    edgesind(ind) = edgesArray(l);
end

rearrange = [linked(:,3),linked(:,4), linked(:,1), linked(:,2), radind, edgesind];
rearrange = sortrows(rearrange);
% %% 
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
dlmwrite([directory,'centers_tracked.txt'], rearrange, 'delimiter',',','-append');

if boundaryType == "annulus"
rearrange2 = rearrange
x = rearrange(:,3)
y = rearrange(:,4)

midx = 6304
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
u = uv(:,1)-400
v = uv(:,2)-400

rearrange2(:,3) = u
rearrange2(:,4) = v
out = fopen([directory,'centers_tracked_original.txt'],'w');
fprintf(out,['frame', ',', 'particleID', ',', 'x' ',' 'y', ',','r',',','edge''\n']);
fclose(out);
dlmwrite([directory,'centers_tracked_original.txt'], rearrange2, 'delimiter',',','-append');
% %% 
end
% 
%if verbose
figure;
image = dir([directory, imname(1:end-4),'warped.tif'])
pic = imread([directory, image(1).name]);
imshow(pic);
hold on;
N = unique(rearrange(:,2));
cm = colormap(parula(size(N,1))); 
for n = 1:length(N)
    ind = find(rearrange(:,2) ==n);
    plot(rearrange(ind,3), rearrange(ind,4),'Color',cm(n,:));
    hold on;
end
%colormap(parula(length(all_particleIDs)));
%end
