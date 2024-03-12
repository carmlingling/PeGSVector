%calculates stress and strain in granular packings taken from esbilili
%datadryad data

topDirectory = ['/eno/cllee3/DATA/esbilili/uniaxial/img3/']
fileNames = 'img*solved.mat'
px2m = 0.00019853
files = dir([topDirectory,fileNames]);
nFrames = length(files);
stress = zeros(nFrames,3);
strain = zeros(nFrames, 1);
for n = 1:nFrames
    load([topDirectory, files(n).name]);
    new = getSigma(new, px2m );
    sigmas = cat(1,new.sigma);
    Sigmas = getGlobalSigma(new, px2m);
    stress(n, 1) = str2num(files(n).name(6:9))
    g2 = extractfield(new, 'g2');
    stress(n, 3) = sum(g2);
    stress(n,2) = (Sigmas(2, 2)-Sigmas(1, 1))./(Sigmas(1,1)+Sigmas(2,2));
end
writematrix(stress, [topDirectory, 'stress.txt'])    

plot(stress)