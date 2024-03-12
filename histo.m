
topDirectory = '/eno/cllee3/DATA/jekollme/20160711/Steps/step09/'
%topDirectory = './DATA/test/'
files = '*joAdjacencyAbs.dlm'

datafiles = dir([topDirectory, files])
W = []
for frame = 1:length(datafiles)
    F = load([topDirectory, datafiles(frame).name]);
    data = F(~isnan(F));
    W = [W;data];
end

length(nonzeros(W))
edges = 10.^(-5:0.1:2);
[N,edges] = histcounts(W,edges,'Normalization','countdensity');
g = histogram('BinEdges',edges,'BinCounts',N);
set(gca, "Xscale", "log")
hold on;
jofiles = '*newtonized2.mat'

jodatafiles = dir([topDirectory, jofiles])
joW = []
for frame = 1:length(jodatafiles)
    F = load([topDirectory, jodatafiles(frame).name]);
    forces = extractfield(F.pres,'forces');
    joW = [joW, forces];
end

length(nonzeros(joW))
edges = 10.^(-5:0.1:2);
[joN,joedges] = histcounts(joW,edges,'Normalization','countdensity');
r = histogram('BinEdges',joedges,'BinCounts',joN);
set(gca, "Xscale", "log")

saveas(gcf,  '/eno/cllee3/DATA/jekollme/20160711/step09histo.jpg')
