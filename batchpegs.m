%%batch process files

toptopDirectory = '/eno/cllee3/DATA/jekollme/20160708/Steps/'
folders = dir([toptopDirectory, 'step*'])

for folder=1:length(folders)
folder =folder
    topDirectory = [toptopDirectory, folders(folder).name,'/']
    imageNames = 'Step*.jpg'
    PeGSModular(topDirectory, imageNames)
end
