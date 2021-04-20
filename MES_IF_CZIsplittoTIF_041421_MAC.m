clear all; clc;
tic
cd(userpath);

Folder = '/Filepath/';

Channels = 3; %How many fluorescent channels are in the image? (2-4, inclusive)

%% Everything else (don't touch) %%
disp('Preparing images for analysis...');
cd(Folder);
srcFiles = dir(Folder);

for dots = 1:length(srcFiles)
    if srcFiles(dots).name(1) == '.' || srcFiles(dots).isdir == 1
        isimage(dots,1) = 1;
    else isimage(dots,1) = 0; 
    end
end
ff = sum(isimage);

START = 1;
FINISH = length(srcFiles);

for f = START+ff:FINISH
    
    time(f,1).ElapsedSeconds = toc;
    
    clc
    cd(Folder);    
    
    filename = srcFiles(f).name;
    progress = ((f-ff)/(FINISH-ff)*100);
    progress2 = sprintf('Analyzing image %d of %d; %0.2f%c complete.',f-ff,(length(srcFiles)-ff),progress,'%');
    clc; disp(progress2);

    if f == START+ff
        cd(Folder); mkdir('Analysis'); cd(Folder);
    else end

if progress < 10
    disp('Estimated time remaining will display after 10% of images are analyzed...');
else
    time(f).AverageSecondsPerLoop = time(f).ElapsedSeconds/((FINISH-START+1)-(FINISH-f));
    time(f).EstimatedTotalSeconds = time(f).AverageSecondsPerLoop*(FINISH-START+1);
    time(f).EstimatedSecondsLeft = time(f).EstimatedTotalSeconds-time(f).ElapsedSeconds;
    time(f).EstimatedMinutesLeft = time(f).EstimatedSecondsLeft/60;
    time(f).EstimatedMinutesElapsed = time(f).ElapsedSeconds/60;
    estimate = sprintf('Run time: %0.2f minutes.',time(f).EstimatedMinutesElapsed);
    estimate2 = sprintf('Estimated time remaining: %0.2f minutes.',time(f).EstimatedMinutesLeft);
    disp(estimate);
    disp(estimate2);
end

AllSites = bfopen(filename);
Res = length(AllSites{1,1}{1,1});
Reshape_lin = Res*Res;
ReshapeV = [Reshape_lin 1];
ZPlanes = (length(AllSites{1,1})/Channels);
Blank = zeros(Res,Res,ZPlanes);
Sites = size(AllSites,1);

for s = 1:Sites
    if f>START+ff || s>1, clearvars Image; else end        
    for i = 1:ZPlanes
        Image.Ch1.planeIdx(i,1) = 1+(Channels*i-Channels);
        Image.Ch1.Image(:,:,i) = AllSites{s,1}{Image.Ch1.planeIdx(i,1),1};

        Image.Ch2.planeIdx(i,1) = 2+(Channels*i-Channels);
        Image.Ch2.Image(:,:,i) = AllSites{s,1}{Image.Ch2.planeIdx(i,1),1};

        if Channels>2 
            Image.Ch3.planeIdx(i,1) = 3+(Channels*i-Channels);
            Image.Ch3.Image(:,:,i) = AllSites{s,1}{Image.Ch3.planeIdx(i,1),1};
        else
        end

        if Channels>3
            Image.Ch4.planeIdx(i,1) = 4+(Channels*i-Channels);
            Image.Ch4.Image(:,:,i) = AllSites{s,1}{Image.Ch4.planeIdx(i,1),1};
        else
        end
    end
    
    Image.Full(:,:,:,1) = Image.Ch1.Image(:,:,:);
    Image.Full(:,:,:,2) = Image.Ch2.Image(:,:,:);
    if Channels>2, Image.Full(:,:,:,3) = Image.Ch3.Image(:,:,:); else end
    if Channels>3, Image.Full(:,:,:,4) = Image.Ch4.Image(:,:,:); else end

    
    if f == START+ff && s == 1
        mkdir('SplitTIFFs'); cd SplitTIFFs;
    else cd(Folder); cd SplitTIFFs;
    end
    
    disp('Saving TIFF image...');
    sitestr=num2str(s);
    if s<10, SAVENAME = strcat(filename(1:end-4),' Site 0',sitestr,'.ome.tiff');
    else SAVENAME = strcat(filename(1:end-4),' Site ',sitestr,'.ome.tiff');
    end
    bfsave(Image.Full(:,:,:,:),SAVENAME);
end    

end
