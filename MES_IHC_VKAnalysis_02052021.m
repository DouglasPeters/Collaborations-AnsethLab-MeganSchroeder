clear all; clc;
cd(userpath);
tic
%% Editables %%
Folder = 'G:\For Doug analysis\VK particle analysis_humans\20x\MalE\';

BotCrop = 16; %Number of pixels to remove from the bottom of the image to avoid the weird camera error.
Fig_Show = 0; %Do you want to show an image of the foci segmentation analysis? (1=yes,0=no)
Fig_Save = 0; %Do you want to save an image of the foci segmentation analysis? (1=yes,0=no)

VK_StdDevThresh = 3; %How far from the average intensity does a pixel value need to be to pass threshold?
VK_MinArea = 5; %Minimum number of pixels for VK spot to be detected.
VK_AreaSplitSize = 20000; %Area threshold where big and small are split.
VK_EccentricityMax = 0.99; %How non-circular can VK spots be? This helps filter out out-of-focus cell edges. (0=line,1=circle)

%% Analysis Pre-Analysis and Metadata (Don't touch) %%

cd(Folder);
srcFiles = dir('*.jpg');

START = 1; %At what image do you want to start analysis?
FINISH = length(srcFiles); %At what image do you want to end analysis?
%Area_Distro = cell(length(srcFiles),1);

warning('off','all')
for g = START:FINISH
    clc
    disp('Calculating mean RGB levels for data set...')
    filename = strcat(Folder,srcFiles(g).name);
    I_RGB = imread(filename);
    [ResY ResX ResZ] = size(I_RGB,[1 2 3]);
    I_Red = I_RGB(1:end-BotCrop,:,1); I_Red_inv = imcomplement(I_Red);
    I_Green = I_RGB(1:end-BotCrop,:,2); I_Green_inv = imcomplement(I_Green);
    I_Blue = I_RGB(1:end-BotCrop,:,3); I_Blue_inv = imcomplement(I_Blue);
    Threshold.ImageAverages.Red(g,1) = mean(I_Red_inv(:));
    Threshold.ImageAverages.Green(g,1) = mean(I_Green_inv(:));
    Threshold.ImageAverages.Blue(g,1) = mean(I_Blue_inv(:));
end

Threshold.AllImageMean.Red = mean(Threshold.ImageAverages.Red(:));
Threshold.AllImageMean.Green = mean(Threshold.ImageAverages.Green(:));
Threshold.AllImageMean.Blue = mean(Threshold.ImageAverages.Blue(:));
Threshold.AllImageSTD.Red = std(Threshold.ImageAverages.Red(:,1))*VK_StdDevThresh;
Threshold.AllImageSTD.Green = std(Threshold.ImageAverages.Green(:,1))*VK_StdDevThresh;
Threshold.AllImageSTD.Blue = std(Threshold.ImageAverages.Blue(:,1))*VK_StdDevThresh;
Threshold.AllImageThresh.Red = Threshold.AllImageMean.Red+Threshold.AllImageSTD.Red;
Threshold.AllImageThresh.Green = Threshold.AllImageMean.Green+Threshold.AllImageSTD.Green;
Threshold.AllImageThresh.Blue = Threshold.AllImageMean.Blue+Threshold.AllImageSTD.Blue;

%% Analysis of Individual Images %%

for f = START:FINISH
    time(f,1).ElapsedSeconds = toc;
    
    clc    
    progress = (f/length(srcFiles)*100);
    progress2 = sprintf('Analyzing image %d of %d; %0.2f%c complete.',f,length(srcFiles),progress,'%');
    disp(progress2)
    filename = strcat(Folder,srcFiles(f).name);
    
    if f == START
        cd(Folder); mkdir('Analysis'); cd(Folder);
        if Fig_Show == 1, figure,
        else end
    else end
    
    if progress < 10,
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
    
    Results(f,1).FileName = srcFiles(f).name;
    Results(f,1).Settings.VKMinArea = VK_MinArea;
    Results(f,1).Settings.VKStdDevThresh = VK_StdDevThresh;
    Results(f,1).Settings.VKAreaSplitSize = VK_AreaSplitSize;
    Results(f,1).Settings.VKEccentricityMax = VK_EccentricityMax;
    
    Results(f,1).Settings.RedThresh = Threshold.AllImageThresh.Red;
    Results(f,1).Settings.GreenThresh = Threshold.AllImageThresh.Green;
    Results(f,1).Settings.BlueThresh = Threshold.AllImageThresh.Blue;
    
    Image.RGB = imread(filename);
    
    [ResY,ResX,ResZ] = size(Image.RGB,[1 2 3]);
    
    Image.Red = Image.RGB(1:end-BotCrop,:,1); Image.RedInv = imcomplement(Image.Red);
    Image.Green = Image.RGB(1:end-BotCrop,:,2); Image.GreenInv = imcomplement(Image.Green);
    Image.Blue = Image.RGB(1:end-BotCrop,:,3); Image.BlueInv = imcomplement(Image.Blue);
    Image.BW = zeros([ResY-BotCrop ResX]);
    
    for y = 1:ResY-BotCrop
        for x = 1:ResX
            if Image.RedInv(y,x)>Results(f,1).Settings.RedThresh && Image.GreenInv(y,x)>Results(f,1).Settings.GreenThresh && Image.BlueInv(y,x)>Results(f,1).Settings.BlueThresh
                Image.BW(y,x) = 1;
            else
                Image.BW(y,x) = 0;
            end
        end
    end
    
    Image.BW2 = bwareaopen(Image.BW,VK_MinArea);
    Image.BW3 = bwpropfilt(Image.BW2,'Eccentricity',[0 VK_EccentricityMax]);
    Image.cc = bwconncomp(Image.BW3,4);
    Image.regionprops = regionprops(Image.cc,'Area','Circularity','Solidity');    
    
    if f > START, clearvars VKAreas_CurrentImage; else end
    
    for a = 1:Image.cc.NumObjects
        VKAreas_CurrentImage(a,1) = Image.regionprops(a).Area;
    end
    
    Results(f,1).ObjectNumber = Image.cc.NumObjects;
    Results(f,1).ObjectAreas(1:length(VKAreas_CurrentImage)) = VKAreas_CurrentImage;
    Results(f,1).MeanArea = mean(VKAreas_CurrentImage,'all');
    
    if f == START
        VKAreas_All(1:length(VKAreas_CurrentImage),1) = VKAreas_CurrentImage;
    else
        VKAreas_All(end+1:end+length(VKAreas_CurrentImage),1) = VKAreas_CurrentImage;
    end
    
    %%% Visualization of segmentation %%%
if Fig_Show == 1
    disp('Generating Figure...');
    disp('Red-bordered objects are included in analysis.');
    disp('Green-bordered objects are those removed due to size/shape/etc.');
    B = bwboundaries(Image.BW,8);
    BB = bwboundaries(Image.BW3,8);
    imshow(Image.RGB(1:end-15,:,:));
    hold on
    for k = 1:size(B,1)
        b = B{k};
        plot(b(:,2),b(:,1),'g','LineWidth',1); %Green-bordered objects were NOT INCLUDED in analysis.
    end
    for kk = 1:size(BB,1)
        bb = BB{kk};
        plot(bb(:,2),bb(:,1),'r','LineWidth',1); %Red-bordered objects were INCLUDED in analysis.
    end
    hold off
    
    if Fig_Save == 1
        ('Saving Figure Image...');
        for slash = 1:length(strfind(filename,'\'))
            if slash == 1, filename2 = extractAfter(filename,'\');
            else
                filename2 = extractAfter(filename2,'\');
            end
        end
        filename2 = filename2(1:end-4);
        filename3 = append(filename2,'.jpg');

        if f==START, cd(Folder); cd Analysis; mkdir VKFociSegmentation; cd VKFociSegmentation;
        else cd(Folder); cd Analysis; cd VKFociSegmentation; end
        ax = gca;
        exportgraphics(ax,filename3,'ContentType','image','Resolution','400');
    else
    end
else 
end    
end

disp('Collating Results and Binning VK Areas...');
Results_SizeSplit.VKAreasAll = VKAreas_All;
Results_SizeSplit.BigFilter = VKAreas_All>=VK_AreaSplitSize;
Results_SizeSplit.SmallFilter = VKAreas_All<VK_AreaSplitSize;
Results_SizeSplit.BigAreas = VKAreas_All(Results_SizeSplit.BigFilter);
Results_SizeSplit.SmallAreas = VKAreas_All(Results_SizeSplit.SmallFilter);
Results_SizeSplit.SumAreaBig = sum(Results_SizeSplit.BigAreas);
Results_SizeSplit.SumAreaSmall = sum(Results_SizeSplit.SmallAreas);
Results_SizeSplit.TotalArea = Results_SizeSplit.SumAreaSmall + Results_SizeSplit.SumAreaBig;
close all

cd(Folder); cd Analysis;
save('AnalysisResults.mat','Results','-v7.3');
save('AnalysisResults_SizeSplit.mat','Results_SizeSplit','-v7.3');
