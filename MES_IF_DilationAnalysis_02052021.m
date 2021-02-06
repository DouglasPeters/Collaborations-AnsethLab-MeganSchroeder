clear all; clc;
tic
cd(userpath);

%% Editables %%

Folder = 'G:\For Doug analysis\IF dilation analysis_D12 gels\Males\Plus CaCl2\'; %This tells the script where to find your images (must keep the trailing slash).

Channels = 2; %How many channels are in the image? (Inputs can be 2, 3, or 4). This script assumes Channel 1 is DAPI and will use it to ID nuclei.

Ch_DAPI = 1; %Which channel is being used for nuclear segmentation (i.e., DAPI)? (1-4)
Ch1_Analysis = 0; %Do you want to analyze fluorescence in Channel 1? (Yes=1, No=0)
Ch2_Analysis = 1; %Do you want to analyze fluorescence in Channel 2? (Yes=1, No=0)
Ch3_Analysis = 0; %Do you want to analyze fluorescence in Channel 3? (Yes=1, No=0)
Ch4_Analysis = 0; %Do you want to analyze fluorescence in Channel 4? (Yes=1, No=0)

Radius_Erode = 1; %Adjust the erosion radius to trim more signal from edges of nuclei (prior to dilation). This can help clean up small unwanted pixels.
Radius_Dilate = 5; %Adjust the dilation radius that is used to identify "NEAR" (i.e., cellular) pixels around nuclei. This can make your NEAR:FAR ratio more or less stringent.

Image_Show = 1; %Do you want to render images of the dilation analysis? (Yes=1, No=0)
Ch_Show = 2; %Which channel do you want to render? (Inputs can be 2, 3, or 4).
Image_Save = 1; %Do you want to save images of the dilation analysis as .jpg files? (Yes=1, No=0)


%% Everything else (don't touch) %%
disp('Preparing images for analysis...');
cd(Folder);
srcFiles = dir(Folder);

if Ch1_Analysis == 1, RESULTS_Ch1 = table(); else end
if Ch2_Analysis == 1, RESULTS_Ch2 = table(); else end
if Ch3_Analysis == 1, RESULTS_Ch3 = table(); else end
if Ch4_Analysis == 1, RESULTS_Ch4 = table(); else end

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

    if f == START+ff,
        cd(Folder); mkdir('Analysis'); cd(Folder);
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

disp('Extracting image information...');
    Results(f-ff,1).FileName = srcFiles(f).name;
    
    I = bfopen(filename);
    Res = length(I{1,1}{1,1});
    Reshape_lin = Res*Res;
    ReshapeV = [Reshape_lin 1];
    ZPlanes = (length(I{1,1})/Channels);
    Blank = zeros(Res,Res,ZPlanes);

    for i = 1:ZPlanes
        Ch1.planeIdx(i,1) = 1+(Channels*i-Channels);
        Ch1.Image(:,:,i) = I{1,1}{Ch1.planeIdx(i,1),1};
        Ch1.Linear(1:Reshape_lin,i) = reshape(Ch1.Image(:,:,i),ReshapeV);
        Ch1.MedianImage = median(Ch1.Image,3);
        
        Ch2.planeIdx(i,1) = 2+(Channels*i-Channels);
        Ch2.Image(:,:,i) = I{1,1}{Ch2.planeIdx(i,1),1};
        Ch2.Linear(1:Reshape_lin,i) = reshape(Ch2.Image(:,:,i),ReshapeV);
        Ch2.MedianImage = median(Ch2.Image,3);
        
        if Channels>2 
            Ch3.planeIdx(i,1) = 3+(Channels*i-Channels);
            Ch3.Image(:,:,i) = I{1,1}{Ch3.planeIdx(i,1),1};
            Ch3.Linear(1:Reshape_lin,i) = reshape(Ch3.Image(:,:,i),ReshapeV);
            Ch3.MedianImage = median(Ch3.Image,3);
        else
        end
        
        if Channels>3
            Ch4.planeIdx(i,1) = 4+(Channels*i-Channels);
            Ch4.Image(:,:,i) = I{1,1}{Ch4.planeIdx(i,1),1};
            Ch4.Linear(1:Reshape_lin,i) = reshape(Ch4.Image(:,:,i),ReshapeV);
            Ch4.MedianImage = median(Ch4.Image,3);
        else
        end
    end
    
    Results(f-ff,1).Ch1IntDen = sum(Ch1.Linear(:));
    Results(f-ff,1).Ch1MeanInt = mean(Ch1.Linear(:));
    
    Results(f-ff,1).Ch2IntDen = sum(Ch2.Linear(:));
    Results(f-ff,1).Ch2MeanInt = mean(Ch2.Linear(:));
    
    if Channels > 2
        Results(f-ff,1).Ch3IntDen = sum(Ch3_lin(:));
        Results(f-ff,1).Ch3MeanInt = mean(Ch3_lin(:));    
    else end
    
    if Channels > 3
        Results(f-ff,1).Ch4IntDen = sum(Ch4_lin(:));
        Results(f-ff,1).Ch4MeanInt = mean(Ch4_lin(:));     
    else end

    erodestrel = strel('disk',Radius_Erode);
    dilatestrel = strel('disk',Radius_Dilate);
    
    if Ch_DAPI == 1, DAPI.Image = Ch1.MedianImage;
    elseif Ch_DAPI == 2, DAPI.Image = Ch2.MedianImage;
    elseif Ch_DAPI == 3, DAPI.Image = Ch3.MedianImage;
    elseif Ch_DAPI == 4, DAPI.Image = Ch4.MedianImage;
    else
        warning('The "Ch_DAPI" variable must be a value between 1 and the total number of channels.');
    end
    
    DAPI.BW1 = imbinarize(DAPI.Image);
    DAPI.BW2 = imerode(DAPI.BW1,erodestrel);
    DAPI.BW3 = imdilate(DAPI.BW2,dilatestrel);
    DAPI.Donut = logical(DAPI.BW3-DAPI.BW2);

if Ch1_Analysis == 1
    Results(f-ff,1).Ch1IntDen = sum(Ch1.MedianImage,'all');
    Ch1.CellularPixels = Ch1.MedianImage(DAPI.BW3);
    Results(f-ff,1).Ch1CellularMean = mean(Ch1.CellularPixels);
    Ch1.ECMPixels = Ch1.MedianImage(~DAPI.BW3);
    Results(f-ff,1).Ch1ECMMean = mean(Ch1.ECMPixels);
    Ch1.NuclearPixels = Ch1.MedianImage(DAPI.BW2);
    Results(f-ff,1).Ch1NuclearMean = mean(Ch1.NuclearPixels);
    Ch1.PerinuclearPixels = Ch1.MedianImage(DAPI.Donut);
    Results(f-ff,1).Ch1PerinuclearMean = mean(Ch1.PerinuclearPixels);
    Results(f-ff,1).CelltoECMRatio = Results(f-ff,1).Ch1CellularMean/Results(f-ff,1).Ch1ECMMean;
    Results(f-ff,1).NucleartoPerinuclearRatio = Results(f-ff,1).Ch1NuclearMean/Results(f-ff,1).Ch1PerinuclearMean;
else end    
    
if Ch2_Analysis == 1,
    Results(f-ff,1).Ch2IntDen = sum(Ch2.MedianImage,'all');
    Ch2.CellularPixels = Ch2.MedianImage(DAPI.BW3);
    Results(f-ff,1).Ch2CellularMean = mean(Ch2.CellularPixels);
    Ch2.ECMPixels = Ch2.MedianImage(~DAPI.BW3);
    Results(f-ff,1).Ch2ECMMean = mean(Ch2.ECMPixels);
    Ch2.NuclearPixels = Ch2.MedianImage(DAPI.BW2);
    Results(f-ff,1).Ch2NuclearMean = mean(Ch2.NuclearPixels);
    Ch2.PerinuclearPixels = Ch2.MedianImage(DAPI.Donut);
    Results(f-ff,1).Ch2PerinuclearMean = mean(Ch2.PerinuclearPixels);
    Results(f-ff,1).CelltoECMRatio = Results(f-ff,1).Ch2CellularMean/Results(f-ff,1).Ch2ECMMean;
    Results(f-ff,1).NucleartoPerinuclearRatio = Results(f-ff,1).Ch2NuclearMean/Results(f-ff,1).Ch2PerinuclearMean;
else end

if Ch3_Analysis == 1,
    Results(f-ff,1).Ch3IntDen = sum(Ch3.MedianImage,'all');
    Ch3.CellularPixels = Ch3.MedianImage(DAPI.BW3);
    Results(f-ff,1).Ch3CellularMean = mean(Ch3.CellularPixels);
    Ch3.ECMPixels = Ch3.MedianImage(~DAPI.BW3);
    Results(f-ff,1).Ch3ECMMean = mean(Ch3.ECMPixels);
    Ch3.NuclearPixels = Ch3.MedianImage(DAPI.BW2);
    Results(f-ff,1).Ch3NuclearMean = mean(Ch3.NuclearPixels);
    Ch3.PerinuclearPixels = Ch3.MedianImage(DAPI.Donut);
    Results(f-ff,1).Ch3PerinuclearMean = mean(Ch3.PerinuclearPixels);
    Results(f-ff,1).CelltoECMRatio = Results(f-ff,1).Ch3CellularMean/Results(f-ff,1).Ch3ECMMean;
    Results(f-ff,1).NucleartoPerinuclearRatio = Results(f-ff,1).Ch3NuclearMean/Results(f-ff,1).Ch3PerinuclearMean;
else end

if Ch4_Analysis == 1,
    Results(f-ff,1).Ch4IntDen = sum(Ch4.MedianImage,'all');
    Ch4.CellularPixels = Ch4.MedianImage(DAPI.BW3);
    Results(f-ff,1).Ch4CellularMean = mean(Ch4.CellularPixels);
    Ch4.ECMPixels = Ch4.MedianImage(~DAPI.BW3);
    Results(f-ff,1).Ch4ECMMean = mean(Ch4.ECMPixels);
    Ch4.NuclearPixels = Ch4.MedianImage(DAPI.BW2);
    Results(f-ff,1).Ch4NuclearMean = mean(Ch4.NuclearPixels);
    Ch4.PerinuclearPixels = Ch4.MedianImage(DAPI.Donut);
    Results(f-ff,1).Ch4PerinuclearMean = mean(Ch4.PerinuclearPixels);
    Results(f-ff,1).CelltoECMRatio = Results(f-ff,1).Ch4CellularMean/Results(f-ff,1).Ch4ECMMean;
    Results(f-ff,1).NucleartoPerinuclearRatio = Results(f-ff,1).Ch4NuclearMean/Results(f-ff,1).Ch4PerinuclearMean;
else end

%% Figure %%

if Image_Show == 1
    disp('Generating Figure...');
    B = bwboundaries(DAPI.BW3,8);
    if Ch_Show ==1, imshow(Ch1.MedianImage.*5);
    elseif Ch_Show ==2, imshow(Ch2.MedianImage.*5);
    elseif Ch_Show ==3, imshow(Ch3.MedianImage.*5);
    elseif Ch_Show ==4, imshow(Ch4.MedianImage.*5);
    else
        imshow(Ch1.MedianImage.*5); 
    end
    hold on
    for k = 1:size(B,1)
        b = B{k};
        plot(b(:,2),b(:,1),'r','LineWidth',1);
    end
    hold off
    drawnow;
    
    if Image_Save == 1
        disp('Saving Figure Image...');
        figfilename = strcat(srcFiles(f).name(1:end-4),' DilationAnalysis.jpg');
        if f-ff==1
            cd(Folder); cd Analysis; mkdir DilationAnalysisImages; cd DilationAnalysisImages;
        else
            cd(Folder); cd Analysis; cd DilationAnalysisImages;
        end
    ax = gca;
    exportgraphics(ax,figfilename,'ContentType','image','Resolution','400');
    else
    end
else
end

end

close all;
cd(Folder); cd Analysis;
save('DilationAnalysisResults.mat','Results');