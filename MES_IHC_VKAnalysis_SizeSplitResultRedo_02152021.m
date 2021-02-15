clear all; clc;
cd(userpath);

ResultsFilepath = 'G:\VK foci_Doug_20210209\Hydrogels\Male\CaCl2\rep 4\Analysis\';

PixConversion = 2.73; %What is pixels:micron ratio?

VK_AreaBinsMicrons = [50 200 800 3200 12800 51200];

for vk = 1:numel(VK_AreaBinsMicrons)
    VK_AreaBinsPix(1,vk) = (sqrt(VK_AreaBinsMicrons(1,vk))*PixConversion)^2;
end

Filename = strcat(ResultsFilepath,'AnalysisResults_SizeSplit.mat');
load(Filename);

disp('Collating Results and Binning VK Areas...');
NumAreaBins = numel(VK_AreaBinsPix);

for b = 1:NumAreaBins
    if b == 1
        Results_SizeSplit.AreaFilter(1:size(Results_SizeSplit.VKAreasAll,1),b) = Results_SizeSplit.VKAreasAll<=VK_AreaBinsPix(b);
    elseif b>1 && b<NumAreaBins
        Results_SizeSplit.AreaFilter(1:size(Results_SizeSplit.VKAreasAll,1),b) = Results_SizeSplit.VKAreasAll<=VK_AreaBinsPix(b) & Results_SizeSplit.VKAreasAll>VK_AreaBinsPix(b-1);
    elseif b == NumAreaBins
        Results_SizeSplit.AreaFilter(1:size(Results_SizeSplit.VKAreasAll,1),b) = Results_SizeSplit.VKAreasAll<=VK_AreaBinsPix(b) & Results_SizeSplit.VKAreasAll>VK_AreaBinsPix(b-1);
        Results_SizeSplit.AreaFilter(1:size(Results_SizeSplit.VKAreasAll,1),b+1) = Results_SizeSplit.VKAreasAll>=VK_AreaBinsPix(b);
    else
    end
end

for b = 1:NumAreaBins+1
    if sum(Results_SizeSplit.AreaFilter(:,b)) > 0
        Results_SizeSplit.IntegratedSplitAreas(1,b) = sum(Results_SizeSplit.VKAreasAll(Results_SizeSplit.AreaFilter(1:size(Results_SizeSplit.VKAreasAll,1),b)));
    else
        Results_SizeSplit.IntegratedSplitAreas(1,b) = 0;
    end
    Results_SizeSplit.PercentofTotalArea(1,b) = Results_SizeSplit.IntegratedSplitAreas(1,b)/Results_SizeSplit.IntegratedVKAreasAllPixels*100;
end

cd(ResultsFilepath);
save('AnalysisResults_SizeSplit_NewBins.mat','Results_SizeSplit','-v7.3');