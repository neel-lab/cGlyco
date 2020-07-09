function PlotheatmapStruct(glycanDB1,glycanDB2,loadpath,colormaploadpath,storepath)
% PlotheatmapStruct: Plot heatmap based on the glycan composition to
% compare two glycan profiles
%
% Input:
%          glycanDB1: Glycan structure data
%          glycanDB2: Glycan structure data
%         Sameglycan: Common glycans shared by two glycan profiles
%           loadpath: Directory to load monosacchride and substructure data
%   colormaploadpath: Directory to load custom colormap data
%          storepath: Directory to store the generated heatmap
%
% Author:Yusen Zhou
% Data Lastly Updated:05/20/2020
[~, glycanListcommon]= FindDiffGlycan(glycanDB1,glycanDB2,loadpath);
if(length(glycanListcommon)>1)
    glycanDB1File    = [glycanDB1 '.mat'];
    MS1Loadpath      = fullfile(loadpath, glycanDB1File);
    load(MS1Loadpath, 'newglycanDB');
    MSFile1relativeA = findrelabundance(glycanListcommon,newglycanDB);
    glycanDB2File    = [glycanDB2 '.mat'];
    MS2Loadpath      = fullfile(loadpath, glycanDB2File);
    load(MS2Loadpath, 'newglycanDB');
    MSFile2relativeA = findrelabundance(glycanListcommon,newglycanDB);
    plotheatmap_struct(MSFile1relativeA,MSFile2relativeA,glycanListcommon,...
        glycanDB1,glycanDB2,colormaploadpath,storepath);
end
end

function MSFilerelativeA = findrelabundance(glycanListcommon,newglycanDB)
MSFilerelativeA = cell(length(glycanListcommon),1);
for i = 1 : length(glycanListcommon)
    ithcomp = glycanListcommon{i,1};
    ismatched = arrayfun(@(x)compareComp(x,ithcomp),newglycanDB.expecGlycan);
    ithglycanA = find(ismatched);
    if(~isempty(ithglycanA))
      MSFilerelativeA{i,1} = newglycanDB.abundance(ithglycanA,1);
    end
end
end

function ismatched = compareComp(glycanComp1,glycanComp2)
ismatched = 0;
if(strcmp(glycanComp1,glycanComp2))
    ismatched = 1;
end
end

function plotheatmap_struct(MSFile1relativeA,MSFile2relativeA,glycanListcommon,...
    MSFile1,MSFile2,colormaploadpath,storepath)
datanum     = length(MSFile1relativeA);
ValueDataS  = zeros(datanum,datanum);
RowLabels   = cell(length(datanum),1);
columLabels = cell(length(datanum),1);
for i = 1:datanum
    ithyaxis = MSFile1relativeA{i,1};
    RowLabels{datanum-i+1,1} = glycanListcommon{i,1};
    columLabels{i,1} = glycanListcommon{i,1};
    for j = 1 : datanum
        jthxaxis = MSFile2relativeA{j,1};
        maxValue = min(100,(ithyaxis/jthxaxis));
        minValue = max(0.01,maxValue);
        ValueDataS(datanum-i+1,j)=log10(minValue);
    end
end
colormappath = fullfile(colormaploadpath, 'mycolor.mat');
load(colormappath, 'mycolor')
F = figure();
if(length(ValueDataS)>20)
    heatmap(ValueDataS, columLabels, RowLabels,[],...
        'Colormap',mycolor,'TickAngle',90,'TickFontsize',8);
else
    heatmap(ValueDataS, columLabels, RowLabels,'%0.2f',...
        'Colormap',mycolor,'TickAngle',90);
end
TitleS = 'log10(Y\_Axis/X\_Axis)';
title(TitleS,'Fontsize',15);
set(gca, 'CLim', [-2, 2]);
colorbar('YTickLabel',{'<-2','-1.5','-1','-0.5','0','0.5','1','1.5','>2'});
set(F,'PaperPositionMode','auto','position', [0,0,1500,1600]);
set(F,'visible','off');
StoreFile = ['Struct_based' MSFile1 '_' MSFile2 '.jpg'];
StorePath = fullfile(storepath, StoreFile);
print(F,'-djpeg',StorePath);
delete(F);
end