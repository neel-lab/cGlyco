function [ResidueTypeExist,StructTypeExist] = Plotheatmap(ResidueMSFileList,StructMSFileList,outputfilename,loadpath,colormaploadpath,storepath)
% Plotheatmap: Plot heatmap based on the monoshaccride and substructure
% analysis
%
% Input:
%  ResidueMSFileList: Monosacchride analysis reuslt data
%   StructMSFileList: Substructure analysis reuslt data
%     outputfilename: Output file name for heatmap
%           loadpath: Directory to load monosacchride and substructure data
%   colormaploadpath: Directory to load custom colormap data
%          storepath: Directory to store the generated heatmap
%
% Author:Yusen Zhou
% Data Lastly Updated:05/20/2020
ResidueDataMat = cell(length(ResidueMSFileList),1);
StructDataMat  = cell(length(StructMSFileList),1);
for i = 1:length(ResidueMSFileList)
    ithfilefullpath = [loadpath ResidueMSFileList{i} '.mat'];
    load(ithfilefullpath);
    ResidueDataMat{i,1} = ResidueResult;
end
for i = 1:length(StructMSFileList)
    ithfilefullpath = [loadpath StructMSFileList{i} '.mat'];
    load(ithfilefullpath);
    StructDataMat{i,1} = StructResult;
end
filespec_user = [loadpath outputfilename '.xlsx'];
try
    Excel=actxGetRunningServer('Excel.Application');
catch
    Excel = actxserver('Excel.Application');
end
Workbooks = Excel.Workbooks;
if exist(filespec_user,'file')
    Workbook = invoke(Workbooks,'Open',filespec_user);
else
    Workbook = invoke(Workbooks, 'Add');
    Workbook.SaveAs(filespec_user);
end
colormappath = [colormaploadpath 'mycolor.mat'];
load(colormappath)
[ResidueTypeExist,sheetnumber] = PlotResidueHeatMap(ResidueDataMat,mycolor,storepath,filespec_user);
StructTypeExist  = PlotStructHeatMap(StructDataMat,mycolor,storepath,sheetnumber,filespec_user);
Workbook.Save;
invoke(Excel,'Quit');
Excel.delete
clear Excel
end

function [ResidueTypeExist,sheetnumber] = PlotResidueHeatMap(ResidueDataMat,mycolor,storepath,filespec_user)
ResidueType = cell(4,1);
ResidueType{1,1} = 'Gal';
ResidueType{2,1} = 'GlcNAc';
% ResidueType{3,1} = 'GalNAc';
ResidueType{3,1} = 'Fuc';
ResidueType{4,1} = 'NeuAc';
ResidueTypeExist = ResidueType;
deleteCell = [];
sheetnumber = 0;
for i = 1 : length(ResidueType)
    MSComparisonList = checkMSList(ResidueDataMat,ResidueType{i,1});
    if(length(MSComparisonList)>1)
        MSnumber    = length(MSComparisonList);
        ValueDataS  = zeros(MSnumber,MSnumber);
        ValueDataR  = zeros(MSnumber,MSnumber);
        SpeciesA    = cell(length(MSnumber),1);
        ResidueA    = cell(length(MSnumber),1);
        RowLabels   = cell(length(MSnumber),1);
        columLabels = cell(length(MSnumber),1);
        for j = 1:MSnumber
            jthyaxis    = MSComparisonList{j}.(ResidueType{i,1}).glycanSpecies;
            SpeciesA{j,1} = jthyaxis;
            RowLabels{j,1} = MSComparisonList{j}.MSID;
            columLabels{j,1} = MSComparisonList{j}.MSID;
            for k = 1:MSnumber
                kthxaxis = MSComparisonList{k}.(ResidueType{i,1}).glycanSpecies;
                maxValue = min(100,(jthyaxis/kthxaxis));
                minValue = max(0.01,maxValue);
                ValueDataS(j,k)=log10(minValue);
            end
        end
        F = figure();
        heatmap(ValueDataS, columLabels, RowLabels,'%0.2f',...
            'Colormap',mycolor, 'Colorbar', 'true', 'TickAngle', 45,'FontSize',15,'TickFontSize',15);
        if(max(max(abs(ValueDataS)))~=0)
            colorrange = ceil(max(max(abs(ValueDataS)))/0.5)*0.5;
        else
            colorrange = 0.5;
        end
        set(gca, 'CLim', [-colorrange, colorrange]);
        colorbar('YLim',[-colorrange, colorrange]);
        TitleS = ['% relative abundance (' ResidueType{i,1} ')'];
        title(TitleS,'Fontsize',18);
        set(F,'PaperPositionMode','auto','visible','off','outerposition',[0,0,1100,800],'position', [0,0,900,600]);
        storepath1 = [storepath,'GlycanS_' ResidueType{i,1} '.jpg'];
        print(F,'-djpeg',storepath1);
        delete(F);
        for j = 1:MSnumber
            jthyaxis      = MSComparisonList{j}.(ResidueType{i,1}).glycanResidue;
            ResidueA{j,1} = jthyaxis;
            for k = 1:MSnumber
                kthxaxis = MSComparisonList{k}.(ResidueType{i,1}).glycanResidue;
                maxValue = min(100,(jthyaxis/kthxaxis));
                minValue = max(0.01,maxValue);
                ValueDataR(j,k)=log10(minValue);
            end
        end
        sheetnumber = sheetnumber+1;
        F = figure();
        heatmap(ValueDataR, columLabels, RowLabels,'%0.2f',...
            'Colormap',mycolor, 'Colorbar', 'true','TickAngle', 45,'FontSize',15,'TickFontSize',15);
         if(max(max(abs(ValueDataR)))~=0)
            colorrange = ceil(max(max(abs(ValueDataR)))/0.5)*0.5;
        else
            colorrange = 0.5;
        end
        set(gca, 'CLim', [-colorrange, colorrange]);
        colorbar('YLim',[-colorrange, colorrange]);
        TitleR = ['% relative abundance (' ResidueType{i,1} ')'];
        title(TitleR,'Fontsize',18);
        set(F,'PaperPositionMode','auto','visible','off','outerposition',[0,0,900,600],'position', [0,0,900,600]);
        storepath2 = [storepath,'GlycanR_' ResidueType{i,1} '.jpg'];
        print(F,'-djpeg',storepath2);
        delete(F);
        Excel = evalin('caller','Excel');
        warning('off','MATLAB:xlswrite:AddSheet');
        B1=cellstr(TitleS);
        C1=cellstr(TitleR);
        xlswrite1(filespec_user,B1,sheetnumber,'B1');
        xlswrite1(filespec_user,C1,sheetnumber,'C1');
        xlswrite1(filespec_user,columLabels,sheetnumber,'A2');
        xlswrite1(filespec_user,SpeciesA,sheetnumber,'B2');
        xlswrite1(filespec_user,ResidueA,sheetnumber,'C2');
        Sheets = Excel.ActiveWorkBook.Sheets;
        sheet1 = get(Sheets, 'Item', sheetnumber);
        invoke(sheet1, 'Activate')
        Sheets.Item(sheetnumber).Name = ResidueType{i,1};
        Shapes = sheet1.Shapes;
        Shapes.AddPicture(storepath1 ,0,1,40,120,900,600);
        Shapes.AddPicture(storepath2 ,0,1,40,760,900,600);
    else
        deleteCell = [deleteCell i];
    end
end
ResidueTypeExist(deleteCell) = '';
end

function StructTypeExist  = PlotStructHeatMap(StructDataMat,mycolor,storepath,sheetnumber,filespec_user)
StructType  = cell(10,1);
Structtitle = cell(10,1);
StructType{1,1}  ='Lewisx';
StructType{2,1}  ='sialylLewisx';
StructType{3,1}  ='Lewisy';
StructType{4,1}  ='LacNAc';
StructType{5,1}  ='SialylLacNAc';
StructType{6,1}  ='Hantigen';
StructType{7,1}  ='Aantigen';
StructType{8,1}  ='Bantigen';
StructType{9,1}  ='Bisecting';
StructType{10,1} ='coreFuc';
Structtitle{1,1}  ='Lewisx';
Structtitle{2,1}  ='sialylLewisx';
Structtitle{3,1}  ='Lewisy';
Structtitle{4,1}  ='LacNAc';
Structtitle{5,1}  ='SialylLacNAc';
Structtitle{6,1}  ='Hantigen';
Structtitle{7,1}  ='Aantigen';
Structtitle{8,1}  ='Bantigen';
Structtitle{9,1}  ='Bisecting';
Structtitle{10,1} ='coreFuc';
StructTypeExist  = StructType;
deleteCell = [];
for i = 1 : length(StructType)
    MSComparisonList = checkMSList(StructDataMat,StructType{i,1});
    if(length(MSComparisonList)>1)
        MSnumber    = length(MSComparisonList);
        ValueDataSp = zeros(MSnumber,MSnumber);
        ValueDataSt = zeros(MSnumber,MSnumber);
        SpeciesA     = cell(length(MSnumber),1);
        StructA     = cell(length(MSnumber),1);
        RowLabels   = cell(length(MSnumber),1);
        columLabels = cell(length(MSnumber),1);
        for j = 1:MSnumber
            jthyaxis = MSComparisonList{j}.(StructType{i,1}).glycanSpecies;
            SpeciesA{j,1} = jthyaxis;
            RowLabels{j,1} = MSComparisonList{j}.MSID;
            columLabels{j,1} = MSComparisonList{j}.MSID;
            for k = 1:MSnumber
                kthxaxis = MSComparisonList{k}.(StructType{i,1}).glycanSpecies;
                maxValue = min(100,(jthyaxis/kthxaxis));
                minValue = max(0.01,maxValue);
                ValueDataSp(j,k)=log10(minValue);
            end
        end
        F = figure();
        heatmap(ValueDataSp, columLabels, RowLabels,'%0.2f',...
            'Colormap',mycolor, 'Colorbar', 'true','TickAngle', 45,'FontSize',15,'TickFontSize',15);
        if(max(max(abs(ValueDataSp)))~=0)
            colorrange = ceil(max(max(abs(ValueDataSp)))/0.5)*0.5;
        else
            colorrange = 0.5;
        end
        set(gca, 'CLim', [-colorrange, colorrange]);
        colorbar('YLim',[-colorrange, colorrange]);
        TitleSp = ['% relative abundance of glycan with (' Structtitle{i,1} ')'];
        title(TitleSp,'Fontsize',18);
        set(F,'PaperPositionMode','auto','visible','off','outerposition',[0,0,900,600],'position', [0,0,900,600]);
        StorePath1 = [storepath,'GlycanSp_' StructType{i,1} '.jpg'];
        print(F,'-djpeg',StorePath1);
        delete(F);
        for j = 1:MSnumber
            jthyaxis = MSComparisonList{j}.(StructType{i,1}).glycanStruct;
            StructA{j,1} = jthyaxis;
            RowLabels{j,1} = MSComparisonList{j}.MSID;
            columLabels{j,1} = MSComparisonList{j}.MSID;
            for k = 1:MSnumber
                kthxaxis = MSComparisonList{k}.(StructType{i,1}).glycanStruct;
                maxValue = min(100,(jthyaxis/kthxaxis));
                minValue = max(0.01,maxValue);
                ValueDataSt(j,k)=log10(minValue);
            end
        end
        sheetnumber = sheetnumber+1;
        F = figure();
        heatmap(ValueDataSt, columLabels, RowLabels,'%0.2f',...
            'Colormap',mycolor, 'Colorbar', 'true','TickAngle', 45,'FontSize',15,'TickFontSize',15);
        if(max(max(abs(ValueDataSt)))~=0)
            colorrange = ceil(max(max(abs(ValueDataSt)))/0.5)*0.5;
        else
            colorrange = 0.5;
        end
        set(gca, 'CLim', [-colorrange, colorrange]);
        colorbar('YLim',[-colorrange, colorrange]);
        TitleSt = ['% relative abundance of glycan with (' Structtitle{i,1} ')'];
        title(TitleSt,'Fontsize',18);
        set(F,'PaperPositionMode','auto','visible','off','outerposition',[0,0,900,600],'position', [0,0,900,600]);
        StorePath2 = [storepath,'GlycanSt_' StructType{i,1} '.jpg'];
        print(F,'-djpeg',StorePath2);
        delete(F);
        
        Excel = evalin('caller','Excel');
        B1=cellstr(TitleSp);
        C1=cellstr(TitleSt);
        xlswrite1(filespec_user,B1,sheetnumber,'B1');
        xlswrite1(filespec_user,C1,sheetnumber,'C1');
        xlswrite1(filespec_user,columLabels,sheetnumber,'A2');
        xlswrite1(filespec_user,SpeciesA,sheetnumber,'B2');
        xlswrite1(filespec_user,StructA,sheetnumber,'C2');
        Sheets = Excel.ActiveWorkBook.Sheets;
        sheet1 = get(Sheets, 'Item', sheetnumber);
        invoke(sheet1, 'Activate')
        Sheets.Item(sheetnumber).Name = StructType{i,1};
        Shapes = sheet1.Shapes;
        Shapes.AddPicture(StorePath1 ,0,1,40,120,900,600);
        Shapes.AddPicture(StorePath2 ,0,1,40,760,900,600);
    else
        deleteCell = [deleteCell i];
    end
end
StructTypeExist(deleteCell) = '';
end

function MSComparisonList = checkMSList(ObjectDataMat,object)
MSComparisonList = '';
for i = 1:length(ObjectDataMat)
    relativeValue = ObjectDataMat{i,1}.(object).glycanSpecies;
    if(relativeValue>0)
        MSComparisonList{end+1}=ObjectDataMat{i,1};
    end
end
end





