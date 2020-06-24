function [newglycanDB,outputfilename] = glycanAbundance(MSrawdata,glycanDB,...
    loadpath,staticloadpath,OverSegmentationFilter,abundanceExcelname,outputfilepath,varargin)
% glycanAbundance: Calculate the relative abundance for each potential
% glycans based on the MSdata and potential glycans, and write the result
% into a Excel file.
%
%[newglycanDB,outputfilename] = glycanAbundance(MSdatamat,glycanDBname,...
%   loadpath,staticloadpath,OverSegmentationFilter,abundanceExcelname,outputfilepath,varargin)
%
% Input:
%  MSdata: spectra 'peak list' along with corresponding FWHM (full width
%  at half maximum.
%
%  glycanDBname: Candidate glycan structure list including 4 fields: i) glycan structure
%  list in SGP format, ii) glycan composition, iii)monoisotopic mass and iv)
%  isotopic distribution.
%
%  loadpath: Directory to load MSdatamat.
%
%  staticloadpath: Directory to load glycanDB.
%
%  OverSegmentationFilter: Maximum distance between two adjacent peaks. 
%
%  abundanceExcelname: Output excel sheet name.
%
%  outputfilepath: Directory to output the excel file
%
% Author: Yusen Zhou 
% Date Lastly Updated: 05/18/20
fitOption = '';
if(length(varargin)==1)
    fitOption = varargin{1};
end

MSdatamatFile = [MSrawdata '.mat'];
MSdatamatfullname = fullfile(loadpath,MSdatamatFile);
load(MSdatamatfullname);
peaklist = MSdata.peaklist;
pfwhh    = MSdata.pfwhh;
glycanDB = [glycanDB '.mat'];
glycanDBfullpath   = fullfile(staticloadpath,glycanDB);
load(glycanDBfullpath);
if(~isempty(fitOption))
    fitFunloadpath = [loadpath fitOption '.mat'];
    load(fitFunloadpath);
    [newglycanDB,matchedpeakindex] =...
        msfraction(peaklist,pfwhh,glycanDB,OverSegmentationFilter,fitPara);
else
    [newglycanDB,matchedpeakindex] =...
        msfraction(peaklist,pfwhh,glycanDB,OverSegmentationFilter);
end

outputfilename = '';
if((~isempty(newglycanDB.abundance)))
    if(~isempty(abundanceExcelname))&&(~isempty(outputfilepath))
        filespec_user = [outputfilepath abundanceExcelname '.xlsx'];
        jpgname1 = [loadpath MSrawdata 'Matched.jpg'];
        jpgname2 = [loadpath MSrawdata 'Residual.jpg'];
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
        write2Excel(filespec_user,newglycanDB.abundance,newglycanDB.expecGlycan,newglycanDB.monoisomw);
        plotMathchedFig(peaklist,matchedpeakindex,jpgname1);%wfMatch   = plotMathchedFig(peaklist,matchedpeakindex,jpgname1);
        plotResidueFig(newglycanDB.Residue,jpgname2);%wfResidue = plotResidueFig(newglycanDB.Residue,jpgname2);
        Workbook.Save;
        invoke(Excel,'Quit');
        Excel.delete
        clear Excel
    end 
    savepath = [loadpath MSrawdata 'glycanDB' '.mat'];
    outputfilename = [MSrawdata 'newglycanDB'];
    save(savepath,'newglycanDB');
end
end

function plotMathchedFig(peaklist,matchedpeakindex,jpgname)%wfResidual = plotMathchedFig(peaklist,matchedpeakindex,jpgname)
h=figure();
ummatchedpeamz = peaklist(:,1);
ummatchedpeamz(matchedpeakindex) = '';
ummatchedpeaInt = peaklist(:,2);
ummatchedpeaInt(matchedpeakindex) = '';
% bar(ummatchedpeamz,ummatchedpeaInt,'r','EdgeColor','r');
% hold on
bar(peaklist(matchedpeakindex,1),peaklist(matchedpeakindex,2),'g','EdgeColor','g');
lowerbound = floor(ummatchedpeamz(1,1)/500)*500;
upperbound = ceil(max(ummatchedpeamz(:,1))/500)*500;
set(gca,'XLim',[lowerbound,upperbound]);
set(gca,'Fontsize',20);
ylabel('% Intensity','fontsize',10)
xlabel('Mass(m/z)','fontsize',10);
set(h,'PaperPositionMode','auto','visible','on','outerposition',[0,0,900,600],'position', [0,0,900,600]);
saveas(h,jpgname);
% wfResidual = webfigure(h);
Excel = evalin('caller','Excel');
Sheets = Excel.ActiveWorkBook.Sheets;
sheet1 = get(Sheets, 'Item', 1);
invoke(sheet1, 'Activate');
Sheets.Item(1).Name = 'glycan';
sheet1.invoke('Pictures').Insert(jpgname);
Sheets.Add([], Sheets.Item(Sheets.Count));
Sheets.Add([], Sheets.Item(Sheets.Count));
delete(h);
end

function plotResidueFig(Residue,jpgname2)%wfResidual = plotResidueFig(Residue,jpgname2)
h2=figure();
matchedpeamz = Residue(:,2);
matchedpeaInt = Residue(:,3);
bar(matchedpeamz,matchedpeaInt,'b','EdgeColor','b');
set(gca,'Fontsize',20);
ylabel('Residue (%)','fontsize',10)
xlabel('Mass(m/z)','fontsize',10);
set(h2,'PaperPositionMode','auto','visible','off','outerposition',[0,0,900,600],'position', [0,0,900,600]);
saveas(h2,jpgname2);
% wfResidual = webfigure(h2);
Excel = evalin('caller','Excel');
Sheets = Excel.ActiveWorkBook.Sheets;
sheet1 = get(Sheets, 'Item', 1);
invoke(sheet1, 'Activate');
Sheets.Item(1).Name = 'glycan';
sheet1.invoke('Pictures').Insert(jpgname2);
Sheets.Add([], Sheets.Item(Sheets.Count));
Sheets.Add([], Sheets.Item(Sheets.Count));
delete(h2);
end

function write2Excel(filespec_user,relativeabundance,expecGlycan,monoisomw)
    A1=cellstr('Composition');
    B1=cellstr('Monoisotopic mass');
    C1=cellstr('Fraction');
    Excel = evalin('caller','Excel');
    xlswrite1(filespec_user,expecGlycan,1,'A2');
    xlswrite1(filespec_user, monoisomw,1,'B2');
    xlswrite1(filespec_user,relativeabundance,1,'C2');
    xlswrite1(filespec_user,A1,1,'A1');
    xlswrite1(filespec_user,B1,1,'B1');
    xlswrite1(filespec_user,C1,1,'C1'); 
end

