function [newglycanDB,outputfilename] = glycanAbundance(MSrawdata,glycanDBFile,...
    MSdir,glycnDBdir,OverSegFilter,outputExcel,outputdir,varargin)
%glycanAbundance: Calculate the relative abundance for each potential
% glycans based on the MSdata and potential glycans, and write the result
% into a Excel file.
%
%[newglycanDB,outputfilename] = glycanAbundance(MSdatamat,glycanDBname,...
%   loadpath,staticloadpath,OverSegmentationFilter,abundanceExcelname,outputfilepath,varargin)
%
% Input:
%  MSrawdata: spectra 'peak list' along with corresponding FWHM (full width
%  at half maximum.
%
%  glycanDBname: Candidate glycan structure list including 4 fields: i) glycan structure
%  list in SGP format, ii) glycan composition, iii)monoisotopic mass and iv)
%  isotopic distribution.
%
%  MSdir: Directory to load MSdatamat.
%
%  glycnDBdir: Directory to load glycanDB.
%
%  OverSegFilter: Maximum distance between two adjacent peaks. 
%
%  outputExcel: Output excel sheet name.
%
%  outputdir: Directory to output the excel file
%
% Author: Yusen Zhou 
% Date Lastly Updated: 06/25/20
fitOption = '';
if(length(varargin)==1)
    fitOption = varargin{1};
end

MSdatamatFile = [MSrawdata '.mat'];
MSdatamatfullname = fullfile(MSdir,MSdatamatFile);
load(MSdatamatfullname,'MSdata');
peaklist = MSdata.peaklist;
FWHM     = MSdata.FWHM;
glycanDBFile = [glycanDBFile '.mat'];
glycanDBfullpath   = fullfile(glycnDBdir,glycanDBFile);
load(glycanDBfullpath,'glycanDB');
if(~isempty(fitOption))
    fitOptionFile = [fitOption, '.mat'];
    fitFunloadpath = fullfile(MSdir,fitOptionFile);
    load(fitFunloadpath,'fitPara');
    [newglycanDB,Residue,matchedpeakindex] =...
        msfraction(peaklist,FWHM,glycanDB,OverSegFilter,fitPara);
else
    [newglycanDB,Residue,matchedpeakindex] =...
        msfraction(peaklist,FWHM,glycanDB,OverSegFilter);
end
MSdata.Residue = Residue;
outputfilename = '';
if((~isempty(newglycanDB.abundance)))
    if(~isempty(outputExcel))&&(~isempty(outputdir))
        outputExcel   = [outputExcel '.xlsx'];
        filespec_user = fullfile(outputdir,outputExcel);
        MSrawdataFile = [MSrawdata 'Matched.jpg'];
        jpgname1 = fullfile(MSdir, MSrawdataFile);
        MSrawdataFile = [MSrawdata 'Residual.jpg'];
        jpgname2 = fullfile(MSdir, MSrawdataFile);
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
        plotMathchedFig(peaklist,matchedpeakindex,jpgname1);
        plotResidueFig(Residue,jpgname2);
        Workbook.Save;
        invoke(Excel,'Quit');
        Excel.delete
        clear Excel
    end 
    MSrawdataFile = [MSrawdata 'glycanDB' '.mat'];
    savepath = fullfile(MSdir, MSrawdataFile);
    outputfilename = [MSrawdata 'newglycanDB'];
    save(savepath,'newglycanDB','MSdata');
end
end

function plotMathchedFig(peaklist,matchedpeakindex,jpgname)
h=figure();
ummatchedpeamz = peaklist(:,1);
ummatchedpeamz(matchedpeakindex) = '';
ummatchedpeaInt = peaklist(:,2);
ummatchedpeaInt(matchedpeakindex) = '';
bar(ummatchedpeamz,ummatchedpeaInt,'r','EdgeColor','r');
hold on
bar(peaklist(matchedpeakindex,1),peaklist(matchedpeakindex,2),'g','EdgeColor','g');
lowerbound = floor(ummatchedpeamz(1,1)/500)*500;
upperbound = ceil(max(ummatchedpeamz(:,1))/500)*500;
set(gca,'XLim',[lowerbound,upperbound]);
set(gca,'Fontsize',20);
ylabel('% Intensity','fontsize',10)
xlabel('Mass(m/z)','fontsize',10);
set(h,'PaperPositionMode','auto','visible','on','outerposition',[0,0,900,600],'position', [0,0,900,600]);
saveas(h,jpgname);
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

function plotResidueFig(Residue,jpgname2)
h2=figure();
matchedpeamz = Residue(:,2);
matchedpeaInt = Residue(:,3);
bar(matchedpeamz,matchedpeaInt,'b','EdgeColor','b');
set(gca,'Fontsize',20);
ylabel('Residue (%)','fontsize',10)
xlabel('Mass(m/z)','fontsize',10);
set(h2,'PaperPositionMode','auto','visible','off','outerposition',[0,0,900,600],'position', [0,0,900,600]);
saveas(h2,jpgname2);
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

