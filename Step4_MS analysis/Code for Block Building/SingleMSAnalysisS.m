function SingleMSAnalysisS(GlycanDB,outputfilename,loadpath)
% SingleMSAnalysisS: Calculate relative abundance of different
% substructures 
%
%
%
% Author: Yusen Zhou
% Date Lastly Updated: 05/18/2020
% Load the glycanDB
DBfilepath = [loadpath GlycanDB '.mat'];
load(DBfilepath)
% calculate 1st MS subStructure relativeabundance(lewis x/a, sialyl-Lewis x/a, 
%              Lewis y/b, H_Antigen, A_Antigen, B_Antigen, Bisecting, coreFuc)
Lx_a       = calculatebyStr(newglycanDB,'Lx_a');
sialylLx_a = calculatebyStr(newglycanDB,'sialylLx_a');
Ly_b       = calculatebyStr(newglycanDB,'Ly_b');
LacNAc     = calculatebyStr(newglycanDB,'LacNAc');
Sia_LacNAc = calculatebyStr(newglycanDB,'Sia_LacNAc');
Hantigen   = calculatebyStr(newglycanDB,'Hantigen');
Aantigen   = calculatebyStr(newglycanDB,'Aantigen');
Bantigen   = calculatebyStr(newglycanDB,'Bantigen');
Bisecting  = calculatebyStr(newglycanDB,'Bisecting');
coreFuc    = calculatebyStr(newglycanDB,'coreFuc');
% Export to the Excel
filespec_user = [loadpath outputfilename '.xlsx'];
jpgnameSp = [loadpath outputfilename 'StructSpAbundance.jpg'];
jpgnameSt = [loadpath outputfilename 'StructStAbundance.jpg'];
try
    Excel=actxGetRunningServer('Excel.Application');
catch
    Excel = actxserver('Excel.Application');
end
Workbooks = Excel.Workbooks;
if exist(filespec_user,'file');
    Workbook = invoke(Workbooks,'Open',filespec_user);
else
    Workbook = invoke(Workbooks, 'Add');
    Workbook.SaveAs(filespec_user);
end

write2Excel(Lx_a,sialylLx_a,Ly_b,LacNAc,Sia_LacNAc,Hantigen,Aantigen,Bantigen,Bisecting,coreFuc,...
    filespec_user);
% plot the bar figure
plotbarfigure_Residue(Lx_a,sialylLx_a,Ly_b,LacNAc,Sia_LacNAc,Hantigen,Aantigen,Bantigen,Bisecting,coreFuc,...
    jpgnameSp,jpgnameSt);
Workbook.Save;
invoke(Excel,'Quit');
Excel.delete 
clear Excel
end

function StructAbundance  = calculatebyStr(newglycanDB,StructStr)
StructAbundance = struct('glycanSpecies',[],'glycanStruct',[]);
glycanStructStr = newglycanDB.glycanexpec;
for i = 1 : length(glycanStructStr)
    ithglycan = glycanStructStr{i};
    ismatch = 0;
    SpAbundance = 0;
    StAbundance = 0;
    if(ischar(ithglycan))
        [ismatch,NumberOfStruct] = calbyDiffStruct(ithglycan,StructStr,ismatch);
        if(ismatch~=0)
            SpAbundance = SpAbundance+newglycanDB.abundance(i);
            StAbundance = StAbundance+newglycanDB.abundance(i)*NumberOfStruct;
        end
    elseif(isa(ithglycan,'cell'))
        for j = 1 : length(ithglycan)
            str = ithglycan{j};
            [ismatch,NumberOfStruct] = calbyDiffStruct(str,StructStr,ismatch);
            if(ismatch~=0)
                SpAbundance = SpAbundance + ((newglycanDB.abundance(i))/length(ithglycan));
                StAbundance = StAbundance + ((newglycanDB.abundance(i)*NumberOfStruct)/length(ithglycan));
            end
        end
    end
end
SpAbundance = roundn(SpAbundance*100,-1);
StAbundance = roundn(StAbundance,-1);
StructAbundance.glycanSpecies = SpAbundance;
StructAbundance.glycanStruct  = StAbundance;
end

function [ismatch,NumberOfStruct] = calbyDiffStruct(ithglycan,StructStr,ismatch)
NumberOfStruct = 0;
if(strcmp(StructStr,'Lx_a'))
    if(~isempty(strfind(ithglycan,'{n{f}{h}}'))||~isempty(regexp(ithglycan,'{n_b\d{f_a\d}{h_b\d}}','once')))
        ismatch = ismatch+1;
        NumberOfStruct = max(length(strfind(ithglycan,'{n{f}{h}}')),length(regexp(ithglycan,'{n_b\d{f_a\d}{h_b\d}}')));
    end
elseif(strcmp(StructStr,'sialylLx_a'))
    if(~isempty(strfind(ithglycan,'{n{f}{h{s}}}'))||~isempty(regexp(ithglycan,'{n_b\d{f_a\d}{h_b\d{s_a3}}}','once')))
        ismatch = ismatch+1;
        NumberOfStruct = max(length(strfind(ithglycan,'{n{f}{h{s}}}')),length(regexp(ithglycan,'{n_b\d{f_a\d}{h_b\d{s_a3}}}')));
    end
elseif(strcmp(StructStr,'Ly_b'))
    if(~isempty(strfind(ithglycan,'{n{f}{h{f}}}'))||~isempty(regexp(ithglycan,'{n_b\d{f_a\d}{h_b\d{f_a2}}}','once')))
        ismatch = ismatch+1;
        NumberOfStruct = max(length(strfind(ithglycan,'{n{f}{h{f}}}')),length(regexp(ithglycan,'{n_b\d{f_a\d}{h_b\d{f_a2}}}')));
    end
elseif(strcmp(StructStr,'LacNAc'))
    % avoid bisecting structure.
    if(~isempty(strfind(ithglycan,'{n{h{n{'))||~isempty(strfind(ithglycan,'{n{h}'))...
            ||~isempty(regexp(ithglycan,'{n_b\d{h_b\d}','once'))||~isempty(regexp(ithglycan,'{n_b\d{h_b\d{n_b\d{','once')))
        ismatch = ismatch+1;
        NumberOfStruct = max(length(strfind(ithglycan,'{n{h{n{'))+length(strfind(ithglycan,'{n{h}'))...
            ,length(regexp(ithglycan,'{n_b\d{h_b\d}'))+length(regexp(ithglycan,'{n_b\d{h_b\d{n_b\d{')));
    end
elseif(strcmp(StructStr,'Sia_LacNAc'))
    if(~isempty(strfind(ithglycan,'{n{h{s}}}'))||~isempty(regexp(ithglycan,'{n_b\d{h_b\d{s_a3}}}','once')))
        ismatch = ismatch+1;
        NumberOfStruct = max(length(strfind(ithglycan,'{n{h{s}}}')),length(regexp(ithglycan,'{n_b\d{h_b\d{s_a3}}}')));
    end
elseif(strcmp(StructStr,'Hantigen'))
    if(~isempty(strfind(ithglycan,'{n{h{f}}}'))||~isempty(regexp(ithglycan,'{n_b\d{h_b\d{f_a2}}}','once')))
        ismatch = ismatch+1;
        NumberOfStruct = max(length(strfind(ithglycan,'{n{h{f}}}')),length(regexp(ithglycan,'{n_b\d{h_b\d{f_a2}}}')));
    end
elseif(strcmp(StructStr,'Aantigen'))
    if(~isempty(strfind(ithglycan,'{n{h{f}{n}}}'))||~isempty(regexp(ithglycan,'{n_b\d{h_b\d{f_a2}{n_a3}}}','once')))
        ismatch = ismatch+1;
        NumberOfStruct = max(length(strfind(ithglycan,'{n{h{f}{n}}}')),length(regexp(ithglycan,'{n_b\d{h_b\d{f_a2}{n_a3}}}')));
    end
elseif(strcmp(StructStr,'Bantigen'))
    if(~isempty(strfind(ithglycan,'{n{h{f}{h}}}'))||~isempty(regexp(ithglycan,'{n_b\d{h_b\d{f_a2}{h_a3}}}','once')))
        ismatch = ismatch+1;
        NumberOfStruct = max(length(strfind(ithglycan,'{n{h{f}{h}}}')),length(regexp(ithglycan,'{n_b\d{h_b\d{f_a2}{h_a3}}}')));
    end
elseif(strcmp(StructStr,'Bisecting'))
    if(~isempty(strfind(ithglycan,'{n{n{h{n}{h'))||~isempty(strfind(ithglycan,'{n_b{n_b4{h_b4{n_b4}{h_a3'))...
            ||~isempty(strfind(ithglycan,'{n{f}{n{h{n}{h'))||~isempty(strfind(ithglycan,'{n_b{f_a6}{n_b4{h_b4{n_b4}{h_a3')))
        ismatch = ismatch+1;
        NumberOfStruct = max(length(strfind(ithglycan,'{n{n{h{n}{h'))+length(strfind(ithglycan,'{n{f}{n{h{n}{h'))...
            ,length(regexp(ithglycan,'{n_b{n_b4{h_b4{n_b4}{h_a3'))+length(regexp(ithglycan,'{n_b{f_a6}{n_b4{h_b4{n_b4}{h_a3')));
    end
elseif(strcmp(StructStr,'coreFuc'))
    if(~isempty(strfind(ithglycan,'{n{f}{n{h{h'))||~isempty(strfind(ithglycan,'{n_b{f_a6}{n_b4{h_b4{h_a3'))...
            ||~isempty(strfind(ithglycan,'{n{f}{n{h{n}{h'))||~isempty(strfind(ithglycan,'{n_b{f_a6}{n_b4{h_b4{n_b4}{h_a3')))
       ismatch = ismatch+1;
       NumberOfStruct = max(length(strfind(ithglycan,'{n{f}{n{h{h'))+length(strfind(ithglycan,'{n{f}{n{h{n}{h'))...
            ,length(regexp(ithglycan,'{n_b{f_a6}{n_b4{h_b4{h_a3'))+length(regexp(ithglycan,'{n_b{f_a6}{n_b4{h_b4{n_b4}{h_a3')));
    end
end
end

function plotbarfigure_Residue(Lx_a,sialylLx_a,Ly_b,LacNAc,Sia_LacNAc,Hantigen,Aantigen,Bantigen,Bisecting,coreFuc,...
    jpgnameSp,jpgnameSt)

Y1 = [Lx_a.glycanSpecies;sialylLx_a.glycanSpecies;...
    Ly_b.glycanSpecies;LacNAc.glycanSpecies;...
    Sia_LacNAc.glycanSpecies;Hantigen.glycanSpecies;...
    Aantigen.glycanSpecies;Bantigen.glycanSpecies;...
    Bisecting.glycanSpecies;coreFuc.glycanSpecies];
figure
bar(Y1);
set(gcf,'PaperPositionMode','auto','visible','off','outerposition',[0,0,1800,1500],'position', [0,0,1500,1200]);
text(0.8,Lx_a.glycanSpecies-0.1*Lx_a.glycanSpecies,num2str(Lx_a.glycanSpecies),'Fontsize',15,'color','w');
text(1.8,sialylLx_a.glycanSpecies-0.1*sialylLx_a.glycanSpecies,num2str(sialylLx_a.glycanSpecies),'Fontsize',15,'color','w');
text(2.8,Ly_b.glycanSpecies-0.1*Ly_b.glycanSpecies,num2str(Ly_b.glycanSpecies),'Fontsize',15,'color','w');
text(3.8,LacNAc.glycanSpecies-0.1*LacNAc.glycanSpecies,num2str(LacNAc.glycanSpecies),'Fontsize',15,'color','w');
text(4.8,Sia_LacNAc.glycanSpecies-0.1*Sia_LacNAc.glycanSpecies,num2str(Sia_LacNAc.glycanSpecies),'Fontsize',15,'color','w');
text(5.8,Hantigen.glycanSpecies-0.1*Hantigen.glycanSpecies,num2str(Hantigen.glycanSpecies),'Fontsize',15,'color','w');
text(6.8,Aantigen.glycanSpecies-0.1*Aantigen.glycanSpecies,num2str(Aantigen.glycanSpecies),'Fontsize',15,'color','w');
text(7.8,Bantigen.glycanSpecies-0.1*Bantigen.glycanSpecies,num2str(Bantigen.glycanSpecies),'Fontsize',15,'color','w');
text(8.8,Bisecting.glycanSpecies-0.1*Bisecting.glycanSpecies,num2str(Bisecting.glycanSpecies),'Fontsize',15,'color','w');
text(9.8,coreFuc.glycanSpecies-0.1*coreFuc.glycanSpecies,num2str(coreFuc.glycanSpecies),'Fontsize',15,'color','w');
set(gca,'xtick',0:1:10);
set(gca,'xticklabel',{'0' 'LeX/a' 'sia_LeX/a' 'Ley/b'  'LacNAc' 'Sia_LacNAc'...
    'H_Antigen' 'A_Antigen' 'B_Antigen' 'Bisecting' 'coreFuc'},'Fontsize',15);
leg = legend('% relative abundance of Glycan (structure)');
set(leg,'Fontsize',20);
print(gcf,'-djpeg',jpgnameSp);
Y2 = [Lx_a.glycanStruct;sialylLx_a.glycanStruct;...
    Ly_b.glycanStruct;LacNAc.glycanStruct;...
    Sia_LacNAc.glycanStruct;Hantigen.glycanStruct;...
    Aantigen.glycanStruct;Bantigen.glycanStruct;...
    Bisecting.glycanStruct;coreFuc.glycanStruct];
figure
bar(Y2);
set(gcf,'PaperPositionMode','auto','visible','off','outerposition',[0,0,1800,1500],'position', [0,0,1500,1200]);
text(0.8,Lx_a.glycanStruct-0.1*Lx_a.glycanStruct,num2str(Lx_a.glycanStruct),'Fontsize',15,'color','w');
text(1.8,sialylLx_a.glycanStruct-0.1*sialylLx_a.glycanStruct,num2str(sialylLx_a.glycanStruct),'Fontsize',15,'color','w');
text(2.8,Ly_b.glycanStruct-0.1*Ly_b.glycanStruct,num2str(Ly_b.glycanStruct),'Fontsize',15,'color','w');
text(3.8,LacNAc.glycanStruct-0.1*LacNAc.glycanStruct,num2str(LacNAc.glycanStruct),'Fontsize',15,'color','w');
text(4.8,Sia_LacNAc.glycanStruct-0.1*Sia_LacNAc.glycanStruct,num2str(Sia_LacNAc.glycanStruct),'Fontsize',15,'color','w');
text(5.8,Hantigen.glycanStruct-0.1*Hantigen.glycanStruct,num2str(Hantigen.glycanStruct),'Fontsize',15,'color','w');
text(6.8,Aantigen.glycanStruct-0.1*Aantigen.glycanStruct,num2str(Aantigen.glycanStruct),'Fontsize',15,'color','w');
text(7.8,Bantigen.glycanStruct-0.1*Bantigen.glycanStruct,num2str(Bantigen.glycanStruct),'Fontsize',15,'color','w');
text(8.8,Bisecting.glycanStruct-0.1*Bisecting.glycanStruct,num2str(Bisecting.glycanStruct),'Fontsize',15,'color','w');
text(9.8,coreFuc.glycanStruct-0.1*coreFuc.glycanStruct,num2str(coreFuc.glycanStruct),'Fontsize',15,'color','w');
set(gca,'xtick',0:1:10);
set(gca,'xticklabel',{'0' 'LeX/a' 'sia_LeX/a' 'Ley/b'  'LacNAc' 'Sia_LacNAc'...
    'H_Antigen' 'A_Antigen' 'B_Antigen' 'Bisecting' 'coreFuc'},'Fontsize',15);
leg = legend('% Average # of GlycanStruct per Glycan');
set(leg,'Fontsize',20);
print(gcf,'-djpeg',jpgnameSt);

Excel = evalin('caller','Excel');
Sheets = Excel.ActiveWorkBook.Sheets;
sheet1 = get(Sheets, 'Item', 3);
invoke(sheet1, 'Activate');
Sheets.Item(3).Name = 'sub-structure';
Shapes = sheet1.Shapes;
Shapes.AddPicture(jpgnameSp ,0,1,40,200,1400,700);
Shapes.AddPicture(jpgnameSt ,0,1,40,950,1400,700);
delete(gcf);
end

function write2Excel(Lx_a,sialylLx_a,Ly_b,LacNAc,Sia_LacNAc,Hantigen,Aantigen,Bantigen,Bisecting,coreFuc,...
    filespec_user)
warning('off','MATLAB:xlswrite:AddSheet');
B1=cellstr('Relative Abundance of GlycanSpecies');
C1=cellstr('Average # of GlycanStruct per Glycan');
StructName = cell(11,1);
SpFraction = cell(10,1);
StFraction = cell(10,1);
StructName{1} = 'GlycanStruct Name';
StructName{2} = 'LewisX_a';
StructName{3} = 'sialylLewisX_a';
StructName{4} = 'LewisY_b';
StructName{5} = 'LacNAc';
StructName{6} = 'Sialyl_LacNAc';
StructName{7} = 'H_Antigen';
StructName{8} = 'A_Antigen';
StructName{9} = 'B_Antigen';
StructName{10} = 'Bisecting';
StructName{11} = 'coreFuc';

SpFraction{1} = Lx_a.glycanSpecies;
SpFraction{2} = sialylLx_a.glycanSpecies;
SpFraction{3} = Ly_b.glycanSpecies;
SpFraction{4} = LacNAc.glycanSpecies;
SpFraction{5} = Sia_LacNAc.glycanSpecies;
SpFraction{6} = Hantigen.glycanSpecies;
SpFraction{7} = Aantigen.glycanSpecies;
SpFraction{8} = Bantigen.glycanSpecies;
SpFraction{9} = Bisecting.glycanSpecies;
SpFraction{10} = coreFuc.glycanSpecies;

StFraction{1} = Lx_a.glycanStruct;
StFraction{2} = sialylLx_a.glycanStruct;
StFraction{3} = Ly_b.glycanStruct;
StFraction{4} = LacNAc.glycanStruct;
StFraction{5} = Sia_LacNAc.glycanStruct;
StFraction{6} = Hantigen.glycanStruct;
StFraction{7} = Aantigen.glycanStruct;
StFraction{8} = Bantigen.glycanStruct;
StFraction{9} = Bisecting.glycanStruct;
StFraction{10} = coreFuc.glycanStruct;
Excel = evalin('caller','Excel');
xlswrite1(filespec_user,StructName,3,'A1');
xlswrite1(filespec_user,B1,3,'B1'); 
xlswrite1(filespec_user,C1,3,'C1'); 
xlswrite1(filespec_user,SpFraction,3,'B2');
xlswrite1(filespec_user,StFraction,3,'C2');
end