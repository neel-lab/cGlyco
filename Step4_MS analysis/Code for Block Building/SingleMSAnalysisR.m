function SingleMSAnalysisR(GlycanDB,withHighmannose,outputfilename,loadpath)
% SingleMSAnalysisR: Calculate relative abundance of different
% monosaccharides 
%
%
%
% Author: Yusen Zhou
% Date Lastly Updated: 05/18/2020
% Load the glycanDB
DBfilepath = [loadpath GlycanDB '.mat'];
load(DBfilepath)
% calculate MS Residue relativeabundance(Gal, GlcNAc, GalNAc, Fuc, NeuAc)
GalAbundance    = calculatebyStr(newglycanDB,'Gal',withHighmannose);
GlcNAcAbundance = calculatebyStr(newglycanDB,'GlcNAc',withHighmannose);
GalNAcAbundance = calculatebyStr(newglycanDB,'GalNAc',withHighmannose);
FucAbundance    = calculatebyStr(newglycanDB,'Fuc',withHighmannose);
NeuAcAbundance  = calculatebyStr(newglycanDB,'NeuAc',withHighmannose);
% Export to the Excel
filespec_user = [loadpath outputfilename '.xlsx'];
jpgnameS = [loadpath outputfilename 'ResidueAbundanceS.jpg'];
jpgnameR = [loadpath outputfilename 'ResidueAbundanceR.jpg'];
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
write2Excel(GalAbundance,GlcNAcAbundance,GalNAcAbundance,FucAbundance,NeuAcAbundance,filespec_user);
% plot the bar figure
plotbarfigure_Residue(GalAbundance,GlcNAcAbundance,GalNAcAbundance,FucAbundance,NeuAcAbundance,...
    jpgnameS,jpgnameR);
Workbook.Save;
invoke(Excel,'Quit');
Excel.delete 
clear Excel
end

function ResidueAbundance  = calculatebyStr(newglycanDB,ResidueStr,withHighmannose)
ResidueAbundance = struct('glycanSpecies',[],'glycanResidue',[]);
glycanStructStr = newglycanDB.glycanexpec;
SAbundance = 0;
RAbundance = 0;
if(strcmp(ResidueStr,'Gal'))
    for i = 1 : length(glycanStructStr)
        if(ischar(glycanStructStr{i}))
            ithstr = glycanStructStr{i};
        elseif(isa(glycanStructStr{i},'cell'))
            ithstr = glycanStructStr{i}{1};
        end
        if(~withHighmannose)
            if(length(strfind(ithstr,'h'))>3)
                numofGal   = length(strfind(ithstr,'h'))-3;
                SAbundance = SAbundance+newglycanDB.abundance(i);
                RAbundance = RAbundance+newglycanDB.abundance(i)*numofGal;
            end
        elseif(withHighmannose)
            isHighmannose = chkhighmannose(ithstr);
            if(~isHighmannose)
                if(length(strfind(ithstr,'h'))>3)
                    numofGal   = length(strfind(ithstr,'h'))-3;
                    SAbundance = SAbundance+newglycanDB.abundance(i);
                    RAbundance = RAbundance+newglycanDB.abundance(i)*numofGal;
                end
            end
        end
    end
elseif(strcmp(ResidueStr,'GlcNAc'))
    count = 0;
    for i = 1 : length(glycanStructStr)
        if(ischar(glycanStructStr{i}))
            ithstr = glycanStructStr{i};
        elseif(isa(glycanStructStr{i},'cell'))
            ithstr = glycanStructStr{i}{1};
        end
        numofHexNAc = length(strfind(ithstr,'n'));
        numofGalNAc = length(strfind(ithstr,'{h_b4{f_a2}{n_a3}'))+...
            length(strfind(ithstr,'{h{f}{n}'));
        numofGlcNAc = numofHexNAc-numofGalNAc;
        if(numofHexNAc>0)
            SAbundance = SAbundance+newglycanDB.abundance(i);
            RAbundance = RAbundance+newglycanDB.abundance(i)*numofGlcNAc;
            count = count+1;
        end
    end
elseif(strcmp(ResidueStr,'GalNAc'))
    for i = 1 : length(glycanStructStr)
        if(ischar(glycanStructStr{i}))
            ithstr = glycanStructStr{i};
        elseif(isa(glycanStructStr{i},'cell'))
            ithstr = glycanStructStr{i}{1};
        end
        numofGalNAc = length(strfind(ithstr,'{h_b4{f_a2}{n_a3}'))+...
            length(strfind(ithstr,'{h{f}{n}'));
        if(numofGalNAc>0)
            SAbundance = SAbundance+newglycanDB.abundance(i);
            RAbundance = RAbundance+newglycanDB.abundance(i)*numofGalNAc;
        end
    end
elseif(strcmp(ResidueStr,'Fuc'))
    for i = 1 : length(glycanStructStr)
        if(ischar(glycanStructStr{i}))
            ithstr = glycanStructStr{i};
        elseif(isa(glycanStructStr{i},'cell'))
            ithstr = glycanStructStr{i}{1};
        end
        numofFuc = length(strfind(ithstr,'f'));
        if(numofFuc>0)
            SAbundance = SAbundance+newglycanDB.abundance(i);
            RAbundance = RAbundance+newglycanDB.abundance(i)*numofFuc;
        end
    end
elseif(strcmp(ResidueStr,'NeuAc'))
    for i = 1 : length(glycanStructStr)
        if(ischar(glycanStructStr{i}))
            ithstr = glycanStructStr{i};
        elseif(isa(glycanStructStr{i},'cell'))
            ithstr = glycanStructStr{i}{1};
        end
        numofNeuAc = length(strfind(ithstr,'s'));
        if(numofNeuAc>0)
            SAbundance = SAbundance+newglycanDB.abundance(i);
            RAbundance = RAbundance+newglycanDB.abundance(i)*numofNeuAc;
        end
    end
end
ResidueAbundance.glycanSpecies = roundn(SAbundance*100,-1);
ResidueAbundance.glycanResidue = roundn(RAbundance,-1);
end

function isHighmannose = chkhighmannose(Str)
isHighmannose = 0;
highmannoselib = {'{n_b{n_b4{h_b4{h_a3{n_b2}}{h_a6{h_a3}}}}}';...
    '{n_b{n_b4{h_b4{h_a3}{h_a6{h_a3}{h_a6}}}}}';...
    '{n_b{n_b4{h_b4{h_a3{n_b2}}{h_a6{h_a3}{h_a6}}}}}';...
    '{n_b{n_b4{h_b4{h_a3{h_a2}}{h_a6{h_a3}{h_a6}}}}}';...
    '{n_b{n_b4{h_b4{h_a3}{h_a6{h_a3{h_a2}}{h_a6}}}}}';...
    '{n_b{n_b4{h_b4{h_a3}{h_a6{h_a3}{h_a6{h_a2}}}}}}';...
    '{n_b{n_b4{h_b4{h_a3}{h_a6{h_a3{h_a2}}{h_a6{h_a2}}}}}}';...
    '{n_b{n_b4{h_b4{h_a3{h_a2}{h_a2}}{h_a6{h_a3}{h_a6}}}}}';...
    '{n_b{n_b4{h_b4{h_a3{h_a2}}{h_a6{h_a3{h_a2}}{h_a6}}}}}';...
    '{n_b{n_b4{h_b4{h_a3{h_a2}}{h_a6{h_a3}{h_a6{h_a2}}}}}}';...
    '{n_b{n_b4{h_b4{h_a3{h_a2}}{h_a6{h_a3{h_a2}}{h_a6{h_a2}}}}}}';...
    '{n_b{n_b4{h_b4{h_a3{h_a2{h_a2}}}{h_a6{h_a3{h_a2}}{h_a6}}}}}';...
    '{n_b{n_b4{h_b4{h_a3{h_a2{h_a2}}}{h_a6{h_a3}{h_a6{h_a2}}}}}}';...
    '{n_b{n_b4{h_b4{h_a3{h_a2{h_a2}}}{h_a6{h_a3{h_a2}}{h_a6{h_a2}}}}}}';...
    '{n{n{h{h{n}}{h{h}}}}}';...
    '{n{n{h{h}{h{h}{h}}}}}';...
    '{n{n{h{h{n}}{h{h}{h}}}}}';...
    '{n{n{h{h{h}}{h{h}{h}}}}}';...
    '{n{n{h{h{h}}{h{h{h}}{h}}}}}';...
    '{n{n{h{h{h{h}}}{h{h{h}}{h}}}}}';...
    '{n{n{h{h{h{h}}}{h{h{h}}{h{h}}}}}}'};
glycanmwarray = arrayfun(@(x)strcmp(x,Str),highmannoselib);
if(sum(glycanmwarray))
    isHighmannose = 1;
end
end

function plotbarfigure_Residue(GalAbundance,GlcNAcAbundance,GalNAcAbundance,FucAbundance,NeuAcAbundance,...
    jpgnameS,jpgnameR)

Y1 = [GalAbundance.glycanSpecies;...
    GlcNAcAbundance.glycanSpecies;...
    GalNAcAbundance.glycanSpecies;...
    FucAbundance.glycanSpecies;...
    NeuAcAbundance.glycanSpecies];
figure
bar(Y1);
set(gcf,'PaperPositionMode','auto','visible','off','outerposition',[0,0,1800,1500],'position', [0,0,1500,1200]);
text(0.8,GalAbundance.glycanSpecies-0.15*GalAbundance.glycanSpecies,num2str(GalAbundance.glycanSpecies),'Fontsize',18,'color','w');
text(1.8,GlcNAcAbundance.glycanSpecies-0.15*GlcNAcAbundance.glycanSpecies,num2str(GlcNAcAbundance.glycanSpecies),'Fontsize',18,'color','w');
text(2.8,GalNAcAbundance.glycanSpecies-0.15*GalNAcAbundance.glycanSpecies,num2str(GalNAcAbundance.glycanSpecies),'Fontsize',18,'color','w');
text(3.8,FucAbundance.glycanSpecies-0.15*FucAbundance.glycanSpecies,num2str(FucAbundance.glycanSpecies),'Fontsize',18,'color','w');
text(4.8,NeuAcAbundance.glycanSpecies-0.15*NeuAcAbundance.glycanSpecies,num2str(NeuAcAbundance.glycanSpecies),'Fontsize',18,'color','w');
set(gca,'xtick',0:1:5);
set(gca,'xticklabel',{'0' 'Gal' 'GlcNAc' 'GalNAc' 'Fuc' 'NeuAc'},'Fontsize',20);
leg = legend('% relative abundance (glycan)');
set(leg,'Fontsize',20);
print(gcf,'-djpeg',jpgnameS);
Y2 = [GalAbundance.glycanResidue;...
    GlcNAcAbundance.glycanResidue;...
    GalNAcAbundance.glycanResidue;...
    FucAbundance.glycanResidue;...
    NeuAcAbundance.glycanResidue];
figure
bar(Y2);
set(gcf,'PaperPositionMode','auto','visible','off','outerposition',[0,0,1800,1500],'position', [0,0,1500,1200]);
text(0.8,GalAbundance.glycanResidue-0.15*GalAbundance.glycanResidue,num2str(GalAbundance.glycanResidue),'Fontsize',18,'color','w');
text(1.8,GlcNAcAbundance.glycanResidue-0.15*GlcNAcAbundance.glycanResidue,num2str(GlcNAcAbundance.glycanResidue),'Fontsize',18,'color','w');
text(2.8,GalNAcAbundance.glycanResidue-0.15*GalNAcAbundance.glycanResidue-0.1,num2str(GalNAcAbundance.glycanResidue),'Fontsize',18,'color','w');
text(3.8,FucAbundance.glycanResidue-0.15*FucAbundance.glycanResidue,num2str(FucAbundance.glycanResidue),'Fontsize',18,'color','w');
text(4.8,NeuAcAbundance.glycanResidue-0.15*NeuAcAbundance.glycanResidue,num2str(NeuAcAbundance.glycanResidue),'Fontsize',18,'color','w');
set(gca,'xtick',0:1:5);
set(gca,'xticklabel',{'0' 'Gal' 'GlcNAc' 'GalNAc' 'Fuc' 'NeuAc'},'Fontsize',20);
leg = legend('Avg. number of mono./glycan');
set(leg,'Fontsize',20);
print(gcf,'-djpeg',jpgnameR);
Excel = evalin('caller','Excel');
Sheets = Excel.ActiveWorkBook.Sheets;
sheet1 = get(Sheets, 'Item', 2);
invoke(sheet1, 'Activate')
Sheets.Item(2).Name = 'monosaccharide';
Shapes = sheet1.Shapes;
Shapes.AddPicture(jpgnameS ,0,1,40,200,1400,700);
Shapes.AddPicture(jpgnameR ,0,1,40,950,1400,700);
delete(gcf);
end

function write2Excel(GalAbundance,GlcNAcAbundance,GalNAcAbundance,FucAbundance,NeuAcAbundance,filespec_user)
B1=cellstr('Relative Abundance of GlycanSpeices');
C1=cellstr('Relative Abundance of GlycanResidue');
ReisdueName      = cell(6,1);
ReisdueFraction1 = cell(5,1);
ReisdueFraction2 = cell(5,1);

ReisdueName{1} = 'GlycanResidue Name';
ReisdueName{2} = 'Gal';
ReisdueName{3} = 'GlcNAc';
ReisdueName{4} = 'GalNAc';
ReisdueName{5} = 'Fuc';
ReisdueName{6} = 'NeuAc';

ReisdueFraction1{1} = GalAbundance.glycanSpecies;
ReisdueFraction1{2} = GlcNAcAbundance.glycanSpecies;
ReisdueFraction1{3} = GalNAcAbundance.glycanSpecies;
ReisdueFraction1{4} = FucAbundance.glycanSpecies;
ReisdueFraction1{5} = NeuAcAbundance.glycanSpecies;

ReisdueFraction2{1} = GalAbundance.glycanResidue;
ReisdueFraction2{2} = GlcNAcAbundance.glycanResidue;
ReisdueFraction2{3} = GalNAcAbundance.glycanResidue;
ReisdueFraction2{4} = FucAbundance.glycanResidue;
ReisdueFraction2{5} = NeuAcAbundance.glycanResidue;

Excel = evalin('caller','Excel');
xlswrite1(filespec_user,ReisdueName,2,'A1');
xlswrite1(filespec_user,B1,2,'B1');
xlswrite1(filespec_user,C1,2,'C1');
xlswrite1(filespec_user,ReisdueFraction1,2,'B2');
xlswrite1(filespec_user,ReisdueFraction2,2,'C2');
end