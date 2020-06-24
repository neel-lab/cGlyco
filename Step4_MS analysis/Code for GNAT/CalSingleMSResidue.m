function CalSingleMSResidue(glycanDB,outputfilename,loadpath)
% CalSingleMSResidue: Calculate relative abundance of different
% monosaccharides
%
%
%
% Author: Yusen Zhou
% Date Lastly Updated: 05/18/2020

MonoDBPat =  {'Glc[^AN]','Man[^AN]','Gal[^AN]','Gul[^AN]','Alt[^AN]','All[^AN]',...
    'Tal[^AN]','Ido[^AN]','GlcNAc','ManNAc','GalNAc','GulNAc','AltNAc','AllNAc',...
    'TalNAc','IdoNAc','GlcN[^A]','ManN[^A]','GalN[^A]','GulN[^A]','AltN[^A]',...
    'AllN[^A]','TalN[^A]','IdoN[^A]','GlcA','ManA','GalA','GulA','AltA',...
    'AllA','TalA','IdoA','QuiNAc','RhaNAc','FucNAc','Qui[^N]','Rha[^N]',...
    '6dAlt','6dTal','Fuc[^N]','Oli','Tyv','Abe','Par','Dig','Col','Ara',...
    'Lyx','Xyl','Rib','Kdn','NeuGc','NeuAc','Neu[^A]','Bac','LDMan','Kdo',...
    'Dha','DDMan','MurNAc','MurNGc','Mur','Api','Fru','Tag','Sor','Psi'};
MonoDB =  {'Glc','Man','Gal','Gul','Alt','All',...
    'Tal','Ido','GlcNAc','ManNAc','GalNAc','GulNAc','AltNAc','AllNAc',...
    'TalNAc','IdoNAc','GlcN','ManN','GalN','GulN','AltN',...
    'AllN','TalN','IdoN','GlcA','ManA','GalA','GulA','AltA',...
    'AllA','TalA','IdoA','QuiNAc','RhaNAc','FucNAc','Qui','Rha',...
    '6dAlt','6dTal','Fuc','Oli','Tyv','Abe','Par','Dig','Col','Ara',...
    'Lyx','Xyl','Rib','Kdn','NeuGc','NeuAc','Neu','Bac','LDMan','Kdo',...
    'Dha','DDMan','MurNAc','MurNGc','Mur','Api','Fru','Tag','Sor','Psi'};
% Load the glycanlist and glycanDB
DBfilepath = [loadpath glycanDB '.mat'];
load(DBfilepath)
% calculate MS Residue relativeabundance(Gal, GlcNAc, GalNAc, Fuc, NeuAc)
specificresidueabundance = cell(length(MonoDBPat),1);
validNum = 0;
for i = 1 : length(MonoDBPat)
    [specificresidueabundance{i},validNum] = calculateresidue(newglycanDB,MonoDBPat{i},validNum);
end
ResidueResult = cell(length(validNum),2);
Index         = 1;
for i =1 : length(specificresidueabundance)
    if(specificresidueabundance{i}.glycanAbundance~=0)
        ResidueResult{Index,1} = MonoDB{i};
        ResidueResult{Index,2} = specificresidueabundance{i};
        Index = Index+1;
    end
end
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
write2Excel(ResidueResult,filespec_user);
% plot the bar figure
plotbarfigure_Residue(ResidueResult,jpgnameS,jpgnameR);
Workbook.Save;
invoke(Excel,'Quit');
Excel.delete 
clear Excel
end

function [specificresidueabundance,validNum] = calculateresidue(newglycanDB,GlycanResidueName,validNum)
specificresidueabundance = struct('glycanAbundance',[],'residueAbundance',[]);
[glycanAbundance,residueAbundance] = addrequiredresidues(newglycanDB,GlycanResidueName);
if(glycanAbundance~=0)
    validNum = validNum+1;
end
specificresidueabundance.glycanAbundance  = glycanAbundance;
specificresidueabundance.residueAbundance = residueAbundance;
end

function [glycanAbundance,residueAbundance] = addrequiredresidues(newglycanDB,GlycanResidueName)
glycanAbundance   = 0;
residueAbundance  = 0;
glycanarrays = newglycanDB.glycanexpec;
for i = 1 : length(glycanarrays)
    ithglyinfo = glycanarrays{i,1};
    if(length(ithglyinfo(:,1))==1)
        ithglycan = ithglyinfo{1,1};
        if(~isempty(regexp(ithglycan,GlycanResidueName,'once')))
            numofresidue     = length(regexp(ithglycan,GlycanResidueName));
            glycanAbundance  = glycanAbundance+newglycanDB.abundance(i);
            residueAbundance = residueAbundance+newglycanDB.abundance(i)*numofresidue;
        end
    elseif(length(ithglyinfo(:,1))>=1)
        glycanstruct1 = ithglyinfo{1,1};
        if(~isempty(regexp(glycanstruct1,GlycanResidueName,'Once')))
            numofresidue     = length(regexp(glycanstruct1,GlycanResidueName));
            glycanAbundance  = glycanAbundance+newglycanDB.abundance(i);
            residueAbundance = residueAbundance+newglycanDB.abundance(i)*numofresidue;
        end
    end
end
glycanAbundance  = roundn(glycanAbundance*100,-1);
residueAbundance = roundn(residueAbundance,-1);
end

function plotbarfigure_Residue(ResidueResult,jpgnameS,jpgnameR)
Y1 = zeros(length(ResidueResult),1);
for i = 1 : length(ResidueResult)
    Y1(i) = ResidueResult{i,2}.glycanAbundance;
end
figure
bar(Y1);
set(gcf,'PaperPositionMode','auto','visible','off','outerposition',[0,0,1800,1500],'position', [0,0,1500,1200]);
startX = 0.8;
Xlabel = cell(1,length(ResidueResult)+1);
Xlabel{1,1} = '0';
for i = 1 : length(ResidueResult)
    text(startX,Y1(i)*0.85,num2str(Y1(i)),'Fontsize',18,'color','w')
    startX = startX+1;
    Xlabel{1,i+1} = ResidueResult{i,1};
end
set(gca,'xtick',0:1:length(ResidueResult));
set(gca,'xticklabel',Xlabel,'Fontsize',20);
leg = legend('% relative abundance (glycan)');
set(leg,'Fontsize',20);
print(gcf,'-djpeg',jpgnameS);
Y2 = zeros(length(ResidueResult),1);
for i = 1 : length(ResidueResult)
    Y2(i) = ResidueResult{i,2}.residueAbundance;
end
figure
bar(Y2);
set(gcf,'PaperPositionMode','auto','visible','off','outerposition',[0,0,1800,1500],'position', [0,0,1500,1200]);
startX = 0.8;
for i = 1 : length(ResidueResult)
    text(startX,Y2(i)*0.85,num2str(Y2(i)),'Fontsize',18,'color','w')
    startX = startX+1;
end
set(gca,'xtick',0:1:length(ResidueResult));
set(gca,'xticklabel',Xlabel,'Fontsize',20);
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
% sheet1.invoke('Pictures').Insert(jpgname);
delete(gcf);
end

function write2Excel(ResidueResult,filespec_user)
ReisdueName = cell(length(ResidueResult)+1,1);
ReisdueFraction1 = cell(length(ResidueResult),1);
ReisdueFraction2 = cell(length(ResidueResult),1);
B1=cellstr('Relative Abundance of GlycanSpeices');
C1=cellstr('Relative Abundance of GlycanResidue');
ReisdueName{1} = 'GlycanResidue Name';
for i = 1 : length(ResidueResult)
    ReisdueName{i+1}    = ResidueResult{i,1};
    ReisdueFraction1{i} = ResidueResult{i,2}.glycanAbundance;
    ReisdueFraction2{i} = ResidueResult{i,2}.residueAbundance;
end

Excel = evalin('caller','Excel');
xlswrite1(filespec_user,ReisdueName,2,'A1');
xlswrite1(filespec_user,B1,2,'B1');
xlswrite1(filespec_user,C1,2,'C1');
xlswrite1(filespec_user,ReisdueFraction1,2,'B2');
xlswrite1(filespec_user,ReisdueFraction2,2,'C2');
end