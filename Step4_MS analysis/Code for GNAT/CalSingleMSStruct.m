function CalSingleMSStruct(glycanDB,outputfilename,loadpath)
% CalSingleMSStruct: Calculate relative abundance of different substructure
%
%
%
% Author: Yusen Zhou
% Date Lastly Updated: 05/18/2020


% Load glycanDB
DBfilepath = [loadpath glycanDB '.mat'];
load(DBfilepath)
% calculate MS subStructure relativeabundance(lewis x/a, sialyl-Lewis x/a, 
%              Lewis y/b, H_Antigen, A_Antigen, B_Antigen, Bisecting, coreFuc)
validNum  = 0;
TerStruct ={'LeX','Lea','sialylLeX','sialylLea','Ley','Leb','T2LacNAc','Sia_T2LacNAc','Hantigen',...
    'Aantigen','Bantigen','Bisecting','coreFuc'};
specificstructabundance = cell(length(TerStruct),1);
for i = 1 : length(TerStruct)
    [specificstructabundance{i},validNum] = calbystruct(newglycanDB,TerStruct{i},validNum);
end
SubStruct = cell(validNum,2);
Index     = 1;
for i =1 : length(specificstructabundance)
    if(specificstructabundance{i}.glycanSpecies~=0)
        SubStruct{Index,1} = TerStruct{i};
        SubStruct{Index,2} = specificstructabundance{i};
        Index = Index+1;
    end
end
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
if exist(filespec_user,'file')
    Workbook = invoke(Workbooks,'Open',filespec_user);
else
    Workbook = invoke(Workbooks, 'Add');
    Workbook.SaveAs(filespec_user);
end
write2Excel(SubStruct,filespec_user);
% plot the bar figure
plotbarfigure_Residue(SubStruct,jpgnameSp,jpgnameSt);
Workbook.Save;
invoke(Excel,'Quit');
Excel.delete 
clear Excel
end

function [specificstructabundance,validNum] = calbystruct(newglycanDB,glycanStruct,validNum)
specificstructabundance = struct('glycanSpecies',[],'glycanStruct',[]);
glycanarrays            = newglycanDB.glycanexpec;
abundance               = newglycanDB.abundance;
SpAbundance = 0;
StAbundance = 0;
for i = 1 : length(glycanarrays)
    ithglyinfo = glycanarrays{i,1};
    ismatch = 0;
    if(length(ithglyinfo(:,1))==1)
        ithglycan = ithglyinfo{1,1};
        if(isempty(ithglycan))
            continue;
        end
        [ismatch,NumberOfStruct] = calbyDiffStruct(ithglycan,glycanStruct,ismatch);
        if(ismatch~=0)
            SpAbundance = SpAbundance+abundance(i);
            StAbundance = StAbundance+abundance(i)*NumberOfStruct;
        end
    elseif(length(ithglyinfo(:,1))>=1)
        for j = 1 : length(ithglyinfo)
            structstr = ithglyinfo{j,1};
            if(isempty(structstr))
                continue;
            end
            NumberOfStruct = calbyDiffStruct(structstr,glycanStruct);
            if(NumberOfStruct>0)
                SpAbundance = SpAbundance + ((abundance(i))/length(ithglyinfo));
                StAbundance = StAbundance + ((abundance(i)*NumberOfStruct)/length(ithglyinfo));
            end
        end
    end
end
if(SpAbundance~=0)
    validNum = validNum+1;
end
SpAbundance = roundn(SpAbundance*100,-1);
StAbundance = roundn(StAbundance*100,-2);
specificstructabundance.glycanSpecies = SpAbundance;
specificstructabundance.glycanStruct  = StAbundance;
end

function NumberOfStruct = calbyDiffStruct(structstr,glycanStruct)
NumberOfStruct = 0;
if(strcmp(glycanStruct,'LeX'))
    StdStruct = '{(b1-?)GlcNAc[{(a1-3)Fuc}]{(b1-4)Gal}}';
    StdStruct2 = '{(b1-?)GlcNAc[{(a1-3)Fuc}]{(b1-4)Gal{(a2-3)NeuAc}}}';
    NumberOfStruct1 = NumofBinA(structstr,StdStruct,MonoIdx,0);
    NumberOfStruct2 = NumofBinA(structstr,StdStruct2,MonoIdx,0);
    NumberOfStruct  = NumberOfStruct1-NumberOfStruct2;
elseif(strcmp(glycanStruct,'Lea'))
    StdStruct = '{(b1-?)GlcNAc[{(a1-3)Fuc}]{(b1-3)Gal}}';
    StdStruct2 = '{(b1-?)GlcNAc[{(a1-3)Fuc}]{(b1-3)Gal{(a2-3)NeuAc}}}';
    NumberOfStruct1 = NumofBinA(structstr,StdStruct,MonoIdx,0);
    NumberOfStruct2 = NumofBinA(structstr,StdStruct2,MonoIdx,0);
    NumberOfStruct  = NumberOfStruct1-NumberOfStruct2;
elseif(strcmp(glycanStruct,'sialylLeX'))
    StdStruct = '{(b1-?)GlcNAc[{(a1-3)Fuc}]{(b1-4)Gal{(a2-3)NeuAc}}}';
    NumberOfStruct = NumofBinA(structstr,StdStruct,MonoIdx,0);
elseif(strcmp(glycanStruct,'sialylLea'))
    StdStruct = '{(b1-?)GlcNAc[{(a1-3)Fuc}]{(b1-3)Gal{(a2-3)NeuAc}}}';
    NumberOfStruct = NumofBinA(structstr,StdStruct,MonoIdx,0);
elseif(strcmp(glycanStruct,'Ley'))
    StdStruct = '{(b1-?)GlcNAc[{(a1-3)Fuc}]{(b1-4)Gal{(a1-2)Fuc}}}';
    NumberOfStruct = NumofBinA(structstr,StdStruct,MonoIdx,0);
elseif(strcmp(glycanStruct,'Leb'))
    StdStruct = '{(b1-?)GlcNAc[{(a1-3)Fuc}]{(b1-3)Gal{(a1-2)Fuc}}}';
    NumberOfStruct = NumofBinA(structstr,StdStruct,MonoIdx,0);
elseif(strcmp(glycanStruct,'LacNAc'))
    StdStruct = '{(b1-3)GlcNAc{(b1-4)Gal}}';
    NumberOfStruct = NumofBinA(structstr,StdStruct,MonoIdx,0);
elseif(strcmp(glycanStruct,'Sia_LacNAc'))
    StdStruct = '{(b1-?)GlcNAc{(b1-4)Gal{(a2-3)NeuAc}}}';
    NumberOfStruct = NumofBinA(structstr,StdStruct,MonoIdx,0);
elseif(strcmp(glycanStruct,'Hantigen'))
    StdStruct1 = '{(b1-?)GlcNAc{(b1-4)Gal{(a1-2)Fuc}}}';
    StdStruct2 = '{(b1-?)GlcNAc[{(a1-3)Fuc}]{(b1-4)Gal{(a1-2)Fuc}}}';
    NumberOfStruct1 = NumofBinA(structstr,StdStruct1,MonoIdx,0);
    NumberOfStruct2 = NumofBinA(structstr,StdStruct2,MonoIdx,0);
    NumberOfStruct  = NumberOfStruct1-NumberOfStruct2;
elseif(strcmp(glycanStruct,'Aantigen'))
    StdStruct = '{(b1-?)GlcNAc{(b1-4)Gal[{(a1-2)Fuc}]{(a1-3)GalNAc}}}';
    NumberOfStruct = NumofBinA(structstr,StdStruct,MonoIdx,0);   
elseif(strcmp(glycanStruct,'Bantigen'))
    StdStruct = '{(b1-?)GlcNAc{(b1-4)Gal[{(a1-2)Fuc}]{(a1-3)Gal}}}';
    NumberOfStruct = NumofBinA(structstr,StdStruct,MonoIdx,0);  
elseif(strcmp(glycanStruct,'Bisecting'))
    StdStruct = '{(b1-4)Man{(b1-4)GlcNAc}}';
    NumberOfStruct = NumofBinA(structstr,StdStruct,MonoIdx,0);  
elseif(strcmp(glycanStruct,'coreFuc'))
    StdStruct = '{(b1-r)GlcNAc{(a1-6)Fuc}}';
    NumberOfStruct = NumofBinA(structstr,StdStruct,MonoIdx,0);  
end
end

function plotbarfigure_Residue(SubStruct,jpgnameSp,jpgnameSt)

Y1 = zeros(length(SubStruct),1);
for i = 1 : length(SubStruct)
    Y1(i) = SubStruct{i,2}.glycanSpecies;
end
figure
bar(Y1);
set(gcf,'PaperPositionMode','auto','visible','off','outerposition',[0,0,1800,1500],'position', [0,0,1500,1200]);
startX = 0.8;
Ylabel = cell(1,length(SubStruct)+1);
Ylabel{1,1} = '0';
for i = 1 : length(SubStruct)
    text(startX,Y1(i)*0.9,num2str(Y1(i)),'Fontsize',15,'color','w')
    startX = startX+1;
    Ylabel{1,i+1} = SubStruct{i,1};
end
set(gca,'xtick',0:1:length(Y1));
set(gca,'xticklabel',Ylabel,'Fontsize',15);
leg = legend('% relative abundance of Glycan (structure)');
set(leg,'Fontsize',20);
print(gcf,'-djpeg',jpgnameSp);

Y2 = zeros(length(SubStruct),1);
for i = 1 : length(SubStruct)
    Y2(i) = SubStruct{i,2}.glycanStruct;
end
figure
bar(Y2);
set(gcf,'PaperPositionMode','auto','visible','off','outerposition',[0,0,1800,1500],'position', [0,0,1500,1200]);
startX = 0.8;
for i = 1 : length(SubStruct)
    text(startX,Y2(i)*0.9,num2str(Y2(i)),'Fontsize',15,'color','w')
    startX = startX+1;
end
set(gca,'xtick',0:1:length(Y2));
set(gca,'xticklabel',Ylabel,'Fontsize',15);
leg = legend('% Average # of GlycanStruct per Glycan');
set(leg,'Fontsize',20);
print(gcf,'-djpeg',jpgnameSt);

Excel = evalin('caller','Excel');
Sheets = Excel.ActiveWorkBook.Sheets;
sheet1 = get(Sheets, 'Item', 3);
invoke(sheet1, 'Activate')
Sheets.Item(3).Name = 'sub-structure';
Shapes = sheet1.Shapes;
Shapes.AddPicture(jpgnameSp ,0,1,40,200,1400,700);
Shapes.AddPicture(jpgnameSt ,0,1,40,950,1400,700);
delete(gcf);
end

function write2Excel(SubStruct,filespec_user)
B1=cellstr('Relative Abundance of GlycanSpecies');
C1=cellstr('Average # of GlycanStruct per Glycan');
StructName = cell(length(SubStruct)+1,1);
SpFraction = cell(length(SubStruct),1);
StFraction = cell(length(SubStruct),1);
StructName{1} = 'GlycanStruct Name';
for i = 1 : length(SubStruct)
    StructName{i+1} = SubStruct{i,1};
    SpFraction{i} = SubStruct{i,2}.glycanSpecies;
    StFraction{i} = SubStruct{i,2}.glycanStruct;
end

Excel = evalin('caller','Excel');
xlswrite1(filespec_user,StructName,3,'A1');
xlswrite1(filespec_user,B1,3,'B1'); 
xlswrite1(filespec_user,C1,3,'C1'); 
xlswrite1(filespec_user,SpFraction,3,'B2');
xlswrite1(filespec_user,StFraction,3,'C2');
end

function [ismatch,NumberOfStruct] = CheckIfMatch(structstr,strPat,target,ismatch)
NumberOfStruct = 0;
strMatch = regexp(structstr,strPat,'Match');
Lvl      = 0;
for i = 1 : length(strMatch)
    if(strcmp(strMatch{i},'['))
        Lvl = Lvl+1;
    elseif(strcmp(strMatch{i},']'))
        Lvl = Lvl-1;
    elseif(strcmp(strMatch{i},target))
        if(Lvl==0)
            ismatch = ismatch+1;
            NumberOfStruct = 1;
        elseif(Lvl==1)
            if(strcmp(strMatch{i-1},'['))
               ismatch = ismatch+1;
               NumberOfStruct = 1; 
            end
        end
    end
end
end