function StructResult = MSStructAnalysisGNAT(MSfilename,newglycanDB)
% CalSingleMSStruct: Calculate relative abundance of different substructure for multi
% glycan profiles.
%
%
%
% Author: Yusen Zhou
% Date Lastly Updated: 05/18/2020

StructResult = struct('MSID',[],'Lewisx',[],'Lewisa',[],'sialylLewisx',[],'sialylLewisa',[],'Lewisy',[],'Lewisb',[],...
    'LacNAc',[],'SialylLacNAc',[],'Hantigen',[],'Aantigen',[],'Bantigen',[],'Bisecting',[],'coreFuc',[]);
% calculate MS subStructure relativeabundance(lewis x/a, sialyl-Lewis x/a, 
%              Lewis y/b, H_Antigen, A_Antigen, B_Antigen, Bisecting, coreFuc)
MonoIdx      = MonoIDgeneration;
LeX          = calbystruct(newglycanDB,'LeX',MonoIdx);
Lea          = calbystruct(newglycanDB,'Lea',MonoIdx);
sialylLeX    = calbystruct(newglycanDB,'sialylLeX',MonoIdx);
sialylLea    = calbystruct(newglycanDB,'sialylLea',MonoIdx);
Ley          = calbystruct(newglycanDB,'Ley',MonoIdx);
Leb          = calbystruct(newglycanDB,'Leb',MonoIdx);
LacNAc       = calbystruct(newglycanDB,'LacNAc',MonoIdx);
Sia_LacNAc   = calbystruct(newglycanDB,'Sia_LacNAc',MonoIdx);
Hantigen     = calbystruct(newglycanDB,'Hantigen',MonoIdx);
Aantigen     = calbystruct(newglycanDB,'Aantigen',MonoIdx);
Bantigen     = calbystruct(newglycanDB,'Bantigen',MonoIdx);
Bisecting    = calbystruct(newglycanDB,'Bisecting',MonoIdx);
coreFuc      = calbystruct(newglycanDB,'coreFuc',MonoIdx);

StructResult.MSID             = MSfilename;
StructResult.Lewisx           = LeX;
StructResult.Lewisa           = Lea;
StructResult.sialylLewisx     = sialylLeX;
StructResult.sialylLewisa     = sialylLea;
StructResult.Lewisy           = Ley;
StructResult.Lewisb           = Leb;
StructResult.LacNAc           = LacNAc;
StructResult.SialylLacNAc     = Sia_LacNAc;
StructResult.Hantigen         = Hantigen;
StructResult.Aantigen         = Aantigen;
StructResult.Bantigen         = Bantigen;
StructResult.Bisecting        = Bisecting;
StructResult.coreFuc          = coreFuc;
end

function specificstructabundance = calbystruct(newglycanDB,glycanStruct,MonoIdx)
specificstructabundance = struct('glycanSpecies',[],'glycanStruct',[]);
glycanarrays            = newglycanDB.glycanexpec;
abundance               = newglycanDB.abundance;
SpAbundance = 0;
StAbundance = 0;
for i = 1 : length(glycanarrays)
    ithglyinfo = glycanarrays{i,1};
    if(length(ithglyinfo(:,1))==1)
        ithglycan = ithglyinfo{1,1};
        ithglycan = regexprep(ithglycan,'\(b1-\)','\(b1-r\)');
        NumberOfStruct = calbyDiffStruct(ithglycan,glycanStruct,MonoIdx);
        if(NumberOfStruct>0)
            SpAbundance = SpAbundance+abundance(i);
            StAbundance = StAbundance+abundance(i)*NumberOfStruct;
        end
     elseif(length(ithglyinfo(:,1))>=1)
        for j = 1 : length(ithglyinfo)
            structstr = ithglyinfo{j,1};
            structstr = regexprep(structstr,'\(b1-\)','\(b1-r\)');
            NumberOfStruct = calbyDiffStruct(structstr,glycanStruct,MonoIdx);
            if(NumberOfStruct>0)
                SpAbundance = SpAbundance + ((abundance(i))/length(ithglyinfo));
                StAbundance = StAbundance + ((abundance(i)*NumberOfStruct)/length(ithglyinfo));
            end
        end
    end
end
SpAbundance = roundn(SpAbundance*100,-1);
StAbundance = roundn(StAbundance*100,-1);
specificstructabundance.glycanSpecies = SpAbundance;
specificstructabundance.glycanStruct  = StAbundance;
end

function NumberOfStruct = calbyDiffStruct(structstr,glycanStruct,MonoIdx)
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