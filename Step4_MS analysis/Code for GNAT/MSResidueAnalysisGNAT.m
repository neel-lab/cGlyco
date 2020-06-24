function ResidueResult = MSResidueAnalysisGNAT(MSfilename,newglycanDB)
% CalSingleMSResidue: Calculate relative abundance of different
% monosaccharides for multiple glycan profiles
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
ResidueResult = struct();
% calculate MS Residue relativeabundance
specificresidueabundance = cell(length(MonoDBPat),1);
validNum = 0;
for i = 1 : length(MonoDBPat)
    [specificresidueabundance{i},validNum] = calculateresidue(newglycanDB,MonoDBPat{i},validNum);
end
ResidueResult.MSID = MSfilename;
for i =1 : length(specificresidueabundance)
    if(specificresidueabundance{i}.glycanSpecies~=0)
        ResidueResult.(MonoDB{i}) = specificresidueabundance{i};
    end
end
end

function [specificresidueabundance,validNum] = calculateresidue(newglycanDB,GlycanResidueName,validNum)
specificresidueabundance = struct('glycanSpecies',[],'glycanResidue',[]);
[glycanAbundance,residueAbundance] = addrequiredresidues(newglycanDB,GlycanResidueName);
if(glycanAbundance~=0)
    validNum = validNum+1;
end
specificresidueabundance.glycanSpecies = glycanAbundance;
specificresidueabundance.glycanResidue = residueAbundance;
end

function [glycanAbundance,residueAbundance] = addrequiredresidues(newglycanDB,GlycanResidueName)
glycanAbundance   = 0;
residueAbundance  = 0;
glycanarrays = newglycanDB.glycanexpec;
abundance    = newglycanDB.abundance;
for i = 1 : length(glycanarrays)
    ithglyinfo = glycanarrays{i,1};
    if(length(ithglyinfo(:,1))==1)
        ithglycan = ithglyinfo{1,1};
        if(~isempty(regexp(ithglycan,GlycanResidueName,'once')))
            numofresidue     = length(regexp(ithglycan,GlycanResidueName));
            glycanAbundance  = glycanAbundance+abundance(i);
            residueAbundance = residueAbundance+abundance(i)*numofresidue;
        end
    elseif(length(ithglyinfo(:,1))>=1)
        glycanstruct1 = ithglyinfo{1,1};
        if(~isempty(regexp(glycanstruct1,GlycanResidueName,'once')))
            numofresidue     = length(regexp(glycanstruct1,GlycanResidueName));
            glycanAbundance  = glycanAbundance+abundance(i);
            residueAbundance = residueAbundance+abundance(i)*numofresidue;
        end
    end
end
glycanAbundance  = roundn(glycanAbundance*100,-1);
residueAbundance = roundn(residueAbundance*100,-1);
end