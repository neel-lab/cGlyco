function ResidueResult = MSResidueAnalysis(MSfilename,newglycanDB,withHighmannose,storepath)
% MSResidueAnalysis: Calculate relative abundance of different
% monosaccharides for multi glycan profiles
%
%
%
% Author: Yusen Zhou
% Date Lastly Updated: 05/18/2020
ResidueResult = struct('MSID',[],'Gal',[],'GlcNAc',[],'GalNAc',[],'Fuc',[],'NeuAc',[]);
% calculate MS Residue relativeabundance(Gal, GlcNAc, GalNAc, Fuc, NeuAc)
GalAbundance    = calculatebyStr(newglycanDB,'Gal',withHighmannose);
GlcNAcAbundance = calculatebyStr(newglycanDB,'GlcNAc',withHighmannose);
GalNAcAbundance = calculatebyStr(newglycanDB,'GalNAc',withHighmannose);
FucAbundance    = calculatebyStr(newglycanDB,'Fuc',withHighmannose);
NeuAcAbundance  = calculatebyStr(newglycanDB,'NeuAc',withHighmannose);

ResidueResult.MSID       = MSfilename;
ResidueResult.Gal        = GalAbundance;
ResidueResult.GlcNAc     = GlcNAcAbundance;
ResidueResult.GalNAc     = GalNAcAbundance;
ResidueResult.Fuc        = FucAbundance;
ResidueResult.NeuAc      = NeuAcAbundance;

Storepath = [storepath MSfilename 'Residue.mat'];
save(Storepath,'ResidueResult');
end

function ResidueAbundance  = calculatebyStr(newglycanDB,ResidueStr,withHighmannose)
ResidueAbundance = struct('glycanSpecies',[],'glycanResidue',[]);
glycanStructStr = newglycanDB.glycanexpec;
abundance       = newglycanDB.abundance;
SAbundance = 0;
RAbundance = 0;
if(strcmp(ResidueStr,'Gal'))
    for i = 1 : length(glycanStructStr)
        if(ischar(glycanStructStr{i}))
            ithstr = glycanStructStr{i};
        elseif(isa(glycanStructStr{i},'cell'))
            ithstr = glycanStructStr{i}{1,1};
        end
        if(~withHighmannose)
            if(length(strfind(ithstr,'h'))>3)
                numofGal   = length(strfind(ithstr,'h'))-3;
                SAbundance = SAbundance+abundance(i);
                RAbundance = RAbundance+abundance(i)*numofGal;
            end
        elseif(withHighmannose)
            isHighmannose = chkhighmannose(ithstr);
            if(~isHighmannose)
                if(length(strfind(ithstr,'h'))>3)
                    numofGal   = length(strfind(ithstr,'h'))-3;
                    SAbundance = SAbundance+abundance(i);
                    RAbundance = RAbundance+abundance(i)*numofGal;
                end
            end
        end
    end
elseif(strcmp(ResidueStr,'GlcNAc'))
    for i = 1 : length(glycanStructStr)
        if(ischar(glycanStructStr{i}))
            ithstr = glycanStructStr{i};
        elseif(isa(glycanStructStr{i},'cell'))
            ithstr = glycanStructStr{i}{1,1};
        end
        numofHexNAc = length(strfind(ithstr,'n'));
        numofGalNAc = length(strfind(ithstr,'{h_b4{f_a2}{n_a3}'))+...
            length(strfind(ithstr,'{h{f}{n}'));
        numofGlcNAc = numofHexNAc-numofGalNAc;
        if(numofHexNAc>0)
            SAbundance = SAbundance+abundance(i);
            RAbundance = RAbundance+abundance(i)*numofGlcNAc;
        end
    end
elseif(strcmp(ResidueStr,'GalNAc'))
    for i = 1 : length(glycanStructStr)
        if(ischar(glycanStructStr{i}))
            ithstr = glycanStructStr{i};
        elseif(isa(glycanStructStr{i},'cell'))
            ithstr = glycanStructStr{i}{1,1};
        end
        numofGalNAc = length(strfind(ithstr,'{h_b4{f_a2}{n_a3}'))+...
            length(strfind(ithstr,'{h{f}{n}'));
        if(numofGalNAc>0)
            SAbundance = SAbundance+abundance(i);
            RAbundance = RAbundance+abundance(i)*numofGalNAc;
        end
    end
elseif(strcmp(ResidueStr,'Fuc'))
    for i = 1 : length(glycanStructStr)
        if(ischar(glycanStructStr{i}))
            ithstr = glycanStructStr{i};
        elseif(isa(glycanStructStr{i},'cell'))
            ithstr = glycanStructStr{i}{1,1};
        end
        numofFuc = length(strfind(ithstr,'f'));
        if(numofFuc>0)
            SAbundance = SAbundance+abundance(i);
            RAbundance = RAbundance+abundance(i)*numofFuc;
        end
    end
elseif(strcmp(ResidueStr,'NeuAc'))
    for i = 1 : length(glycanStructStr)
        if(ischar(glycanStructStr{i}))
            ithstr = glycanStructStr{i};
        elseif(isa(glycanStructStr{i},'cell'))
            ithstr = glycanStructStr{i}{1,1};
        end
        numofNeuAc = length(strfind(ithstr,'s'));
        if(numofNeuAc>0)
            SAbundance = SAbundance+abundance(i);
            RAbundance = RAbundance+abundance(i)*numofNeuAc;
        end
    end
end
ResidueAbundance.glycanSpecies = roundn(SAbundance*100,-1);
ResidueAbundance.glycanResidue = roundn(RAbundance*100,-1);
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