function glycanMW = glycanMolWt(gly,varargin)
%glycanMolWt calculate glycan mass
%
%   glycanMW = glycanMolWt(glycanstring)  takes an input of the glycan 
%       string and outputs its molecular weight. Default options for the 
%       calculation of molecular weight are methylated, Na ionized,
%       monoisotopic and free reducing end.
%
%   glycanMW = glycanMolWt(glycanstring,options)  uses customized options
%       structure for the calcuation of molecular weight. The fields in the
%       options include 'methylation' ([true]/false),'acetylation'([true]/false),'ion' ('none','Na','K','H+','H-'), 
%       'mono' ([true]/false), 'redEnd' ([true]/false).
%   
%   Example1:
%        options.ion='none'; 
%        options.methylation=false;
%        glycanMolWt('shn',options); for 'NeuAcHexHexNAc'
%        ans=674.2381  -- this is the same values as GlycanMass
%        Linkage specificity mentioned in the text is ignored as are the brackets
%
%    Example 2:
%        glycanMolWt('sa2,3hb1,3{s2,6}n',options); for NeuAc2,3Hex2,3[NeuAc2,6]HexNAc
%        ans=965.3336
% 
%    Example 3:
%        glycanMolWt('h(s)n',options); Hex(NeuAc)HexNAc
%        ans=674.2381
%
%    Example 4:
%        slex='fshn'; forFucNeuAcHexHeN
%        glycanMolWt(strcat(slex,'hn'),options)
%        ans=1185.428231
%   
%   Example 5: 
%      glycanstring='hhhhhnn'; % Hex5HexNAc2
%      glycanMolWt(glycanstring);
%     
% See also mwconstant,glycanFormula.

% Author: Gang Liu, Yusen Zhou
% Date Lastly Updated: 05/18/2020

newMolSymbol = [];
options = struct('methylation',true,'acetylation',false,'ion','Na','reductive',false);
nvarargin = numel(varargin);
if rem(nvarargin,2)
    error(message('IncorrectNumberOfArguments'));
end
okargs    = {'optionsstr','newMol'};
for j=1:2:nvarargin
    pname = varargin{j};
    pval = varargin{j+1};
    k = find(strncmpi(pname,okargs,length(pname)));
    if isempty(k)
        error(message('UnknownParameterName'));
    elseif length(k)>1
        error(message('AmbiguousParameterName'));
    else
        switch(k)
            case 1 % optionsstr
                options = pval;
            case 2 % newMol
                newMolSymbol = pval;
        end
    end
end

if(isfield(options,'methylation'))
    ismethyl = options.methylation;
else
    ismethyl = true;
    % error('options is not setup properly');
end

if(isfield(options,'acetylation'))
    isacetyl = options.acetylation;
else
    isacetyl = false;
    % error('options is not setup properly');
end

if(isfield(options,'mono'))
    ismono = options.mono;
else
    ismono = true;
    % error('options is not setup properly');
end

if(isfield(options,'redEnd'))
    isredEnd = options.redEnd;
else
    isredEnd = true;
end

if(isfield(options,'reductive'))
    isreductive = options.reductive;
else
    isreductive = false;
end

if(isfield(options,'ion'))
    ion = options.ion;
else
    ion = 'Na';
end

% handle input argument
% Note: % all masses are internal glycan fragments
if((ismethyl)&&(ismono))
    Hexmwunit    = mwconstant.HexMonoMethyl;  
    HexNACmwunit  = mwconstant.HexNACMonoMethyl;  
    NeuAcmwunit    = mwconstant.NeuACMonoMethyl;  
    NeuGcmwunit    =  mwconstant.NeuGCMonoMethyl;
    DeoHexmwunit  = mwconstant.DeoHexMonoMethyl; 
    Penmwunit  = mwconstant.PenMonoMethyl;
    Sulfmwunit = mwconstant.SulfMonoMethyl;
    Phosmwunit = mwconstant.PhosMonoMethyl;
    KDNmwunit = mwconstant.KDNMonoMethyl;
    HexAmwunit  = mwconstant.HexAMonoMethyl;
    if(~isempty(newMolSymbol))
      newMolmwunit = mwconstant.newMonoMethyl;
    end
    
    
elseif((ismethyl)&&(~ismono))
    Hexmwunit =mwconstant.HexAvgMethyl; 
    HexNACmwunit =mwconstant.HexNACAvgMethyl;  
    NeuAcmwunit =mwconstant.NeuACAvgMethyl;  
    NeuGcmwunit =mwconstant.NeuGCAvgMethyl;
    DeoHexmwunit =mwconstant.DeoHexAvgMethyl; 
    Penmwunit = mwconstant.PenAvgMethyl;
    Sulfmwunit  = mwconstant.SulfAvgMethyl;
    Phosmwunit  = mwconstant.PhosAvgMethyl;
    KDNmwunit  = mwconstant.KDNAvgMethyl;
    HexAmwunit =mwconstant.HexAAvgMethyl; 
    if(~isempty(newMolSymbol))
      newMolmwunit = mwconstant.newAvgMethyl;
    end
    
elseif((isacetyl)&&(ismono))
    Hexmwunit    = mwconstant.HexMonoAcetyl;  
    HexNACmwunit  = mwconstant.HexNACMonoAcetyl;  
    NeuAcmwunit    = mwconstant.NeuACMonoAcetyl;  
    NeuGcmwunit    =  mwconstant.NeuGCMonoAcetyl;
    DeoHexmwunit  = mwconstant.DeoHexMonoAcetyl; 
    Penmwunit  = mwconstant.PenMonoAcetyl;
    Sulfmwunit = mwconstant.SulfMonoAcetyl;
    Phosmwunit = mwconstant.PhosMonoAcetyl;
    KDNmwunit = mwconstant.KDNMonoAcetyl;
    HexAmwunit  = mwconstant.HexAMonoAcetyl; 
    if(~isempty(newMolSymbol))
      newMolmwunit = mwconstant.newMonoAcetyl;
    end
    
elseif((isacetyl)&&(~ismono))
    Hexmwunit =mwconstant.HexAvgAcetyl; 
    HexNACmwunit =mwconstant.HexNACAvgAcetyl;  
    NeuAcmwunit =mwconstant.NeuACAvgAcetyl;  
    NeuGcmwunit =mwconstant.NeuGCAvgAcetyl;
    DeoHexmwunit =mwconstant.DeoHexAvgAcetyl; 
    Penmwunit = mwconstant.PenAvgAcetyl;
    Sulfmwunit  = mwconstant.SulfAvgAcetyl;
    Phosmwunit  = mwconstant.PhosAvgAcetyl;
    KDNmwunit  = mwconstant.KDNAvgAcetyl;
    HexAmwunit =mwconstant.HexAAvgAcetyl; 
    if(~isempty(newMolSymbol))
      newMolmwunit = mwconstant.newAvgAcetyl;
    end
    
elseif((~ismethyl)&&(~isacetyl)&&(ismono))
    Hexmwunit =mwconstant.HexMonoUnde ;  
    HexNACmwunit =mwconstant.HexNACMonoUnde ; 
    NeuAcmwunit =mwconstant.NeuACMonoUnde ;  
    NeuGcmwunit =mwconstant.NeuGCMonoUnde ;
    DeoHexmwunit =mwconstant.DeoHexMonoUnde;  
    Penmwunit=mwconstant.PenMonoUnde ;
    Sulfmwunit=mwconstant.SulfMonoUnde ;
    Phosmwunit=mwconstant.PhosMonoUnde ;
    KDNmwunit=mwconstant.KDNMonoUnde ;
    HexAmwunit =mwconstant.HexAMonoUnde ;  
    if(~isempty(newMolSymbol))
      newMolmwunit = mwconstant.newMonoUnde;
    end
    
elseif((~ismethyl)&&(~isacetyl)&&(~ismono))
    Hexmwunit =mwconstant.HexAvgUnde;  
    HexNACmwunit =mwconstant.HexNACAvgUnde; 
    NeuAcmwunit =mwconstant.NeuACAvgUnde; 
    NeuGcmwunit =mwconstant.NeuGCAvgUnde;
    DeoHexmwunit =mwconstant.DeoHexAvgUnde; 
    Penmwunit =mwconstant.PenAvgUnde;
    Sulfmwunit=mwconstant.SulfAvgUnde;
    Phosmwunit=mwconstant.PhosAvgUnde;
    KDNmwunit=mwconstant.KDNAvgUnde;
    HexAmwunit =mwconstant.HexAAvgUnde;  
    if(~isempty(newMolSymbol))
      newMolmwunit = mwconstant.newAvgUnde;
    end
     
end

if(isredEnd)
    if(isreductive)
        if(ismethyl)
            glycanMW=mwconstant.RedredEndMethyl;
        elseif(isacetyl)
            glycanMW=mwconstant.RedredEndAcetyl;
        else
            glycanMW=mwconstant.RedredEnd;
        end
    elseif(~isreductive)
        if(ismethyl)
            glycanMW=mwconstant.redEndMethyl;
        elseif(isacetyl)
            glycanMW=mwconstant.redEndAcetyl;
        else
            glycanMW=mwconstant.redEnd;
        end
    end
else
    glycanMW=0;
end

% if ion options
if(strcmpi(ion,'none'))
    % do nothing
elseif(strcmpi(ion,'Na'));    
    if(ismono)
       glycanMW=glycanMW+mwconstant.NaionMono;
    else
      glycanMW= glycanMW+mwconstant.NaionAvg;
    end    
elseif(strcmpi(ion,'K'));    
    if(ismono)
       glycanMW=glycanMW+mwconstant.KionMono;
    else
      glycanMW= glycanMW+mwconstant.KionAvg;
    end 
elseif(strcmpi(ion,'H+'));    
    if(ismono)
       glycanMW=glycanMW+mwconstant.HionMono;
    else
      glycanMW= glycanMW+mwconstant.HionAvg;
    end
elseif(strcmpi(ion,'H-'));    
    if(ismono)
       glycanMW=glycanMW-mwconstant.HionMono;
    else
      glycanMW= glycanMW-mwconstant.HionAvg;
    end  
elseif(strcmpi(ion,'Na+H'));    
    if(ismono)
       glycanMW=glycanMW+mwconstant.NaionMono+mwconstant.HionMono;
    else
      glycanMW= glycanMW+mwconstant.NaionAvg+mwconstant.HionAvg;
    end 
elseif(strcmpi(ion,'K+H'));    
    if(ismono)
       glycanMW=glycanMW+mwconstant.KionMono+mwconstant.HionMono;
    else
      glycanMW= glycanMW+mwconstant.KionAvg+mwconstant.HionAvg;
    end 
end

%Hex
nHex=findstr('h',gly);
glycanMW=glycanMW+length(nHex)*Hexmwunit;

% HexNAC
nHexNAc=findstr('n',gly);
glycanMW=glycanMW+length(nHexNAc)*HexNACmwunit;

% Sailic Acid
nNeuAc=findstr('s',gly);
glycanMW=glycanMW+length(nNeuAc)*NeuAcmwunit;

%NEUGC
nNeuGc=findstr('g',gly);
glycanMW=glycanMW+length(nNeuGc)*NeuGcmwunit;

%DeoxyHex
nFuc=findstr('f',gly);
glycanMW=glycanMW+length(nFuc)*DeoHexmwunit;

%Xylose
nXyl=findstr('x',gly);
glycanMW=glycanMW+length(nXyl)*Penmwunit;

%Sulfation
nSO3=findstr('z',gly);                  % sulfation adds SO3 to MW
glycanMW=glycanMW+length(nSO3)*Sulfmwunit;

%Phosphorylation
nPO3=findstr('p',gly);                  % phosphoration adds PO3H to MW
glycanMW=glycanMW+length(nPO3)*Phosmwunit;

%KDN
nKDN=findstr('k',gly);
glycanMW=glycanMW+length(nKDN)*KDNmwunit;

%HEXA
nHexA=findstr('u',gly);
glycanMW=glycanMW+length(nHexA)*HexAmwunit;

if(~isempty(newMolSymbol))
    %User-designed Mol
    nMol=findstr(newMolSymbol,gly);
    glycanMW=glycanMW+length(nMol)*newMolmwunit;
end

end