function [MSdata,fitPara]=peakList(fname,varargin)
% [MSdata,firstPeak,wf]=peakList(msrawdatafilename,varargin)
% peakList: loads raw MS data and use msprocess code to convert raw peak 
% data to a peak list and their half height width column (FWHM)
%
% MSdata = peakList(fname) use the default value of msrawdatapath
% and options to convert 'msrawdatafilename' file, a '.txt' file, into matlab form 
% and then preprocesses it by using a four-step method to finally get the peak 
% list and their FWHM.
%
% MSdata = peakList(fname,'msrawdatapathstr',msrawdatapath) find the 'msrawdatafilename' 
% file in msrawdatapath and converts it into matlab form and then preprocesses it 
% by using a four-step method to finally get the peak list and FWHM.
%
% MSdata = peakList(fname,'optionsstr',options) sets up your 
% own options and converts msrawdatafilename file into matlab form and then 
% preprocesses it using a four-step method to finally get the peak list and
% FWHM. 
%
% See also MSprocess
% This requires the MATLAB Bioinformatics toolbox

% Author: Yusen Zhou & Gang Liu
% Date Lastly Updated: 6/25/20


% default value
dir = 'D:\Work\cGlyco\Testcase\Lung';
options.MA                     = 0.9;
options.heightFilter           = 0.003;
options.unit                   = 'Da';
min                            = 1400;
max                            = 5000;
SN                             = 1;

nvarargin = numel(varargin);
if rem(nvarargin,2)
    error(message('IncorrectNumberOfArguments'));
end
okargs    = {'dir','optionsstr','min','max','SN'};
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
            case 1 % msrawdatapathstr
                dir = pval;
            case 2 % optionsstr
                options.MA           = pval.MA;
                options.heightFilter = pval.heightFilter;
                options.unit         = pval.unit;
            case 3 % min
                min = pval;
            case 4 % max
                max = pval;
        end
    end
end

% read MS raw data
fitPara = cell(4,1);
filename  = [fname '.txt'];
msrawdatafullfilename = fullfile(dir,filename);
data      = readMS(msrawdatafullfilename);
% Preprocess of MS data
[peaklist,FWHM] = msprocess(data,min,max,options);
[noisepeaklist,peaklist,FWHM] = modifyPeaks(peaklist,FWHM,SN);
try
    % Non-linear
    fitFun = fit(noisepeaklist(:,1),noisepeaklist(:,2),'exp2');
    fitPara{1} = fitFun.a;
    fitPara{2} = fitFun.b;
    fitPara{3} = fitFun.c;
    fitPara{4} = fitFun.d;
catch
end
MSdata   = struct('peaklist',peaklist,'FWHM',FWHM);
filepath = fullfile(dir,fname);
save(filepath,'MSdata');
end

function [noisepeaklist,peaklist,pfwhh] = modifyPeaks(peaklist,pfwhh,SN)
peakindex    = 0;
count        = 0;
startnum     = 1;
selectpeakindex = [];
while (peakindex <= length(peaklist)-3)
    peakindex   = peakindex+1;
    if(peakindex>startnum+49)
        startnum = peakindex;
    end
    % Take median of every 50 peaks as local noise level
    if(startnum+49<=length(peaklist))
        localnoise = median(peaklist(startnum:startnum+49,2));
    else
        localnoise = median(peaklist(startnum:length(peaklist),2));
    end
    % Consider three successive peaks
    localratio1 = peaklist(peakindex,2)/localnoise;
    localratio2 = peaklist(peakindex+1,2)/localnoise;
    localratio3 = peaklist(peakindex+2,2)/localnoise;
    if(localratio1<1.5&&localratio2<1.5&&localratio3<1.5)
        selectpeakindex(count+1) = peakindex;
        count = count+1;
    else
        if(localratio1>2)
            candidatepeak = peakindex;
        elseif(localratio2>2)
            candidatepeak = peakindex+1;
        else
            candidatepeak = peakindex+2;
        end
        maxIntensity  = peaklist(candidatepeak,2);
        breakpoint   = 0;
        % the intensity of last peak included should be >=15% of the most
        % intensity peak with local area
        while(breakpoint==0)
            candidatepeak = candidatepeak+1;
            if(candidatepeak<=length(peaklist))
                if((peaklist(candidatepeak,2)/maxIntensity)<0.15||peaklist(candidatepeak,2)<=localnoise)
                    breakpoint = 1;
                else
                    maxIntensity = max(maxIntensity,peaklist(candidatepeak,2));
                end
            else
                breakpoint = 1;
            end
        end  
        peakindex = candidatepeak-1;
    end
end
selectpeakindex(count+1)       = length(peaklist)-2;
selectpeakindex(count+2)       = length(peaklist)-1;
selectpeakindex(count+3)       = length(peaklist);
noisepeaklist                  = peaklist(selectpeakindex,:);
if(SN)
    peaklist(selectpeakindex,:)    = '';
    pfwhh(selectpeakindex,:)       = '';
end
end