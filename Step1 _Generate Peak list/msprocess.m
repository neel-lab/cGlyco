function [peaklist,pfwhh]= msprocess(msrawdata,min,max,varargin)
%MSPROCESS convert raw peak data to a peak list and their half height width column
% using a four-step method. The method starts with background adjustment,
% followed by the normalization, signal noise removal, and ends with peak
% finding.
%
% [p,pfwhh]= MSPROCESS(msrawdata,min,max) uses default options which
%   are 0.3 for oversegmentationfilter, 5 for heightfilter, true for
%   showplot. 'min' and 'max' value are set up for the processed mass range.
%   The peak list p is a matrix with two columns of peak
%   locations, and peak intensities. The matrix pfwhh has two columns
%   indicating the left and right location of the full width height.
%
% [p,pfwhh]= MSPROCESS(msrawdata,min, max, option) uses options to assign the
%   parameters for mass spectra processing. 
%

% See also readMS.

% Author: Gang Liu, Yusen Zhou
% Date Lastly Updated: 5/18/20
if(length(varargin)==1)
    options=varargin{1};
    if(isfield(options,'unit'))
        unit = options.unit;
    else
        unit = 'Da';
    end
    if(strcmp(unit,'Da'))&&(isfield(options,'MA'))
        overSegmentationFilter = 0.9 * options.MA;
    elseif(strcmp(unit,'ppm'))&&(isfield(options,'MA'))
        overSegmentationFilter = 0.6;
    else
        overSegmentationFilter = 0.6;
    end
    
    if(isfield(options,'heightFilter'))
        heightFilter = options.heightFilter;
    else
        heightFilter = 5;
    end
    
    if(isfield(options,'showPlot'))
        showPlot = options.showPlot;
    else
        showPlot = true;
    end
    if(isfield(options,'span'))
        span = options.span;
    else
        span = 10;
    end
else
    overSegmentationFilter = 0.6;
    heightFilter = 0.01;
    showPlot     = false;
    span  = 10;
end

chosenindexMin  = find((msrawdata(:,1)-min)>=0);
startindex      = chosenindexMin(1);
if(isnumeric(max))
    chosenindexMax  = find((msrawdata(:,1)-max)<=0);
    endindex        = chosenindexMax(length(chosenindexMax));
elseif(ischar(max))
    endindex = chosenindexMin(length(chosenindexMin));
end
mz           = msrawdata(startindex:endindex,1);
intensity    = msrawdata(startindex:endindex,2);
% msviewer(mz,intensity);

% first step baseline correction
winsize = 200;
[bsCorrectedInt] = msbackadj(mz, intensity,...
    'WindowSize', winsize,...
    'RegressionMethod', 'spline',...
    'ShowPlot',false,'Quantile',0.10);
%msviewer(mz,bsCorrectedInt,'ylabel','Baseline Corrected');

% remove noise signal smooth signal using mssgolay (polynomial filters)
noiseremInt = mssgolay(mz,bsCorrectedInt,...
    'span',span,'showplot',false);

% second step normalization
normInt = msnorm(mz,noiseremInt,...
    'LIMITS',[0 inf],'MAX',100);

%  if(false)
%    msviewer(mz,normInt,'ylabel','Normalized');
%  end

% peak finding with wavelets denoising
[peaklist,pfwhh]= mspeaks(mz,normInt,...
    'DENOISING',false,'HEIGHTFILTER',heightFilter,...
    'OverSegmentationFilter', overSegmentationFilter,...
    'SHOWPLOT',showPlot);%'FWHHFILTER',0.1
end