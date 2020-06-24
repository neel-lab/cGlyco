function  [peaktotalarea,matchedpeakindex]= findpeakarea(md,isotopicindex,isopeaknum,peaklist,pfwhh,matchedpeakindex)
% findoverlaparea: find peak area when there is only one glycan
%
% See also msfraction.
%
% Author: Yusen Zhou
% Data Lastly Updated: 05/20/2020

dist = md(:,2);
isotopicindexarea = abs(pfwhh(isotopicindex,1)-pfwhh(isotopicindex,2))...
    *peaklist(isotopicindex,2);

peakarea =isotopicindexarea;

peakindex = isotopicindex;
totaldist = dist(1,1);
for i = 2 : isopeaknum
    peakindex = peakindex+1;
    try
        if((peaklist(peakindex,1)-peaklist(peakindex-1,1))>2)
            break
        end
    catch
        break
    end
    listnumber = length(matchedpeakindex);
    if(isempty(find((matchedpeakindex-peakindex)==0, 1)))
        matchedpeakindex(listnumber+1) = peakindex;
    end
    totaldist = totaldist+dist(i,1);
    peakindexarea =abs(pfwhh(peakindex,1)-pfwhh(peakindex,2))...
        *peaklist(peakindex,2);
    peakarea  =  peakarea + peakindexarea;
end
peaktotalarea{1} = peakarea/totaldist;
end
