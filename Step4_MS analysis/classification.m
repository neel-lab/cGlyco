function newglycanDB = classification(MSfilename,glycanDB,loadpath,storepath)
% classification: Monosaccharide and substructure comparison analysis
%
% Input:
%   MSfilename: MS data for glycan profiles.
%     glycanDB: glycan structure data with relative abundance information
%     loadpath: Directory to load glycanDB.
%    storepath: Directory to output results.
%
% Author: Yusen Zhou
% Date Lastly Updated: 05/18/2020


glycanDBname = [glycanDB '.mat'];
glycanDBfullpath   = fullfile(loadpath,glycanDBname);
load(glycanDBfullpath,'newglycanDB');
ResidueResult = MSResidueAnalysisGNAT(MSfilename,newglycanDB);
StructResult = MSStructAnalysisGNAT(MSfilename,newglycanDB);
newglycanDB.MonoAnalysis   = ResidueResult;
newglycanDB.StructAnalysis = StructResult;
MSfilename = [MSfilename 'FinalResult' '.mat'];
savepath = fullfile(storepath, MSfilename);
save(savepath,'newglycanDB');
end