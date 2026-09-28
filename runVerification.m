function runVerification
%RUNVERIFICATION Run all four legacy suites and the new starfish regressions.
% Each suite runs with a clean MATLAB path and its own current directory,
% so the QC-ring and tree mmsCases helpers cannot shadow each other.
root=fileparts(mfilename('fullpath')); oldPath=path; oldDir=pwd;
cleanup=onCleanup(@() restore(oldPath,oldDir)); %#ok<NASGU>
cases={...
 'EllipticBSAM/tests', 'testEbNeumannSemantics; runAllTests'; ...
 'EllipticBSAMTree/tests', 'testEbtNeumannSemantics; runAllTreeTests'; ...
 'ParabolicBSAM/tests', 'runParabolicRegressions; runParabolicTests'; ...
 'DiffuseDomainBSAM/tests', 'runDdRegressions; runDdTests'; ...
 'tests', 'testStarfish'};
for j=1:size(cases,1)
    restoredefaultpath;
    runOne(fullfile(root,cases{j,1}),cases{j,2});
    cd(root);
end
fprintf('\nALL BSAM PORTFOLIO VERIFICATION SUITES PASSED\n');
end

function runOne(folder,command)
cd(folder); eval(command);
end

function restore(oldPath,oldDir)
path(oldPath); cd(oldDir);
end
