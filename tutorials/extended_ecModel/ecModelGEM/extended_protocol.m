%% Moving to adapter file location
projectLoc = fullfile(findGECKOroot,'tutorials', 'extended_ecModel','ecModelGEM');
cd(projectLoc);

%% Start the project
%startGECKOproject();
adapterLocation = fullfile('.\ecMycobacteriumGEMAdapter.m'); 
ModelAdapterManager.setDefault(adapterLocation);

ModelAdapter = ModelAdapterManager.getDefault();
params = ModelAdapter.getParameters();

%% Load model
model = loadConventionalGEM(); %

%% Begin expansion
[ecModel, noUniprot] = makeEcModel(model); 

%% Query database for BLAST: uniprotGEM.fasta
% Give the extended model first and the conventional GEM second 
% uniprotGEM.fasta file based on the intersection
%subsetGenomeBasedOnGEM(ecModel,model); %commented out because the file is
%already generated

%% Retrieve complex data
% This will generate the 'ComplexPortal.json' file.
%complexInfo = getComplexData(0); %commented out because the file is
%already generated
%%
%save('./data/complexInfo.mat', 'complexInfo'); % saved for faster build
load('complexInfo.mat')


%% Local database based on Complex Portal data: ComplexPortalDB.fasta
% 
%[noSeq] = getComplexLocalDB(complexInfo); % the file is already generated

%% BLAST
% Correct order for the files is essential
% Give a name for the database that contain protein IDs
organismID = 'database';
% Give the location of the database from Complex Portal
fastaFile = fullfile(projectLoc,'data','ComplexPortalDB.fasta');
% Give a name for the query that contains the generated uniprotGEM.fasta
% file
modelIDs = 'query';
% Give the location of the query
refFastaFiles = fullfile(projectLoc,'data','uniprotGEM.fasta');
% Set developMode true
developMode = true;

[blastStructure, blastReport] = getBlast(organismID, fastaFile, modelIDs, refFastaFiles, developMode);

%% fakeComplexInfo
minAligLen = 150; % Set minimum alignment length for matches
fakeComplexInfo = getBlastComplexInfo(blastStructure, minAligLen, complexInfo, ModelAdapter);



%% Apply fakeComplexInfo
%This will give only full matches because the fakeComplexInfo stores
%information as a one-protein complex to collect stoichiometry
[complexModel, foundComplex, proposedComplex] = applyComplexData(ecModel, fakeComplexInfo, ModelAdapter);

%% Extract the EC numbers 
[ecModel, invalidECs, invalidPos] =     getECfromGEM(ecModel);
noEC = cellfun(@isempty, ecModel.ec.eccodes); 


%% perform search in BRENDA
kcatList_fuzzy = fuzzyKcatMatching(ecModel);
save(fullfile(findGECKOroot,'tutorials', 'extended_ecModel','ecModelGEM', 'data', 'kcatList_fuzzy.mat'), 'kcatList_fuzzy');

load(fullfile(findGECKOroot,'tutorials', 'extended_ecModel','ecModelGEM', 'data', 'kcatList_fuzzy.mat'))
kcat = table(kcatList_fuzzy.rxns,kcatList_fuzzy.kcats, kcatList_fuzzy.origin, kcatList_fuzzy.wildcardLvl,...
            'VariableNames',{'Rxn','KcatValues','Origin','WildcardLevel'});

%% Paste the kcat list into the model
ecModel = selectKcatValue(ecModel, kcatList_fuzzy, 'max');
% Apply the kcat for isozymes w/o evidence
ecModel = getKcatAcrossIsozymes(ecModel);

%% Merge kcat data into the model
ecModel = applyKcatConstraints(ecModel);
ecModel = setProtPoolSize(ecModel); % use the initialized parameters
 
%% Test the model with unlimited nutrient supply
ecModel = setParam(ecModel, 'obj', params.bioRxn, 1);
sol_un = solveLP(ecModel); 
bioRxnIdx = getIndexes(ecModel, params.bioRxn, 'rxns'); 
