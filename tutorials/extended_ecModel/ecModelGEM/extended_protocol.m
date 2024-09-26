%% Moving to adapter file location
projectLoc = 'C:\Users\blank\Documents\VU\GECKO\tutorials\extended_ecModel\ecModelGEM';
cd(projectLoc);

%% Start the project
%startGECKOproject();
adapterLocation = fullfile('.\ecMycobacteriumGEMAdapter.m'); 
ModelAdapterManager.setDefault(adapterLocation);

ModelAdapter = ModelAdapterManager.getDefault();
params = ModelAdapter.getParameters();

%% Load model
model = loadConventionalGEM(); %load the GEM if the location is given in
%the adapter
%model = importModel(fullfile('C:\Users\blank\Documents\VU\GECKO\tutorials\extended_ecModel\ecModelGEM\models\iEK1008.xml'));

%% Begin expansion
[ecModel, noUniprot] = makeEcModel(model); %make the enzyme constrained
%model
%saveEcModel(ecModel,'ecMycoGEM_stage1.yml'); %save the enzyme constrained
%model to 'models' to load later
%ecModel=loadEcModel('ecMycoGEM_stage1.yml'); %load the enzme constrained model

%% Query database for BLAST: uniprotGEM.fasta
% Give the extended model first and the GEM second and it writes a 
% uniprotGEM.fasta file based on the intersection
subsetGenomeBasedOnGEM(ecModel,model);

%% Retrieve complex data
% Call getComplexData function like this:variable = getComplexData(0)
% and give '[]' in the adapter
% This will generate the 'ComplexPortal.json' file.
complexInfo = getComplexData(0); % Downloads all complex data for all
% organisms

%% IF THERE IS NO COMPLEXINFO THEN THE APPLICATION SHOULD BE READ IN FROM THE FILE
%% AND IT SHOULD ALSO BE ABLE TO HANDLE WHEN THERE IS COMPLEXINFO

%% Local database based on Complex Portal data: ComplexPortalDB.fasta
% 
[noSeq] = getComplexLocalDB(complexInfo);

%% BLAST
%PAY ATTENTION TO GIVE THE CORRECT FILES IN ORDER
%Give a name for the database that contain protein IDs
organismID = 'database';
%give the location of the database
fastaFile = fullfile(projectLoc,'data','ComplexPortalDB.fasta');
%give a name for the query that contains the metabolic genes
modelIDs = 'query';
%give the location of the query
refFastaFiles = fullfile(projectLoc,'data','uniprotGEM.fasta');
%Set developMode true
developMode = true;

[blastStructure, blastReport] = getBlast(organismID, fastaFile, modelIDs, refFastaFiles, developMode);

%% fakeComplexInfo
minAligLen = 150;
fakeComplexInfo = getBlastComplexInfo(blastStructure, minAligLen, complexInfo, ModelAdapter);

%% Making sure fakeComplexInfo only contains genes from the model
modelGenes = ecModel.ec.genes;
complexGenes = fakeComplexInfo.geneName;

for i = 1:numel(complexGenes)
    [Lia, Locb] = ismember(complexGenes(i), modelGenes);
    fakeComplexInfo(i).complexID(~Lia) = [];
    fakeComplexInfo(i).name(~Lia) = [];
    fakeComplexInfo(i).species(~Lia) = [];
    fakeComplexInfo(i).geneName(~Lia) = [];
    fakeComplexInfo(i).protID(~Lia) = [];
    fakeComplexInfo(i).stochiometry(~Lia) = [];
    fakeComplexInfo(i).defined(~Lia) = [];
end

%% Histogram
% Convert structure to table
fakeComplexTable = struct2table(fakeComplexInfo);

% Access the 'stochiometry' field from the table
stoichiometryData = fakeComplexTable.stochiometry;

% Filter values that are greater than 1
plottedStochiometry = stoichiometryData(stoichiometryData > 1);

% Get the histogram data (bins and counts)
[counts, edges] = histcounts(stoichiometryData);

% Plot using the 'bar' function instead of 'histogram'
binCenters = edges(1:end-1) + diff(edges)/2;  % Calculate bin centers
bar(binCenters, counts, 'FaceColor', [0.5, 0.5, 0.5]);

title('Stoichiometry frequencies based on Complex Portal', 'FontSize', 14)
xlabel('Stoichiometry values', 'FontSize', 12)
ylabel('Frequency', 'FontSize', 12)

% Display the frequency on top of each bar
for i = 1:length(counts)
    if counts(i) > 0  % Only show text for non-zero counts
        text(binCenters(i), counts(i), num2str(counts(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end

savepath = 'C:/Users/blank/Documents/VU/Research_project/thesis/histogram_stochiometry.pdf';

exportgraphics(gcf, savepath, 'Resolution',600);
%% Apply fakeComplexInfo
[complexModel, foundComplex, proposedComplex] = applyComplexData(ecModel, fakeComplexInfo, ModelAdapter);

%% Writing the matrix into a file to look into it
% Step 1: Extract the rxnEnzMat matrix
%complexMatrixToFile = complexModel.ec.rxnEnzMat;
%ecMatrixToFile = ecModel.ec.rxnEnzMat;


% Step 2: Write the matrix to an Excel file
%writematrix(complexMatrixToFile, 'nonzerocomplexrxnEnzMat.xlsx');
%writematrix(ecMatrixToFile, 'nonzeroecrxnEnzMat.xlsx');
%% 
% Adding pairs for comparison
% Query against Database
% blastStructure(1).pairs = strcat(blastStructure(1).fromGenes, '-', blastStructure(1).toGenes);
% Database against query
% blastStructure(2).pairs = strcat(blastStructure(2).toGenes, '-', blastStructure(2).fromGenes);

% checking unique values in pairs before ismembering

% Find the unique pairs
% uniquePairs = unique(blastStructure(1).pairs, 'stable');
% uniqueFromGenes = unique(blastStructure(1).fromGenes, 'stable');


% Get the number of unique pairs
% numUniquePairs = numel(uniquePairs);

% Display the count
% disp(['Number of unique pairs in blastStructure(1) before ismember: ', num2str(numUniquePairs)]);

% uniquePairs = unique(blastStructure(2).pairs, 'stable');

% Get the number of unique pairs
% numUniquePairs = numel(uniquePairs);

% Display the count
% disp(['Number of unique pairs in blastStructure(2) before ismember: ', num2str(numUniquePairs)]);

%% This section tries to find the intersection, but with that, we lose the
%% information of the bidirectional BLAST
%{
% only keep the information if a match is found
Lia1 = ismember(blastStructure(1).pairs, blastStructure(2).pairs);

blastStructure(1).fromGenes = blastStructure(1).fromGenes(Lia1);
blastStructure(1).toGenes = blastStructure(1).toGenes(Lia1);
blastStructure(1).evalue = blastStructure(1).evalue(Lia1);
blastStructure(1).identity = blastStructure(1).identity(Lia1);
blastStructure(1).aligLen = blastStructure(1).aligLen(Lia1);
blastStructure(1).bitscore = blastStructure(1).bitscore(Lia1);
blastStructure(1).ppos = blastStructure(1).ppos(Lia1);
blastStructure(1).pairs = blastStructure(1).pairs(Lia1);
% use the filtered blastStructure(1) for comparison
Lia2 = ismember(blastStructure(2).pairs, blastStructure(1).pairs);

blastStructure(2).fromGenes = blastStructure(2).fromGenes(Lia2);
blastStructure(2).toGenes = blastStructure(2).toGenes(Lia2);
blastStructure(2).evalue = blastStructure(2).evalue(Lia2);
blastStructure(2).identity = blastStructure(2).identity(Lia2);
blastStructure(2).aligLen = blastStructure(2).aligLen(Lia2);
blastStructure(2).bitscore = blastStructure(2).bitscore(Lia2);
blastStructure(2).ppos = blastStructure(2).ppos(Lia2);
blastStructure(2).pairs = blastStructure(2).pairs(Lia2);
%}
%%
% Find the unique pairs
% uniquePairs = unique(blastStructure(1).pairs, 'stable');

% Get the number of unique pairs
% numUniquePairs = numel(uniquePairs);

% Display the count
% disp(['Number of unique pairs in blastStructure(1): ', num2str(numUniquePairs)]);

% uniquePairs = unique(blastStructure(2).pairs, 'stable');

% Get the number of unique pairs
% numUniquePairs = numel(uniquePairs);

% Display the count
% disp(['Number of unique pairs in blastStructure(2): ', num2str(numUniquePairs)]);


% checking for duplicates
%[uniquePairs1, ~, idx1] = unique(blastStructure(1).pairs);
%[uniquePairs2, ~, idx2] = unique(blastStructure(2).pairs);

%duplicates1 = blastStructure(1).pairs(setdiff(1:numel(blastStructure(1).pairs), idx1));
%duplicates2 = blastStructure(2).pairs(setdiff(1:numel(blastStructure(2).pairs), idx2));

%disp('Duplicate pairs in blastStructure(1):');
%disp(duplicates1);

%disp('Duplicate pairs in blastStructure(2):');
%disp(duplicates2);

%% tried batches, didnt work
%{
batchIDs = strjoin(keys(1:idsToProcess), ',');

url = ['https://rest.uniprot.org/uniprotkb/accessions?accessions=' batchIDs '&format=fasta'];
options = weboptions('ContentType', 'binary', 'RequestMethod', 'get');

try
    % Fetch compressed data for all IDs
    compressedFastaData = webread(url, options);

    % Save the compressed data to a temporary file
    tempGzipFile = [tempname, '.gz'];
    fid = fopen(tempGzipFile, 'wb');  % Open for writing binary data
    fwrite(fid, compressedFastaData);
    fclose(fid);

    % Decompress the gzip file
    gunzip(tempGzipFile);
    decompressedFile = tempGzipFile(1:end-3);  % Remove .gz extension

    % Read the decompressed FASTA data
    fid = fopen(decompressedFile, 'r');
    fastaData = fread(fid, '*char')';
    fclose(fid);

    % Append the FASTA data to the output file
    fid = fopen(outputFilePath, 'a');  % Open in append mode
    fprintf(fid, '%s', fastaData);
    fclose(fid);

    % Clean up temporary files
    delete(tempGzipFile);
    delete(decompressedFile);

catch ME
    warning('Error fetching sequences from UniProt: %s', ME.message);
end

%}

%% Model-specific ComplexPortal database stored at
%C:\Users\blank\Documents\VU\GECKO\tutorials\extended_ecModel\ecModelGEM\data\ComplexPortal.json
%[ecModel, foundComplex, proposedComplex] = applyComplexData(ecModel, complexInfo);

%% Extract the EC numbers 
%[ecModel, invalidECs, invalidPos] = getECfromGEM(ecModel);
%noEC = cellfun(@isempty, ecModel.ec.eccodes); 
%ecModel = getECfromDatabase(ecModel, noEC); 

%% perform search in BRENDA
%kcatList_fuzzy = fuzzyKcatMatching(ecModel);
%save('./data/kcatList_fuzzy.mat', 'kcatList_fuzzy');

%load('./data/kcatList_fuzzy.mat')
%kcat = table(kcatList_fuzzy.rxns, kcatList_fuzzy.origin, kcatList_fuzzy.wildcardLvl,...
           % 'VariableNames',{'Rxn','Origin','WildcardLevel'});