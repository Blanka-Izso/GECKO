function blastFakeComplexInfo = getBlastComplexInfo(blastStructure, minLen, complexInfo, modelAdapter)
% getBlastComplexInfo
%   Creates a fakeComplexInfo structure based on bidirectional blast
%   results for applyComplexData() function. 
% Input:
%   blastStructure:     structure containing one bidirectional homology
%                       measurements (it only works with one to one blast)
%   complexInfo:        a structure as generated by getComplexData.
%   minLen:             a minimum alignment lenght to filter out hits below 
%                       it from the blastStructure   
%   modelAdapter:       a loaded modeladapter (Optional, will otherwise use
%                       the default model adapter.)
% Output:
%   blastFakeComplexInfo: structure following complexInfo structure except 
%                       each protein is handled as a separate complex
% Usage:
%  blastFakeComplexInfo = getBlastComplexInfo(blastStructure, minLen,
%  complexInfo, modelAdapter)

if nargin < 4 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end

if nargin<3 || isempty(complexInfo)
    error('No complexInfo provided. Please call getComplexData() before usage')
else
    complexData = complexInfo;
end

if isempty(blastStructure)
    error('BlastStructure was not given. Please call getBlast() before usage')
end

if isempty(minLen)
    error('Minimum alignment length was not given for filtering.')
end



%% Using BLAST info, filtering
% Handling the fasta headers to only contain protein ids
for i = 1:2
    blastStructure(i).fromGenes = cellfun(@(x) extractProteinID(x), blastStructure(i).fromGenes, 'UniformOutput', false);
    blastStructure(i).toGenes = cellfun(@(x) extractProteinID(x), blastStructure(i).toGenes, 'UniformOutput', false);
end

% The extracting function
function proteinID = extractProteinID(header)
    parts = strsplit(header, '|');
    %if the header is a fasta header it will have at least 3 parts
    if length(parts) >= 3
        proteinID = parts{2};
    else %if not valid fasta header format return original
        proteinID = header;
    end
end

%filter based on thresholds
% getBlast() removes hits with an E-value higher than 10e-5, so we do not 
% need to filter based on E-values
% Define the minimum alignment length for filtering

%Remove all gene matches that are below the cutoff values
for i = 1:numel(blastStructure)
    indexes = blastStructure(i).aligLen >= minLen;
    blastStructure(i).fromGenes(~indexes) = [];
    blastStructure(i).toGenes(~indexes) = [];
    blastStructure(i).evalue(~indexes) = [];
    blastStructure(i).identity(~indexes) = [];
    blastStructure(i).aligLen(~indexes) = [];
    blastStructure(i).bitscore(~indexes) = [];
    blastStructure(i).ppos(~indexes) = [];
end

% Combine struct to make the trace back and application of complex data 
% possible

% Get the lengths of the two structures for the 2 directions
n1 = length(blastStructure(1).fromGenes);
n2 = length(blastStructure(2).fromGenes);

% Initialize the combined structure
combinedBlastStructure = struct();

% Create the 'direction' field for both structures
direction1 = [blastStructure(1).fromId '-' blastStructure(1).toId];
direction2 = [blastStructure(2).fromId '-' blastStructure(2).toId];

combinedBlastStructure.direction = [repmat({direction1}, n1, 1); ...
                                    repmat({direction2}, n2, 1)];

% Concatenate the remaining fields
% metGeneID contains the ORF IDs, the metabolic genes
% databaseProteinID contains the protein IDs, the uniprot IDs from the
% fasta headers based on the complex portal database
combinedBlastStructure.metGeneID            = [blastStructure(1).fromGenes; blastStructure(2).toGenes];
combinedBlastStructure.databaseProteinID    = [blastStructure(1).toGenes; blastStructure(2).fromGenes];
combinedBlastStructure.evalue               = [blastStructure(1).evalue; blastStructure(2).evalue];
combinedBlastStructure.identity             = [blastStructure(1).identity; blastStructure(2).identity];
combinedBlastStructure.aligLen              = [blastStructure(1).aligLen; blastStructure(2).aligLen];
combinedBlastStructure.bitscore             = [blastStructure(1).bitscore; blastStructure(2).bitscore];
combinedBlastStructure.ppos                 = [blastStructure(1).ppos; blastStructure(2).ppos];
%adding pairs to find unique pairs and keep the ones with the lowest
%evalue
%combinedBlastStructure.pairs = strcat(combinedBlastStructure.metGeneID, '-', combinedBlastStructure.databaseProteinID);

% Display the result
%disp(combinedBlastStructure);
%% Filtering again to gain one to one links based on the genes
% Keeping only the unique geneIDs with the lowest e-value

% Convert the combined structure to a table
blastTable = struct2table(combinedBlastStructure);

% Sort the table by 'metGeneID' and 'evalue' (ascending order), default is
% ascending order - lowest to highest
% metGeneIDs are grouped together in ascending order of the evalue
sortedTable = sortrows(blastTable, {'metGeneID', 'evalue'});

% Find unique metGeneID and their first occurrences (which will be the lowest e-value due to sorting)
[~, uniqueIdx] = unique(sortedTable.metGeneID, 'first');

% Select rows corresponding to unique metGeneID with the lowest e-value
filteredTable = sortedTable(uniqueIdx, :);

% Convert the table back to a structure
filteredBlastStructure = table2struct(filteredTable);

%%

% Ddding stochiometry and other info from complexInfo
% Initialize new fields in filteredBlastStructure
[filteredBlastStructure(1:length(filteredBlastStructure)).complexID] = deal('');
[filteredBlastStructure(1:length(filteredBlastStructure)).name] = deal('');
[filteredBlastStructure(1:length(filteredBlastStructure)).species] = deal('');
[filteredBlastStructure(1:length(filteredBlastStructure)).stochiometry] = deal(0);
[filteredBlastStructure(1:length(filteredBlastStructure)).defined] = deal(3);  % Initialize with 3 as placeholder

%looping through the filteredBlastStructure
for i=1:numel(filteredBlastStructure)
    %getting the current uniprot ID from filteredBlastStructure
    currProtID = filteredBlastStructure(i).databaseProteinID;
    %Looping through complexInfo to find matching protID
    for j=1:numel(complexData)
        %Getting the proteins and their stochiometries
        protIDs = complexData(j).protID;
        %
        stochiometries = complexData(j).stochiometry;
        %Checking if the current protein ID is in this complex
        matchIndex = find(strcmp(protIDs, currProtID));

        if ~isempty(matchIndex)
            %disp('matchindex is not empty')
            %disp('i ')
            %disp(i)
            %disp('j ')
            %disp(j)
            %if match is found, assign the corresponding fields
            % Assign strings: complexID, name, species
            filteredBlastStructure(i).complexID = complexData(j).complexID;
            filteredBlastStructure(i).name = complexData(j).name;
            filteredBlastStructure(i).species = complexData(j).species;
            %if match found, assign the corresponding stochiometry
            %disp(complexInfo(j).complexID);
            filteredBlastStructure(i).stochiometry = stochiometries(matchIndex);

            % Assign defined (numerical field)
            filteredBlastStructure(i).defined = complexData(j).defined;
            break; %Stopping to look for the same protein if we found it in one complex

        else 
        % If no match is found, keep initialized values (already set)
        % Optional: You could explicitly reset if needed
        %disp('I havent found anything bruhuhu')
        filteredBlastStructure(i).complexID = '';
        filteredBlastStructure(i).name = '';
        filteredBlastStructure(i).species = '';
        filteredBlastStructure(i).stochiometry = 0; 
        filteredBlastStructure(i).defined = 3; 
        end
    end
end

%%
%Structuring the blast structure as a fakeComplexInfo
fakeComplexInfo = struct('complexID', [], 'name', [], 'species', [], ...
                          'geneName', [], 'protID', [], ...
                          'stochiometry', [], 'defined', []);

% Assign all elements from filteredBlastStructure with non-zero stochiometry
% wrapping the filteredBlastStructure field in cells so applyComplexData
% can handle the structure

% Initialize index counter for valid entries
validIdx = 1;
for i = 1:length(filteredBlastStructure)
    if filteredBlastStructure(i).stochiometry > 0
        fakeComplexInfo(validIdx).complexID   = {filteredBlastStructure(i).complexID};
        fakeComplexInfo(validIdx).name        = {filteredBlastStructure(i).name};
        fakeComplexInfo(validIdx).species     = {filteredBlastStructure(i).species};
        fakeComplexInfo(validIdx).geneName    = {filteredBlastStructure(i).metGeneID};
        fakeComplexInfo(validIdx).protID      = {filteredBlastStructure(i).metGeneID};
        fakeComplexInfo(validIdx).stochiometry = filteredBlastStructure(i).stochiometry;
        fakeComplexInfo(validIdx).defined     = filteredBlastStructure(i).defined;
        % Increment the next valid index
        validIdx = validIdx + 1;
    end
end

blastFakeComplexInfo = fakeComplexInfo;
end
