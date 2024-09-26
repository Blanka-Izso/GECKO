function [noUniProtSeq] = getComplexLocalDB(complexinfo,modeladapter)
% getComplexLocalDB
%   Creates a 'ComplexPortalDB.fasta' file based on Complex Portal data
%   and gives back the IDs that do not have a sequence in the UniProt 
%   database. If no argument is given, an attempt will be made to read in
%   the ComplexPortal.json file if found
%
% Input:
%   complexinfo:    output of 'getComplexData()' function
%   modeladapter:   
% Output:
%   noUniProtSeq:   list of IDs that does not have a sequence in the
%                   UniProt database


% if modelAdapter is not given, set it to default
if nargin < 2 || isempty(modeladapter)
    modeladapter = ModelAdapterManager.getDefault();
    if isempty(modeladapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end

% Get the parameters from the modelAdapter
params      = modeladapter.getParameters();
% Set the filePath to the given path given in the adapter
filePath    = fullfile(params.path,'data');
% Set the path to where the file should be written
fastaPath = fullfile(filePath,'ComplexPortalDB.fasta');

% If no input argument is provided, or the 'complexinfo' input is empty
% and complexinfo is not a file, load the complexinfo with getComplexData
% if a file is found, it loads it
if nargin < 1 || isempty(complexinfo)
    complexFilePath = fullfile(params.path, 'data', 'ComplexPortal.json');
    if ~isfile(complexFilePath)
        complexinfo = getComplexData(0); %if file not found, load data
    else %if file found, read in the file
        jsonStr = fileread(complexFilePath);
        complexinfo = jsondecode(jsonStr);
        disp(['File named ''ComplexPortal.json'' found at: '...
            complexFilePath 'and loaded.'])
    end
end

%

% Initialize a container for goint through the UniProt IDs
complexPortalProtIDs = containers.Map('KeyType', 'char', 'ValueType', 'double');
wrongID = 0; %count the IDs that should not be included

% Iterate over each element of complexInfo
for i = 1:length(complexinfo)
    % Get the current protID cell array
    currentProtIDs = complexinfo(i).protID;
    
    % Iterate over each UniProt ID in the cell array
    for j = 1:length(currentProtIDs)
        uniprotID = currentProtIDs{j};
        
        % if contains pro, cpx, ebi, dont put it in the map, just increment the
        % counter
        if contains(uniprotID, '-PRO_') || contains(uniprotID,  'CPX-') || contains(uniprotID, 'EBI-')
            wrongID = wrongID + 1;
        else
            if isKey(complexPortalProtIDs, uniprotID) %if it doesnt contain pro and is already in the map, count
                complexPortalProtIDs(uniprotID) = complexPortalProtIDs(uniprotID) + 1;
            else
            % Otherwise, add it to the map with a count of 1
                complexPortalProtIDs(uniprotID) = 1;
            end
        end
    end
end

% DISPLAYING THE MAP
%{
keys = complexPortalProtIDs.keys;
values = complexPortalProtIDs.values;

% Display the contents
for i = 1:length(keys)
    fprintf('%s: %d\n', keys{i}, values{i});
end
%}


keys = complexPortalProtIDs.keys;
missingIDs = {}; %to store failing ids
foundID = 0; %count how many ids are being written in the file

if exist(fastaPath, 'file')
    error(['A ComplexPortalDB.fasta already exists at: ' fastaPath])
else
    disp('Fetching sequences from UniProt database. This might take a while.\n')
    for i = 1:numel(keys) % this is going through each id and requests thes equences 
        currentID = keys{i};
        url = ['https://rest.uniprot.org/uniprotkb/accessions?accessions=' num2str(currentID) '&format=fasta'];
        options = weboptions('ContentType', 'binary', 'RequestMethod', 'get');
    
        try
            compressedFastaData = webread(url, options);
    
            %Check if response is empty
            if isempty(compressedFastaData)
                %warning('Empty response for ID: %s', currentID);
                missingIDs{end+1} = currentID; % Store the ID in the missing list
                continue;
            end
            tempGzipFile = [tempname, '.gz'];
            fid = fopen(tempGzipFile, 'w');
            fwrite(fid, compressedFastaData);
            fclose(fid);
    
            % Decompress the gzip file
            gunzip(tempGzipFile);
            decompressedFile = tempGzipFile(1:end-3);  % Remove .gz extension
    
            % Read the decompressed FASTA data
            fid = fopen(decompressedFile, 'r');
            fastaData = fread(fid, '*char')';
            fclose(fid);
    
            % Verify if FASTA data was actually retrieved
            if isempty(fastaData)
                %warning('No FASTA data found for ID: %s', currentID);
                missingIDs{end+1} = currentID; % Store the ID in the missing list
                delete(tempGzipFile);
                delete(decompressedFile);
                continue;
            end
    
            %disp(options)
    
            %disp('first few bites: \n')
            %disp(compressedFastaData(1:min(10, numel(compressedFastaData))))
    
            % Write the FASTA data to the output file
            fid = fopen(fastaPath, 'a');
            fprintf(fid, '%s', fastaData);
            fclose(fid);
    
            %increment the counter
            foundID = foundID + 1;
    
            % Clean up temporary files
            delete(tempGzipFile);
            delete(decompressedFile);
    
        catch ME
            disp(ME.identifier)
            disp(ME.stack)
            disp(ME.message)
            disp(currentID)
    
            % Store the ID in the missing list
            missingIDs{end+1} = currentID;
            %warning('Error fetching sequence from uniProt for %s: %s', currentID, ME.message);
        end
    end
    
    disp(['The ComplexPortalDB.fasta file is stored at: ' fastaPath]);
    disp(['Number of IDs processed and written to the FASTA file: ' num2str(foundID)]);
    % If there are missing IDs, display them
    if ~isempty(missingIDs)
        disp('Sequence is not found for the following ID(s): \n');
        disp(missingIDs);
        remove(complexPortalProtIDs, missingIDs)
    end
    % the output gets the list of IDs that could not fetch a sequence
    noUniProtSeq = missingIDs;
end
noUniProtSeq = missingIDs;
end

