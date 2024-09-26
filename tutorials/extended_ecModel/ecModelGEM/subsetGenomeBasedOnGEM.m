function subsetGenomeBasedOnGEM(ecmodel, model, modelAdapter)
% subsetGenomeBasedOnGEM
%   Creates a fasta file subsetting the loaded genome
%   based on the metabolic genes present in the conventional GEM. 
%   The result is a 'uniprotGEM.fasta' file
%   containing the ORF IDs of the metabolic genes present in the conventional
%   GEM and their corresponging sequences. If a file named
%   'uniprotGEM.fasta' already exists the function will give an error
%   The file will be used for the BLAST algorithm
%
% Input: 
%       ecmodel:        extended model created by makeEcModel
%       model:          loaded conventional GEM
%       modelAdapter:   the default value is what is given in the adapter
% Output:
%       uniprotGEM.fasta: a .fasta file containing metabolic enzyme ORF 
%                       IDs and their sequence
% Usage:
%       subsetGenomeBasedOnGem(ecmodel, model, modelAdapter)

% if modelAdapter is not given, set it to default
if nargin < 3 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end

% Get the parameters from the modelAdapter
params      = modelAdapter.getParameters();
% Set the filePath to the given path given in the adapter
filePath    = fullfile(params.path,'data');

% Do the comparison and write the fasta file
% Checks whether a gene in the extended model is present in the
% conventional GEM writes the IDs and the corresponding sequences into a
% .fasta file
[Lia, ~] = ismember(ecmodel.ec.genes,model.genes);
GemGene = ecmodel.ec.genes(Lia);
GemSeq = ecmodel.ec.sequence(Lia);
%assert(length(GemGene) == length(GemSeq), 'Mismatch between genes and sequences.'); %Checking their length, should be equal

fastaPath = fullfile(filePath,'uniprotGEM.fasta');

if exist(fastaPath, 'file')
    error(['A uniprotGEM.fasta file already exists at: ' fastaPath])
else
    fid = fopen(fastaPath, 'w');
    for i = 1:length(GemSeq)
        fprintf(fid, '>%s\n%s\n', GemGene{i}, GemSeq{i});
    end
    fclose(fid);
    disp(['The uniprotGEM.fasta file is stored at: ', fastaPath])
end
end

