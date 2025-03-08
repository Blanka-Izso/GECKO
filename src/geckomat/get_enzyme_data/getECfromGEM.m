function [model, invalidEC, invalidECpos] = getECfromGEM(model, ecRxns)
% getECfromGEM
%   Use the model.eccodes to populates the model.ec.eccodes field. EC
%   numbers that are not formatted as four numbers separated by periods,
%   possibly with trailing wildcards. Examples: 1.2.3.4 or 1.2.3.- while
%   invalid EC numbers are 1.2.3 or 1_2_3_4. Multiple EC numbers are
%   separated by ; for instance 1.2.3.4;1.2.3.5 not 1.2.3.4|1.2.3.5.
%
% Input:
%   model           an ecModel in GECKO 3 format (with ecModel.ec structure)
%   ecRxns          logical of length model.ec.rxns that specifies for
%                   which reactions the existing model.ec.eccodes entry
%                   should be kept and not modified by this function
%                   (optional, by default all model.ec.eccodes entries
%                   are populated by this function)
%
% Output:
%   model           ecModel with populated model.ec.eccodes
%   invalidEC       incorrectly formatted EC numbers
%   invalidECpos    position of invalidEC in model.eccodes
%
% Usage:
%   [model, invalidEC, invalidECpos] = getECfromGEM(model, ecRxns)

if ~isfield(model,'eccodes')
    error('The model has no model.eccodes field.')
end

%Need to remove the prefix of GECKO light rxn names in the ec structure
if ~model.ec.geckoLight
    rxnNames = model.ec.rxns;
else
    rxnNames = extractAfter(model.ec.rxns, 4);
end

[~,rxnIdxs] = ismember(rxnNames,model.rxns);


% Capture each EC number and validate their validity separately
eccodes = model.eccodes;

validECPattern = '^(\d+|-)\.(\d+|-)\.(\d+|-)\.(\d+|-)$';

invalidEC = cell(size(eccodes));
invalidECpos = [];

for i = 1:numel(eccodes)
    if isempty(eccodes{i})
        continue;
    end

    codes = strsplit(eccodes{i}, ';');

    validCodes = {};
    invalidCodes = {};

    for j = 1:numel(codes)
        if ~isempty(regexp(codes{j}, validECPattern, 'once'))
            validCodes{end+1} = codes{j};

        else
            invalidCodes{end+1} = codes{j};
        end
    end


    if ~isempty(validCodes)
        eccodes{i} = strjoin(validCodes, ';');
    else
        eccodes{i} = '';
    end

    if ~isempty(invalidCodes)
        invalidEC{i} = strjoin(invalidCodes, ';');
        invalidECpos = [invalidECpos; i];
    else
        invalidEC{i} = '';
    end

    invalidEC = invalidEC(~cellfun(@isempty, invalidEC));
    invalidEC = invalidEC(:);

    if nargin<2 || all(ecRxns)
        model.ec.eccodes = eccodes(rxnIdxs);
    else
        if ~isfield(model.ec,'eccodes')
            model.ec.eccodes(1:numel(model.ec.rxns),1) = {''};
        end
        model.ec.eccodes(ecRxns) = eccodes(rxnIdxs(ecRxns));
    end
end
if ~isempty(invalidECpos)
    if nargout < 2
        fprintf('Skipped incorrectly formatted EC numbers, rerun getECfromGEM with all outputs to get a list of invalid entries.\n');
    else
        fprintf('Skipped incorrectly formatted EC numbers.\n');
    end
end
end

%}
