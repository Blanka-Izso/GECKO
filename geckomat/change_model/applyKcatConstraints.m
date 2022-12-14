function model = applyKcatConstraints(model,updateRxns)
% applyKcatConstraints
%   Applies kcat-derived enzyme constraints to an ec-model. Existing enzyme
%   constraints are first removed (unless updateRxns is provided), and new
%   constraints are defined based on the content of model.ec.kcat.
%
% Input:
%   model       ec-model that was generated by makeEcModel, or loaded from
%               an earlier run. Not compatible with ec-models generated by
%               earlier GECKO versions (pre 3.0).
%   updateRxns  if not all enzyme constraints should be updated, this can
%               be given as either a logical vector of length
%               model.ec.rxns, a vector of model.ec.rxns indices, or a
%               (cell array of) string(s) with model.ec.rxns identifiers.
%
% Output:
%   model       ec-model where reactions are constrained by enzyme usage
%               if a kcat value was provided for the reaction-enzyme pair
%               in model.ec.kcat
%
% Usage: model = applyKcatConstraints(model,updateRxns);

if nargin<2
    updateRxns = true(numel(model.ec.rxns),1);
elseif isnumeric(updateRxns)
    updateRxnsLog = false(numel(model.ec.rxns),1);
    updateRxnsLog(updateRxns) = true;
    updateRxns = updateRxnsLog;
elseif iscellstr(updateRxns) || ischar(updateRxns) || isstring(updateRxns)
    updateRxnsIds = convertCharArray(updateRxns);
    updateRxns = ismember(model.ec.rxns,updateRxnsIds);
end
    
if ~isfield(model,'ec')
    error(['No model.ec structure could be found: the provided model is'...
           ' not a valid GECKO3 ec-model. First run makeEcModel(model).'])
end
if all(model.ec.kcat==0)
    warning('No kcat values are provided in model.ec.kcat, model remains unchanged.')
    return
end

%Clear existing incorporation of enzyme usage
% TODO: modify for geckoLight 
%if ~model.ec.geckoLight
protMetIdx = startsWith(model.mets,'prot_') & ~strcmp(model.mets,'prot_pool');
metabolRxn = unique(model.ec.rxns(updateRxns));
metabolRxn = ismember(model.rxns,metabolRxn);
model.S(protMetIdx,metabolRxn) = 0;

if ~model.ec.geckoLight %For normal GECKO formulation, where each enzyme is explicitly considered
    %Make vector where column 1 = rxn; 2 = enzyme; 3 = subunit copies; 4 = kcat
    newKcats=zeros(numel(updateRxns)*10,4);
    updateRxns=find(updateRxns);
    kcatFirst=0;
    for i=1:numel(updateRxns)
        enzymes   = find(model.ec.rxnEnzMat(i,:));
        kcatLast  = kcatFirst+numel(enzymes);
        kcatFirst = kcatFirst+1;
        newKcats(kcatFirst:kcatLast,1) = updateRxns(i);
        newKcats(kcatFirst:kcatLast,2) = enzymes;
        newKcats(kcatFirst:kcatLast,3) = model.ec.rxnEnzMat(i,enzymes);
        newKcats(kcatFirst:kcatLast,4) = model.ec.kcat(i);
        kcatFirst = kcatLast;
    end
    newKcats(kcatLast+1:end,:)=[];
    
    newKcats(:,4) = newKcats(:,4) * 3600; %per second -> per hour
    newKcats(:,4) = newKcats(:,4).^-1; %Inverse: hours per reaction
    newKcats(:,4) = newKcats(:,4)*1000; %In umol instead of mmol
    newKcats(:,4) = newKcats(:,3).*newKcats(:,4); %Multicopy subunits.
    %Unit of enzyme usage is umol/gDW/h, while metabolic flux is in mmol/gDW/h.
    %This prevents very low fluxes.
    
    %Replace rxns and enzymes with their location in model
    [~,newKcats(:,1)] = ismember(model.ec.rxns(newKcats(:,1)),model.rxns);
    [~,newKcats(:,2)] = ismember(strcat('prot_',model.ec.enzymes(newKcats(:,2))),model.mets);
    linearIndices     = sub2ind(size(model.S),newKcats(:,2),newKcats(:,1));
    model.S(linearIndices) = -newKcats(:,4); %Substrate = negative
    
else %GECKO light formulation, where prot_pool represents all usages
    newKcats=zeros(numel(updateRxns),1);
    for i=1:numel(updateRxns)
        %TODO:The current code is incorrect, it takes the sum of all
        %enzymes that are associated. But there could be isoenzymes, so not
        %all enzymes are necessarily required. Not sure how this can easily
        %be parsed (without running expandModel).
        enzymes     = find(model.ec.rxnEnzMat(i,:));
        MW          = sum(model.ec.mw(enzymes));
        newKcats(i) = MW/model.ec.kcat(i);
        %TODO:Correct units.
    end
    %All draw directly from protein pool
    prot_pool_idx = find(ismember(model.mets,'prot_pool'));

    [~,rxnIds] = ismember(model.ec.rxns(updateRxns),model.rxns);
    model.S(prot_pool_idx,rxnIds) = -newKcats; %Substrate = negative
end
end
