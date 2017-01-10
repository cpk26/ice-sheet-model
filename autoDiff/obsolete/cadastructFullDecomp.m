function [ outStruct ] = cadastructFullDecomp( cadaObj )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

outStruct = {};

tmpStruct = cadastructDecomp(cadaObj);

% Get the field names of the structure.
fields = fieldnames(tmpStruct, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);

for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('outStruct.%s = tmpStruct.%s.func.value;', thisField, thisField);
    eval(commandLine);
end


end

