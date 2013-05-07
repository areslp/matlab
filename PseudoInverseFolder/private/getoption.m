function res = getoption(options, field, default)
% Retreive the data from the options structure
% If the field does not exist, return the default value
field = lower(field);
if ~isfield(options, field) || isempty(options.(field))
    res = default;
else
    res = options.(field);
end
end