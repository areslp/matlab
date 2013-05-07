function s = lowercase(s)
% convert all field of a structure to lower case
field = fieldnames(s);
for k=1:length(field)
    f = field{k};
    if ~strcmp(lower(f),f)
        s.(lower(f)) = s.(f);
        s = rmfield(s, f);
    end
end
end % lowercase