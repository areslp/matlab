%EPS2PDF  Convert an eps file to pdf format using ghostscript
%
% Examples:
%   eps2pdf source dest
%   eps2pdf(source, dest, crop)
%
% This function converts an eps file to pdf format. If the output pdf file
% already exists, the eps file is appended as a new page on the end of the
% eps file.
%
% This function requires that you have ghostscript installed on your
% system. Ghostscript can be downloaded from: http://www.ghostscript.com
%
%IN:
%   source - filename of the source eps file to convert. The filename is
%            assumed to already have the extension ".eps".
%   dest - filename of the destination pdf file. The filename is assumed to
%          already have the extension ".pdf".
%   crop - boolean indicating whether to crop the borders off the pdf.
%          Default: true.

% Copyright (C) Oliver Woodford 2009

% Suggestion of appending pdf files provided by Matt C at:
% http://www.mathworks.com/matlabcentral/fileexchange/23629

% $Id: eps2pdf.m,v 1.3 2009/07/29 20:15:24 ojw Exp $

function eps2pdf(source, dest, crop)
% Set crop option
if nargin < 3 || ~isequal(crop, 0)
    options =  '-dEPSCrop ';
else
    options = '';
end
% Construct the options string for ghostscript
options = [options '-q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile="' dest '"'];
% Check if the output file exists
if exist(dest, 'file') == 2
    % File exists, so append current figure to the end
    tmp_nam = tempname;
    % Copy the file
    copyfile(dest, tmp_nam);
    % Add the output file names
    options = [options ' "' tmp_nam '" "' source '"'];
    try
        % Convert to pdf using ghostscript
        ghostscript(options);
    catch
        % Delete the intermediate file
        delete(tmp_nam);
        rethrow(lasterror);
    end
    % Delete the intermediate file
    delete(tmp_nam);
else
    % File doesn't exist
    % Add the output file names
    options = [options ' "' source '"'];
    % Convert to pdf using ghostscript
    ghostscript(options);
end
return

