%EXPORT_FIG  Exports figures suitable for publication
%
% Examples:
%   im = export_fig
%   [im alpha] = export_fig
%   export_fig filename
%   export_fig filename -format1 -format2
%   export_fig ... -nocrop
%   export_fig ... -a2
%   export_fig ... -zbuffer
%   export_fig(..., handle)
%
% This function saves a figure or single axes to one or more vector and/or
% bitmap file formats, and/or outputs a rasterized version to the
% workspace, with the following properties:
%   - Figure/axes reproduced as it appears on screen
%   - Cropped borders
%   - Embedded fonts (vector formats)
%   - Improved line and grid line styles (vector formats)
%   - Anti-aliased graphics (bitmap formats)
%   - Transparent background supported (pdf, eps, png)
%   - Semi-transparent patch objects supported (png only)
%   - Append to file (pdf only)
%   - Vector formats: pdf, eps
%   - Bitmap formats: png, tif, jpg, bmp, export to workspace 
%   
% This function is especially suited to exporting figures for use in
% publications and presentations, because of the high quality and
% portability of media produced.
%
% Note that the background color and figure dimensions are reproduced
% (the latter approximately, and ignoring cropping) in the output file. For
% transparent background (and semi-transparent patch objects), set the
% figure (and axes) 'Color' property to 'none'; pdf, eps and png are the
% only file formats to support a transparent background, whilst the png
% format alone supports transparency of patch objects. 
%
% When exporting to vector format (pdf & eps), this function requires that
% ghostscript is installed on your system. You can download this from:
%   http://www.ghostscript.com
% When exporting to eps it additionally requires pdftops, from the Xpdf
% suite of functions. You can download this from:
%   http://www.foolabs.com/xpdf
%
%IN:
%   filename - string containing the name (optionally including full or
%              relative path) of the file the figure is to be saved as. If
%              a path is not specified, the figure is saved in the current
%              directory. If no name and no output arguments are specified,
%              the default name, 'export_fig_out', is used. If neither a
%              file extension nor a format are specified, a ".png" is added
%              and the figure saved in that format. If the file already
%              exists, it is overwritten, except for pdf files, for which
%              the figure is appended as a new page.
%   -format1, -format2, etc. - strings containing the extensions of the
%                              file formats the figure is to be saved as.
%                              Valid options are: '-pdf', '-eps', '-png',
%                              '-tif', '-jpg' and '-bmp'. All combinations
%                              of formats are valid.
%   -nocrop - optional string indicating the borders of the output are not
%             to be cropped.
%   -a1, -a2, -a3, -a4 - optional sting indicating the amount of
%                        anti-aliasing to use for bitmap outputs. -a1 means
%                        no anti-aliasing, -a4 is the maximum amount
%                        (default).
%   -zbuffer - optional string indicating the zbuffer renderer should be
%              used in place of the opengl renderer (default) for bitmaps.
%              Note that this renderer doesn't support transparent patches.
%   handle - The handle of the figure or axes to be saved. Default: gcf.
%
%OUT:
%   im - MxNxC uint8 image array of the figure.
%   alpha - MxN single array of alphamatte values in range [0,1], for the
%           case when the background is transparent.
%
%   See also PRINT, SAVEAS.

% Copyright (C) Oliver Woodford 2008-2009

% The idea of using ghostscript is inspired by Peder Axensten's SAVEFIG
% (fex id: 10889) which is itself inspired by EPS2PDF (fex id: 5782).
% The idea for using pdftops came from the MATLAB newsgroup (id: 168171).
% The idea of editing the EPS file to change line styles comes from Jiro
% Doke's FIXPSLINESTYLE (fex id: 17928).
% The idea of changing dash length with line width came from comments on
% fex id: 5743, but the implementation is mine :)
% The idea of anti-aliasing bitmaps came from Anders Brun's MYAA (fex id:
% 20979).
% The idea of appending figures in pdfs came from Matt C in comments on the
% FEX (id: 23629)

% $Id: export_fig.m,v 1.20 2009/07/31 20:22:13 ojw Exp $

function [im alpha] = export_fig(varargin)
% Parse the input arguments
[fig options] = parse_args(nargout, varargin{:});
% Isolate the subplot, if it is one
cls = strcmp(get(fig, 'Type'), 'axes');
if cls
    % Given a handle of a single set of axes
    fig = isolate_subplot(fig);
else
    old_mode = get(fig, 'InvertHardcopy');
end
% Hack the font units where necessary (due to a font rendering bug in
% print?). This may not work perfectly in all cases. Also it can change the
% figure layout if reverted, so use a copy.
if isbitmap(options) && options.aa_factor > 1
    fontu = findobj(fig, 'FontUnits', 'normalized');
    if ~isempty(fontu)
        % Some normalized font units found
        if ~cls
            fig = copyfig(fig);
            set(fig, 'Visible', 'off');
            fontu = findobj(fig, 'FontUnits', 'normalized');
            cls = true;
        end
        set(fontu, 'FontUnits', 'points');
    end
end
% Set to print exactly what is there
set(fig, 'InvertHardcopy', 'off');
% Do the bitmap formats first
if isbitmap(options)
    % Get the background colour
    tcol = get(fig, 'Color');
    if isequal(tcol, 'none') && (options.png || options.alpha)
        % Get out an alpha channel
        % Set the background colour to black
        set(fig, 'Color', 'k');
        % Print large version to array
        B = print2array(fig, options.aa_factor, options.zbuffer);
        % Downscale the image
        B = downsize(single(B), 0, options.aa_factor);
        % Set background to white
        set(fig, 'Color', 'w');
        % Print large version to array
        A = print2array(fig, options.aa_factor, options.zbuffer);
        % Downscale the image
        A = downsize(single(A), 255, options.aa_factor);
        % Set the background colour back to normal
        set(fig, 'Color', 'none');
        % Compute the alpha map
        alpha = sum(B - A, 3) / (255*3) + 1;
        A = alpha;
        A(A==0) = 1;
        A = uint8(B ./ A(:,:,[1 1 1]));
        clear B
        % Crop the background
        if options.crop
            [alpha v] = crop_background(alpha, 0);
            A = A(v(1):v(2),v(3):v(4),:);
        end
        if options.png
            % Save the png
            imwrite(A, [options.name '.png'], 'Alpha', alpha);
            % Clear the png bit
            options.png = false;
        end
        % Return only one channel for greyscale
        if isbitmap(options)
            A = check_greyscale(A);
        end
        if options.alpha
            % Store the image
            im = A;
            % Clear the alpha bit
            options.alpha = false;
        end
        % Get the non-alpha image
        if isbitmap(options)
            alph = alpha(:,:,ones(1, size(A, 3)));
            A = uint8(single(A) .* alph + 255 * (1 - alph));
            clear alph
        end
        if options.im
            % Store the new image
            im = A;
        end
    else
        % Print large version to array
        if isequal(tcol, 'none')
            set(fig, 'Color', 'w');
            A = print2array(fig, options.aa_factor, options.zbuffer);
            set(fig, 'Color', 'none');
            tcol = 255;
        else
            A = print2array(fig, options.aa_factor, options.zbuffer);
            tcol = tcol * 255;
            if ~isequal(tcol, round(tcol))
                tcol = squeeze(A(1,1,:));
            end
        end
        % Crop the background
        if options.crop
            A = crop_background(A, tcol);
        end
        % Downscale the image
        A = downsize(A, tcol, options.aa_factor);
        % Return only one channel for greyscale
        A = check_greyscale(A);
        % Outputs
        if options.im
            im = A;
        end
        if options.alpha
            im = A;
            alpha = zeros(size(A, 1), size(A, 2), 'single');
        end
    end
    % Save the images
    for a = {'png', 'tif', 'bmp'}
        if options.(a{1})
            imwrite(A, [options.name '.' a{1}]);
        end
    end
    % Save jpeg with higher quality than default
    if options.jpg
        imwrite(A, [options.name '.jpg'], 'Quality', 95);
    end
end
% Now do the vector formats
if isvector(options)
    % Generate some filenames
    tmp_nam = [tempname '.eps'];
    if options.pdf
        pdf_nam = [options.name '.pdf'];
    else
        pdf_nam = [tempname '.pdf'];
    end
    try
        % Generate an eps
        print2eps(tmp_nam, fig);
        % Generate a pdf
        eps2pdf(tmp_nam, pdf_nam, options.crop);
    catch
        % Delete the eps
        delete(tmp_nam);
        rethrow(lasterror);
    end
    % Delete the eps
    delete(tmp_nam);
    if options.eps
        try
            % Generate an eps from the pdf
            pdf2eps(pdf_nam, [options.name '.eps']);
        catch
            if ~options.pdf
                % Delete the pdf
                delete(pdf_nam);
            end
            rethrow(lasterror);
        end
        if ~options.pdf
            % Delete the pdf
            delete(pdf_nam);
        end
    end
end
if cls
    % Close the created figure
    close(fig);
else
    % Reset the hardcopy mode
    set(fig, 'InvertHardcopy', old_mode);
end
return

function [fig options] = parse_args(nout, varargin)
% Parse the input arguments
% Set the defaults
fig = get(0, 'CurrentFigure');
options = struct('name', 'export_fig_out', ...
                 'crop', true, ...
                 'zbuffer', false, ...
                 'pdf', false, ...
                 'eps', false, ...
                 'png', false, ...
                 'tif', false, ...
                 'jpg', false, ...
                 'bmp', false, ...
                 'im', nout == 1, ...
                 'alpha', nout == 2, ...
                 'aa_factor', 4);

% Go through the other arguments
for a = 1:nargin-1
    if ishandle(varargin{a})
        fig = varargin{a};
    elseif ischar(varargin{a}) && ~isempty(varargin{a})
        if varargin{a}(1) == '-'
            switch lower(varargin{a}(2:end))
                case 'nocrop'
                    options.crop = false;
                case 'zbuffer'
                    options.zbuffer = true;
                case 'pdf'
                    options.pdf = true;
                case 'eps'
                    options.eps = true;
                case 'png'
                    options.png = true;
                case {'tif', 'tiff'}
                    options.jpg = true;
                case {'jpg', 'jpeg'}
                    options.jpg = true;
                case 'bmp'
                    options.bmp = true;
                case {'a1', 'a2', 'a3', 'a4'}
                    options.aa_factor = str2double(varargin{a}(3));
            end
        else
            name = varargin{a};
            if numel(name) > 3 && name(end-3) == '.' && any(strcmpi(name(end-2:end), {'pdf', 'eps', 'png', 'tif', 'jpg', 'bmp'}))
                options.(lower(name(end-2:end))) = true;
                name = name(1:end-4);
            end
            options.name = name;
        end
    end
end

% Set the default format
if ~isvector(options) && ~isbitmap(options)
    options.png = true;
end
return

function fh = isolate_subplot(ah, vis)
% Isolate the axes in a figure on their own
% Tag the axes so we can find them in the copy
old_tag = get(ah, 'Tag');
set(ah, 'Tag', 'AxesToCopy');
% Create a new figure exactly the same as the old one
fh = copyfig(ancestor(ah, 'figure')); %copyobj(ancestor(ah, 'figure'), 0);
if nargin < 2 || ~vis
    set(fh, 'Visible', 'off');
end
% Reset the axes tag
set(ah, 'Tag', old_tag);
% Get all the axes
axs = findobj(fh, 'Type', 'axes');
% Find the axes to save
ah = findobj(axs, 'Tag', 'AxesToCopy');
if numel(ah) ~= 1
    close(fh);
    error('Too many axes found');
end
I = true(size(axs));
I(axs==ah) = false;
% Set the axes tag to what it should be
set(ah, 'Tag', old_tag);
% Keep any legends which overlap the subplot
ax_pos = get(ah, 'OuterPosition');
ax_pos(3:4) = ax_pos(3:4) + ax_pos(1:2);
for ah = findobj(axs, 'Tag', 'legend', '-or', 'Tag', 'Colorbar')'
    leg_pos = get(ah, 'OuterPosition');
    leg_pos(3:4) = leg_pos(3:4) + leg_pos(1:2);
    % Overlap test
    if leg_pos(1) < ax_pos(3) && leg_pos(2) < ax_pos(4) &&...
       leg_pos(3) > ax_pos(1) && leg_pos(4) > ax_pos(2)
        I(axs==ah) = false;
    end
end
% Delete all axes except for the input axes and associated items
delete(axs(I));
return

function fh = copyfig(fh)
% Is there a legend?
if isempty(findobj(fh, 'Type', 'axes', 'Tag', 'legend'))
    % Safe to copy using copyobj
    fh = copyobj(fh, 0);
else
    % copyobj will change the figure, so save and then load it instead
    tmp_nam = [tempname '.fig'];
    hgsave(fh, tmp_nam);
    fh = hgload(tmp_nam);
    delete(tmp_nam);
end
return

function A = downsize(A, padval, factor)
% Downsample an image
if factor == 1
    % Nothing to do
    return
end
try
    % Faster, but requires image processing toolbox
    A = imresize(A, 1/factor, 'bilinear');
catch
    % No image processing toolbox - resize manually
    % Lowpass filter - use Gaussian as is separable, so faster
    switch factor
        case 4
            % sigma: 1.7
            filt = single([0.0148395 0.0498173 0.118323 0.198829 0.236384 0.198829 0.118323 0.0498173 0.0148395]);
        case 3
            % sigma: 1.35
            filt = single([0.025219 0.099418 0.226417 0.297892 0.226417 0.099418 0.025219]);
        case 2
            % sigma: 1.0
            filt = single([0.054489 0.244201 0.40262 0.244201 0.054489]);
    end
    padding = floor(numel(filt) / 2);
    if numel(padval) == 3 && padval(1) == padval(2) && padval(2) == padval(3)
        padval = padval(1);
    end
    if numel(padval) == 1
        B = repmat(single(padval), [size(A, 1) size(A, 2)] + (2 * padding));
    end
    for a = 1:size(A, 3)
        if numel(padval) == 3
            B = repmat(single(padval(a)), [size(A, 1) size(A, 2)] + (2 * padding));
        end
        B(padding+1:end-padding,padding+1:end-padding) = A(:,:,a);
        A(:,:,a) = conv2(filt, filt', B, 'valid');
    end
    clear B
    % Subsample
    A = A(2:factor:end,2:factor:end,:);
end
return

function A = check_greyscale(A)
% Check if the image is greyscale
if size(A, 3) == 3 && ...
        all(reshape(A(:,:,1) == A(:,:,2), [], 1)) && ...
        all(reshape(A(:,:,2) == A(:,:,3), [], 1))
    A = A(:,:,1); % Save only one channel for 8-bit output
end
return

function [A v] = crop_background(A, bcol)
% Map the foreground pixels
M = A(:,:,1) ~= bcol(1);
for a = 2:size(A, 3)
    M = M | A(:,:,a) ~= bcol(min(a, end));
end
% Crop the background
N = any(M, 1);
M = any(M, 2);
v = [find(M, 1) find(M, 1, 'last') find(N, 1) find(N, 1, 'last')];
A = A(v(1):v(2),v(3):v(4),:);
return

function b = isvector(options)
b = options.pdf || options.eps;
return

function b = isbitmap(options)
b = options.png || options.tif || options.jpg || options.bmp || options.im || options.alpha;
return