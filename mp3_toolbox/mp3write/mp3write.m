function mp3write(varargin)
%MP3WRITE Write MP3 (".mp3") sound file.
%    MP3WRITE(Y,FS,NBITS,MP3FILE,ENCODING) writes data Y to a MP3
%    file specified by the file name MP3FILE, with a sample rate
%    of FS Hz and with NBITS number of bits. Stereo data should 
%    be specified as a matrix with two columns. 
%    ENCODING must be specified as an integer number from 1 to 5
% 
%   1 = Fixed bit rate 128kbs encoding.
%   2 = Fixed bit rate jstereo 128kbs encoding, high quality (recommended).
%   3 = Average bit rate 112kbs encoding.
%   4 = Fast encode, low quality.
%   5 = Variable bitrate.
%
%   Y,FS and NBITS are mandatory fields. If MP3FILE is not defined the file
%   name will be 'Default_name.mp3'. If ENCODING is not defined encoding
%   type '2' will be used by deault.
%
%    See also MP3READ, WAVREAD, WAVWRITE.
if length(varargin) < 3 | length(varargin) > 5
    error('Unsopported number of argument inputs') 
end
Y = varargin{1};
FS = varargin{2};
NBITS = varargin{3};
if NBITS~=8 & NBITS~=16 & NBITS~=24 & NBITS~=32
    error('Unsopported bit depth')
end
if length(varargin) >= 4
    MP3FILE = varargin{4};
    if ischar(MP3FILE) ~= 1
    error('File name is not a string') 
    end
else
    MP3FILE = 'Default_name.mp3';
    disp('File name = Default_name.mp3')
end

if isempty(findstr(MP3FILE,'.mp3'))
    MP3FILE = strcat(MP3FILE,'.mp3');
end

if length(varargin) == 5
    ENCODING = varargin{5};
else
    ENCODING = '2';
    disp('Fixed bit rate, joint-stereo, 128 kb/s encoding')
end

s = which('mp3write.m');
ww = findstr('mp3write.m',s);
lame = s(1:ww-2);
wavwrite(Y,FS,NBITS,strcat(lame,'\temp.wav'));
tmpfile = strcat(lame,'\temp.wav');
MP3FILE = strcat(pwd,'\',MP3FILE);
ENCODING =  num2str(ENCODING);
switch ENCODING
    case {'1'}
        cmd = [lame,'\lame', ' --quiet', ' ', tmpfile, ' ',MP3FILE];
    case {'2'}
        cmd = [lame,'\lame', ' --quiet', ' -b 128 ', tmpfile, ' ',MP3FILE];
    case {'3'}
        cmd = [lame,'\lame', ' --quiet', ' --abr 112 ', tmpfile, ' ',MP3FILE];
    case {'4'}
        cmd = [lame,'\lame', ' --quiet', ' -f ', tmpfile, ' ',MP3FILE];
    case {'5'}
        cmd = [lame,'\lame', ' --quiet', ' -h ', ' -V ', tmpfile, ' ',MP3FILE];
    otherwise
        error('Encoding parameters not suported') 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Data Encoding  using "Lame.exe"%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dos(cmd);
% Delete temporary file
delete(tmpfile);