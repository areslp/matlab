function success = cvx_grbgetkey( kcode, overwrite )

% CVX_GRBGETKEY   Retrieves and saves a Gurobi/CVX license.
%
% This function is used to install Gurobi license keys for use in CVX. It
% is called with your Gurobi license code as a string argument; e.g.
%
%     cvx_grbgetkey xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
%     cvx_grbgetkey( 'xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx' )
% 
% If you have not yet obtained a Gurobi license code, please visit the page
% 
%     <a href="matlab: web('http://www.gurobi.com/documentation/5.5/quick-start-guide/node5','-browser');">http://www.gurobi.com/documentation/5.5/quick-start-guide/node5</a>
% 
% for information on your various options (trial, academic, commercial).
% Once you have received notice that your license has been created, visit
% the Gurobi license page
% 
%     <a href="matlab: web('http://www.gurobi.com/download/licenses/current','-browser');">http://www.gurobi.com/download/licenses/current</a>
% 
% to retrieve the 36-character license code.
%
% The retrieved Gurobi license will be stored in your MATLAB preferences
% directory (see the PREFDIR command). If a license file already exists at,
% this location, and its expiration date has not yet passed, CVX_GRBGETKEY
% will refuse to retrieve a new license. To override this behavior, call
% the function with an -OVERWRITE argument:
%
%     cvx_grbgetkey xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx -overwrite
%     cvx_grbgetkey( 'xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx', '-overwrite' )
% 
% Note that this utility is meant to retrieve licenses for the version of
% Gurobi that is bundled with CVX. While it will retrieve full licenses as
% well, it is strongly recommended that you move such licenses to one of the
% standard Gurobi locations, discussed here:
%
%    <a href="matlab: web(''http://www.gurobi.com/documentation/5.5/quick-start-guide/node8'',''-browser'');">http://www.gurobi.com/documentation/5.5/quick-start-guide/node8</a>
%
% Note that in order to use Gurobi with CVX, *both* a Gurobi license and a
% CVX Professional license are required; see this page for details:
%
%
%

success = true;
downloaded = false;
line = '---------------------------------------------------------------------------';
jprintf({
    ''
    line
    'CVX/Gurobi license key installer'
    line
});

%%%%%%%%%%%%%%%%%%%%%%%%%
% Process the arguments %
%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 2 && ischar( kcode ) && size( kcode, 1 ) == 1 && kcode(1) == '-',
    tmp = kcode;
    kcode = overwrite;
    overwrite = tmp;
end
emsg = [];
if nargin < 1,
    emsg = '';
elseif isempty( kcode ) || ~ischar( kcode ) || ndims( kcode ) > 2 || size( kcode, 1 ) ~= 1, %#ok
    emsg = 'Invalid license code: must be a string.';
elseif ~regexp( kcode, '^[0-9a-f]{8,8}-[0-9a-f]{4,4}-[0-9a-f]{4,4}-[0-9a-f]{4,4}-[0-9a-f]{12,12}$', 'once' )
    emsg = sprintf( 'Invalid license code: %s', kcode );
elseif nargin < 2 || isempty( overwrite ),
    overwrite = false;
elseif isequal( overwrite, '-overwrite' ),
    overwrite = true;
elseif ischar( overwrite ) && size( overwrite, 1 ) == 1,
    emsg = sprintf( 'Invalid argument: %s', overwrite );
else
    emsg = 'Invalid second argument.';
end
if ischar( emsg ),
    if ~isempty( emsg ),
        fprintf( '*** %s\n\n', emsg );
    end
    jprintf({
        'This function is used to install Gurobi license keys for use in CVX. It'
        'is called with your license code as an argument; e.g.'
        ''
        '    %s xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx'
        '    %s( ''xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx'' )'
        ''
        'If you have not yet obtained a license code, please visit the page'
        ''
        '    <a href="matlab: web(''http://www.gurobi.com/documentation/5.5/quick-start-guide/node5'',''-browser'');">http://www.gurobi.com/documentation/5.5/quick-start-guide/node5</a>'
        ''
        'for information on your various options (trial, academic, commercial). Once'
        'a license has been created, you may retrieve its 36-character code by'
        'logging into your Gurobi account and visiting the page'
        ''
        '    <a href="matlab: web(''http://www.gurobi.com/download/licenses/current'',''-browser'');">http://www.gurobi.com/download/licenses/current</a>'
        ''
    }, mfilename, mfilename );
    success = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check to see if Gurobi is present %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if success,
    [ fs, ps, mpath, mext ] = cvx_version; %#ok
    if fs == '\', fsre = '\\'; else fsre = fs; end
    gname = [ mpath, fs, 'gurobi' ];
    if ~exist( gname, 'dir' ),
        jprintf({
            'This function is meant to be used only with the version of Gurobi that is'
            'bundled with CVX; but your CVX installation does not include Gurobi.'
            'To rectify this, you may either download a CVX/Gurobi bundle from'
            ''
            '    <a href="matlab: web(''http://cvxr.com/cvx/download'',''-browser'');">http://cvxr.com/cvx/download</a>'
            ''
            'or download the full Gurobi package from Gurobi directly:'
            ''
            '    <a href="matlab: web(''http://www.gurobi.com'',''-browser'');">http://www.gurobi.com</a>'
            ''
            'In either case, you will need both a CVX Professional license and a Gurobi'
            'license to proceed.'
            });
        success = false;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ensure that this platform is compatible with the bundled Gurobi solver %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if success,
    gname = [ gname, fs, mext(4:end) ];
    if ~exist( gname, 'dir' ),
        mismatch = false;
        switch mext,
            case 'mexmaci', pform = '32-bit OSX';
            case 'glx',     pform = '32-bit Linux';
            case 'a64',     pform = '64-bit Linux'; mismatch = true;
            case 'maci64',  pform = '64-bit OSX'; mismatch = true;
            case 'w32',     pform = '32-bit Windows'; mismatch = true;
            case 'w64',     pform = '64-bit Windows'; mismatch = true;
            otherwise,      pform = [];
        end
        if mismatch,
            jprintf({
                'The %s version of Gurobi is missing, perhaps because you downloaded a CVX'
                'package for a different MATLAB platform. Please visit'
                ''
                '    <a href="matlab: web(''http://cvxr.com/cvx/download'',''-browser'');">http://cvxr.com/cvx/download</a>'
                ''
                'and download and install the correct package.'
                }, pform );
        elseif isempty( pform ),
            fprintf( 'CVX/Gurobi is not supported on the %s platform.\n', computer );
        else
            fprintf( 'CVX/Gurobi is not supported on %s. Consider using the 64-bit version of MATLAB.\n', pform );
        end
        success = false;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Confirm the existence of the grbgetkey utility %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if success,
    gname = [ gname, fs, 'grbgetkey' ];
	if fs == '\',
		gname = [ gname, '.exe' ];
	end
	if ~exist( gname, 'file' ),
		jprintf({
			'Your CVX package is missing the file'
			''
			'    %s'
			''
			'which is necessary to complete this task. Please reinstall CVX.'
            }, gname );
        success = false;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the preferences directory, if necessary %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if success,
	fdir = regexprep( prefdir, [ 'matlab', fsre, 'R\d\d\d\d\w$' ], 'matlab' );
	if ~exist( fdir, 'dir' ), %#ok
		[success,msg] = mkdir( fdir );
		if ~success, %#ok
			jprintf({
				'This function needs to write to the directory'
				''
				'     %s'
				''
				'which does not exist. An attempt to create it resulted in this error:'
				''
				'    %s'
				''
				'Please rectify this problem and try again.'
                }, fdir, msg );
            success = false;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Guard against overwriting an existing license %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if success,
    % We actually use the parent directory of 'prefdir' if that directory has
    % the standard version format: e.g., if prefdir is ~/joe/.matlab/R2012b,
    % we store our data in ~/joe/.matlab, so that all versions of MATLAB can
    % see the same CVX data.
    fname = [ fdir, fs, 'cvx_gurobi.lic' ];
    if exist( fname, 'file' ),
        if overwrite,
            msg = [];
        else
            msg = 'This license may or may not be current.';
            fid = fopen( fname );
            if fid ~= 0,
                fstr = fread( fid, Inf, 'uint8=>char' )';
                fclose( fid );
                matches = regexp( fstr, 'EXPIRATION=\d\d\d\d-\d\d-\d\d', 'once' );
                if matches && floor(datenum(fstr(matches+11:matches+20),'yyyy-mm-dd'))<floor(datenum(clock))
                    msg = [];
                else
                    msg = 'This license has not yet expired.';
                end
            end
        end
        if ~isempty( msg )
            jprintf({
                'An existing license file has been found at location:'
                ''
                '    %s'
                ''
                '%s If you wish to overwrite this file,'
                'please re-run this function with an "-overwrite" argument: that is,'
                ''
                '   %s %s -overwrite'
                '   %s( ''%s'', ''-overwrite'' )'
                ''
                'Otherwise, please move the existing license and try again.'
                }, fname, msg, mfilename, kcode, mfilename, kcode );
            success = false;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a temporary destination directory %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tdir = [];
if success,
    for k = 1 : 10,
        tdir = tempname;
        if ~exist( tdir, 'file' ), break; end
    end
    [success,msg] = mkdir(tdir);
    if ~success,
        jprintf({
            'This function attempted to create the temporary directory'
            ''
            '    %s'
            ''
            'but the following error occurred:'
            ''
            '   %s'
            ''
            'Please rectify the problem and try again.';
        }, tdir, msg );
        success = false;
        tdir = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
% Download the license %
%%%%%%%%%%%%%%%%%%%%%%%%

if success
    result = '';
	fprintf( 'Contacting the Gurobi Optimization license server...' );
	[status,result]=system( sprintf( '%s --path=%s %s', gname, tdir, kcode ) );
	fprintf( 'done.\n' );
    if any( strfind( result, 'Unable to determine hostname' ) ) || any( strfind( result, 'not recognized as belonging to an academic domain' ) ),
        jprintf({
            'The attempt to retrieve the license key failed with the following error'
            'while trying to verify your academic license eligibility:'
            ''
            '    %s'
            'For information about this error, please consult the Gurobi documentation'
            ''
            '    <a href="matlab: web(''http://www.gurobi.com/documentation/5.5/quick-start-guide/node5'',''-browser'');">http://www.gurobi.com/documentation/5.5/quick-start-guide/node5</a>'
            ''
            'Once you have rectified the error, please try again.'
        }, regexprep(result,'.*---------------------\n+(.*?)\n+$','$1') );
        success = false;
    elseif any( strfind( result, 'already issued for host' ) )
        matches = regexp( fstr, 'already issued for host ''([^'']+)''', 'once' );
        if ~isempty( matches ),
            matches = sprintf( ' (%s)', matches{1}{1} );
        else
            matches = '';
        end
        jprintf({
            'This license has already been issued for a different host%s.'
            'Please acquire a new license for this host from Gurobi.'
            }, matches );
        success = false;
    elseif ~any( strfind( result, 'License key saved to file' ) ),
        jprintf({
            'The attempt to retrieve the license key failed with the following error:'
            ''
            '    %s'
            'For information about this error, please consult the Gurobi documentation'
            ''
            '    <a href="matlab: web(''http://www.gurobi.com/documentation/5.5/quick-start-guide/node5'',''-browser'');">http://www.gurobi.com/documentation/5.5/quick-start-guide/node5</a>'
            ''
            'Once you have rectified the error, please try again.'
        }, regexprep(result,'.*---------------------\n+(.*?)\n+$','$1') );
        fprintf( '%s\n', line );
        fprintf( result(1+(result(1)==10):end) );
        success = false;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the license file to determine its expiration %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if success,
    tname = [ tdir, fs, 'gurobi.lic' ];
    [fid,msg] = fopen( tname, 'r' );
    if fid == 0,
        jprintf({
            'An unexpected error occured: the Gurobi license file could not be'
            'opened. The file was expected to be in the temporary location'
            ''
            '    %s'
            ''
            'but the attempt to read it resulted in the following error:'
            ''
            '    %s'
            ''
            'Please attempt to rectify the problem and try again; if necessary,'
            'please contact <a href="mailto:cvx@cvxr.com">CVX support</a>.'
            }, tname, msg );
        success = false;
        tdir = [];
    else
        fstr = fread( fid, Inf, 'uint8=>char' )';
        fclose( fid );
        matches = regexp( fstr, 'EXPIRATION=\d\d\d\d-\d\d-\d\d', 'once' );
        if matches,
            if floor(datenum(fstr(matches+11:matches+20),'yyyy-mm-dd'))<floor(datenum(clock))
                jprintf({
                    'The license was successfully downloaded, but it expired on %s.'
                    'Please contact Gurobi for a new license.'
                    }, fstr(matches+11:matches+20) );
                success = false;
            else
                jprintf({
                    'Download successful. The license can be found at'
                    ''
                    '    %s'
                    ''
                    'The license will expire at the end of the day on %s.'
                }, fname, fstr(matches+11:matches+20) ); 
            end
        end
        if success,
            matches = regexp( fstr, '\nAPPNAME=CVX\n', 'once' );
            if ~any( matches ),
                [ matches, mends ] = regexp( fstr, '\nAPPNAME=\w+', 'start', 'end' );
                if any( matches ),
                    jprintf({
                        line
                        'ERROR: This license is reserved for the application "%s" and will'
                        '*not* work with CVX. We strongly recommend that you move this license to'
                        'another location, in accordance with this documentation:'
                        ''
                        '    <a href="matlab: web(''http://www.gurobi.com/documentation/5.5/quick-start-guide/node8'',''-browser'');">http://www.gurobi.com/documentation/5.5/quick-start-guide/node8</a>'
                        ''            
                        'To use Gurobi with CVX, you must either obtain a full Gurobi license or a'
                        'CVX specific license.'
                    }, fstr(matches(1)+9:mends) );
                    success = false;
                else
                    jprintf({
                        line
                        'WARNING: This license is not a CVX-specific license. It will still work with'
                        'CVX; however, if you also wish to use it with other applications, then we'
                        'strongly recommend that you copy it to a standard Gurobi location. Consult'
                        'the following Gurobi documentation for more information:'
                        ''
                        '    <a href="matlab: web(''http://www.gurobi.com/documentation/5.5/quick-start-guide/node8'',''-browser'');">http://www.gurobi.com/documentation/5.5/quick-start-guide/node8</a>'
                        ''
                        'For instance, to move this file to the standard home directory location,'
                        'copy and paste this command into your MATLAB command line:'
                        ''
                        '    movefile %s ~/gurobi.lic'
                    }, fname );
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Move the license to its proper location %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if success,
    [success,msg] = movefile( tname, fname, 'f' );
    if ~success,
        jprintf({
            'The attempt to move the Gurobi license file from its temporary location'
            ''
            '    %s'
            ''
            'to its final location'
            ''
            '    %s'
            ''
            'resulted in the following error:'
            ''
            '    %s'
            ''
            'Please rectify the problem and try again; or move the file manually. Once'
            'the license is in place, please run CVX_SETUP to configure CVX to use the'
            'Gurobi solver with the new license.'
        }, tname, fname, msg );
        success = false;
    end
end

%%%%%%%%%%%%
% Clean up %
%%%%%%%%%%%%

if ~isempty(tdir),
    [success,message] = rmdir( tdir, 's' ); %#ok
end

if success,
    jprintf({
        line
'Now that the license has been retrieved, please run CVX_SETUP to configure'
'CVX to use the Gurobi solver with the new license.'
    });
end

jprintf({
    line
    ''
});

if nargout == 0,
    clear success
end

function jprintf( strs, varargin )
strs = strs(:)';
[ strs{2,:} ] = deal( '\n' );
fprintf( cat(2,strs{:}), varargin{:} );

% Copyright 2013 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
