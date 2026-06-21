function varargout = lockfile(varargin)
% lockfile  Save or load an RF02 version lock JSON file.
%
%   Save mode (2 args):
%     emk.setup.lockfile(snap, filePath)
%
%   Load mode (1 string arg):
%     snap = emk.setup.lockfile(filePath)
%
%   Arguments (save mode):
%     snap     - struct returned by emk.setup.snapshot()
%     filePath - output path for the JSON file (parent dir must exist)
%
%   Arguments (load mode):
%     filePath - path to an existing lock JSON file
%
%   The JSON is written with PrettyPrint=true for readability and version control.
%
%   Errors: invalidInput, writeError, fileNotFound
%   See also: emk.setup.snapshot, emk.setup.verifyLock

    if nargin == 2
        % Save mode: lockfile(snap, filePath)
        snap     = varargin{1};
        filePath = string(varargin{2});
        if ~isstruct(snap)
            error("emk:setup:lockfile:invalidInput", ...
                "First argument must be a struct (from emk.setup.snapshot).");
        end
        saveLock_(snap, filePath);
        if nargout > 0
            varargout{1} = snap;
        end

    elseif nargin == 1 && (ischar(varargin{1}) || isStringScalar(varargin{1}))
        % Load mode: snap = lockfile(filePath)
        filePath     = string(varargin{1});
        varargout{1} = loadLock_(filePath);

    else
        error("emk:setup:lockfile:invalidInput", ...
            "Usage: lockfile(snap, filePath) to save, " + ...
            "or snap = lockfile(filePath) to load.");
    end
end

% --------------------------------------------------------------------------

function saveLock_(snap, filePath)
    filePath = char(filePath);
    fid      = fopen(filePath, 'w');
    if fid == -1
        error('emk:setup:lockfile:writeError', ...
            'Cannot open file for writing: %s', filePath);
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '%s\n', jsonencode(snap, 'PrettyPrint', true));
    logInfo('lockfile: saved -> %s', filePath);
end

function snap = loadLock_(filePath)
    if ~isfile(char(filePath))
        error("emk:setup:lockfile:fileNotFound", ...
            "Lock file not found: %s", filePath);
    end
    raw  = fileread(char(filePath));
    snap = jsondecode(raw);
    logInfo("lockfile: loaded <- %s", filePath);
end
