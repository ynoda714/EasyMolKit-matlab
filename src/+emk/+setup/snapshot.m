function snap = snapshot()
% snapshot  Capture current environment versions for RF02 version lock.
%
%   snap = emk.setup.snapshot()
%
%   Returns a struct suitable for serialising to lock_snapshot.json:
%     snap.matlab     (string) - MATLAB release, e.g. "R2026a"
%     snap.python     (string) - Python version, e.g. "3.10.6"
%     snap.rdkit      (string) - RDKit version, e.g. "2024.03.6"
%     snap.commit     (string) - Git short hash or "n/a"
%     snap.toolboxes  (struct) - installed toolboxes as name -> version
%     snap.timestamp  (string) - ISO-8601 capture time
%
%   Errors: none (non-throwing by design; missing fields fall back to "unavailable"/"n/a").
%   See also: emk.setup.lockfile, emk.setup.verifyLock

    snap.matlab    = string(version("-release"));
    snap.python    = capturePython_();
    snap.rdkit     = captureRdkit_();
    snap.commit    = captureCommit_();
    snap.toolboxes = captureToolboxes_();
    snap.timestamp = char(datetime("now", "Format", "yyyy-MM-dd'T'HH:mm:ss"));

    logInfo("snapshot: MATLAB=%s  Python=%s  RDKit=%s  commit=%s", ...
        snap.matlab, snap.python, snap.rdkit, snap.commit);
end

% --------------------------------------------------------------------------

function v = capturePython_()
    try
        pe = pyenv();
        v  = string(pe.Version);
        if strlength(v) == 0; v = "unavailable"; end
    catch
        v = "unavailable";
    end
end

function v = captureRdkit_()
    try
        rdmod = py.importlib.import_module("rdkit");
        v     = string(py.getattr(rdmod, "__version__"));
    catch
        v = "unavailable";
    end
end

function h = captureCommit_()
    try
        [status, out] = system("git rev-parse --short HEAD");
        if status == 0
            h = string(strtrim(out));
        else
            h = "n/a";
        end
    catch
        h = "n/a";
    end
end

function tb = captureToolboxes_()
    tb = struct();
    try
        info = ver();
        for i = 1:numel(info)
            fname = matlab.lang.makeValidName(info(i).Name);
            tb.(fname) = string(info(i).Version);
        end
    catch
        % Non-critical -- return empty struct if ver() fails
    end
end
