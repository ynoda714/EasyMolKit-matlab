function result = verifyLock(lockRef)
% verifyLock  Verify current environment against a saved RF02 lock.
%
%   result = emk.setup.verifyLock(filePath)    -- load lock from JSON file
%   result = emk.setup.verifyLock(lockStruct)  -- compare against struct directly
%
%   Critical fields (mismatch sets result.pass = false):
%     matlab, python, rdkit
%
%   Non-critical fields (mismatch appended to result.warnings only):
%     commit
%
%   Returns:
%     result.pass      (logical) - true if all critical fields match
%     result.details   (struct)  - per critical field: .expected .actual .match
%     result.warnings  (string)  - pipe-separated non-critical mismatches ("" if none)
%
%   Errors: invalidInput
%   See also: emk.setup.snapshot, emk.setup.lockfile

    if ischar(lockRef) || isStringScalar(lockRef)
        locked = emk.setup.lockfile(lockRef);
    elseif isstruct(lockRef)
        locked = lockRef;
    else
        error("emk:setup:verifyLock:invalidInput", ...
            "Input must be a file path (string/char) or a lock struct.");
    end

    current = emk.setup.snapshot();

    criticals = ["matlab", "python", "rdkit"];
    noncrit   = ["commit"];

    result.pass     = true;
    result.details  = struct();
    result.warnings = "";

    for i = 1:numel(criticals)
        f   = criticals(i);
        exp = getField_(locked,  f);
        act = getField_(current, f);
        ok  = strcmp(exp, act);

        result.details.(f) = struct("expected", exp, "actual", act, "match", ok);

        if ok
            logInfo("verifyLock: OK       [%s] %s", f, act);
        else
            result.pass = false;
            logWarn("verifyLock: MISMATCH [%s] expected=%s  actual=%s", f, exp, act);
        end
    end

    for i = 1:numel(noncrit)
        f   = noncrit(i);
        exp = getField_(locked,  f);
        act = getField_(current, f);
        if ~strcmp(exp, act)
            logWarn("verifyLock: WARN [%s] expected=%s  actual=%s (non-critical)", f, exp, act);
            if strlength(result.warnings) == 0
                result.warnings = f;
            else
                result.warnings = result.warnings + " | " + f;
            end
        end
    end

    if result.pass
        logInfo("verifyLock: PASS -- environment matches locked snapshot.");
    else
        logWarn("verifyLock: FAIL -- environment differs from locked snapshot.");
    end
end

% --------------------------------------------------------------------------

function v = getField_(s, fieldName)
    fname = char(fieldName);
    if isfield(s, fname)
        v = string(s.(fname));
    else
        v = "n/a";
    end
end
