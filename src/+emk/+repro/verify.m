function result = verify(metrics, criteria)
% verify  Check reproduction metrics against RF03 acceptance criteria.
%
%   result = emk.repro.verify(metrics, criteria)
%
%   Arguments:
%     metrics  : struct -- numeric metric values
%                  e.g. metrics.rmse_cv = 1.017; metrics.r2_cv = 0.762;
%     criteria : struct -- acceptance bounds per metric field.
%                  Each field is a sub-struct with "upper" and/or "lower":
%                    criteria.rmse_cv = struct("upper", 1.20);
%                    criteria.r2_cv   = struct("lower", 0.75);
%                  For classification:
%                    criteria.auc     = struct("lower", 0.80);
%                    criteria.accuracy = struct("lower", 0.75, "upper", 1.0);
%
%   Returns:
%     result.pass    (logical) - true if ALL criteria pass
%     result.details (struct)  - per metric: .value .pass .criteria
%     result.report  (string)  - formatted multiline summary for logging/saving
%
%   Example:
%     crit.rmse_cv = struct("upper", 1.20);
%     crit.r2_cv   = struct("lower", 0.75);
%     met.rmse_cv  = 1.017;
%     met.r2_cv    = 0.762;
%     r = emk.repro.verify(met, crit);
%     disp(r.report);
%
%   Errors: invalidInput, missingMetric
%   See also: emk.setup.lockfile, emk.setup.verifyLock

    if ~isstruct(metrics) || ~isstruct(criteria)
        error("emk:repro:verify:invalidInput", ...
            "Both 'metrics' and 'criteria' must be structs.");
    end

    fields  = fieldnames(criteria);
    allPass = true;
    lines   = {};

    for i = 1:numel(fields)
        f   = fields{i};
        cri = criteria.(f);

        if ~isfield(metrics, f)
            error("emk:repro:verify:missingMetric", ...
                "metrics.%s is required by criteria but missing.", f);
        end

        val    = metrics.(f);
        passed = true;
        parts  = {};

        if isfield(cri, 'upper')
            if val <= cri.upper
                parts{end+1} = sprintf('<= %.4g [OK]', cri.upper);
            else
                passed = false;
                parts{end+1} = sprintf('<= %.4g [FAIL got %.4g]', cri.upper, val);
            end
        end

        if isfield(cri, 'lower')
            if val >= cri.lower
                parts{end+1} = sprintf('>= %.4g [OK]', cri.lower);
            else
                passed = false;
                parts{end+1} = sprintf('>= %.4g [FAIL got %.4g]', cri.lower, val);
            end
        end

        if ~passed; allPass = false; end

        statusStr    = 'PASS'; if ~passed; statusStr = 'FAIL'; end
        boundsStr    = strjoin(parts, ' & ');
        lines{end+1} = sprintf('  %-20s = %+.4f  %s  %s', f, val, statusStr, boundsStr);

        result.details.(f) = struct('value', val, 'pass', passed, 'criteria', cri);
    end

    overallStr = 'PASS'; if ~allPass; overallStr = 'FAIL'; end
    allLines   = [{'RF03 Verification:'} ; lines(:) ; ...
                  {sprintf('  => Overall: %s', overallStr)}];
    result.pass   = allPass;
    result.report = strjoin(allLines, newline);

    if allPass
        logInfo("repro.verify: PASS");
    else
        logWarn("repro.verify: FAIL -- see result.report for details");
    end
end
