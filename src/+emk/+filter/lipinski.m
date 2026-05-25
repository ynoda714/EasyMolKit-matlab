function tbl = lipinski(tbl, varargin)
% lipinski  Annotate a descriptor table with Lipinski Rule of Five (Ro5) results.
%
%   tbl = emk.filter.lipinski(tbl)
%   tbl = emk.filter.lipinski(tbl, MaxViolations=1)
%
%   Applies Lipinski's Rule of Five to each row of a descriptor table and
%   appends two annotation columns:
%     Pass_Ro5      (logical) - true when violations <= MaxViolations
%     Violations_Ro5 (double) - number of violated Ro5 criteria (0-4)
%
%   The four Ro5 criteria (Lipinski et al. 1997):
%     MW   (MolWt)         <= 500 Da
%     LogP                 <= 5
%     HBD  (NumHDonors)    <= 5
%     HBA  (NumHAcceptors) <= 10
%
%   The function is pure MATLAB and does not require Python or RDKit.
%   Required columns are provided by emk.mol.toTable with the default
%   Properties set (MolWt, LogP, NumHDonors, NumHAcceptors are all included
%   in the default 10-descriptor output).
%
%   Arguments:
%     tbl           - (table) Descriptor table with columns MolWt, LogP,
%                     NumHDonors, NumHAcceptors.
%     MaxViolations - (double, name-value, default=0) Maximum number of
%                     allowed Ro5 violations for Pass_Ro5 to be true.
%                     Range [0, 4]. Common choices:
%                       0 - strict Ro5 (all 4 criteria must pass)
%                       1 - relaxed Ro5 (allows one violation; independent
%                           of Veber et al. 2002, which uses TPSA/RotBonds)
%
%   Returns:
%     tbl - (table) Input table with two new columns appended at the right:
%           Pass_Ro5       (logical) - true if violations <= MaxViolations
%           Violations_Ro5 (double)  - count of violated criteria, 0-4
%           All original rows are preserved (no rows removed).
%           Use tbl(tbl.Pass_Ro5, :) to extract passing rows.
%
%   NaN handling:
%     MATLAB evaluates NaN > threshold as false, so a NaN descriptor value
%     is treated as a NON-violation for that criterion.  A row containing
%     NaN in MolWt, LogP, NumHDonors, or NumHAcceptors will NOT accumulate
%     a violation for that criterion, and may appear to pass Ro5.
%     WARNING: molecules with unknown (NaN) descriptors may spuriously
%     receive Pass_Ro5 = true.  Remove them first with rmmissing(tbl) if
%     unknown values should be excluded from the passing set.
%
%   Validation order:
%     1. tbl must be a MATLAB table  (invalidInput)
%     2. MaxViolations must be an integer in [0,4]  (invalidMaxViol)
%     3. Required columns must be present  (missingColumns)
%
%   Errors:
%     emk:filter:lipinski:invalidInput    - tbl is not a table, or bad arg
%     emk:filter:lipinski:missingColumns  - required column(s) absent
%     emk:filter:lipinski:invalidMaxViol  - MaxViolations not integer in [0,4]
%
%   Example:
%     mols = {emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O"), ...  % aspirin  MW=180
%             emk.mol.fromSmiles("c1ccc2ccccc2c1")};             % naphthalene MW=128
%     tbl  = emk.mol.toTable(mols);
%     tbl  = emk.filter.lipinski(tbl);
%     passing = tbl(tbl.Pass_Ro5, :);
%
%     % Relaxed filter (allow 1 violation):
%     tbl = emk.filter.lipinski(tbl, MaxViolations=1);
%
%   See also: emk.mol.toTable, emk.descriptor.calculate,
%             docs/algorithm_guide.md section 10

    % --- Parse optional name-value arguments ---
    maxViol = 0;
    if nargin > 1
        if mod(numel(varargin), 2) ~= 0
            error("emk:filter:lipinski:invalidInput", ...
                "Optional arguments must be name-value pairs.");
        end
        for k = 1:2:numel(varargin)
            argName = varargin{k};
            if ~(ischar(argName) || isStringScalar(argName))
                error("emk:filter:lipinski:invalidInput", ...
                    "Argument name must be a string scalar, got: %s", class(argName));
            end
            if strcmpi(argName, "MaxViolations")
                maxViol = varargin{k+1};
            else
                error("emk:filter:lipinski:invalidInput", ...
                    "Unknown argument: '%s'. Supported: 'MaxViolations'.", argName);
            end
        end
    end

    % --- Validate input table FIRST (primary argument) ---
    % Table check comes before option validation so that the most
    % user-visible argument (tbl) produces actionable errors first.
    if ~istable(tbl)
        error("emk:filter:lipinski:invalidInput", ...
            "Input must be a MATLAB table, got: %s", class(tbl));
    end

    % --- Validate MaxViolations ---
    if ~(isnumeric(maxViol) && isscalar(maxViol) && ...
         isfinite(maxViol)  && maxViol >= 0 && maxViol <= 4 && ...
         floor(maxViol) == maxViol)
        error("emk:filter:lipinski:invalidMaxViol", ...
            "MaxViolations must be an integer in [0, 4], got: %s", mat2str(maxViol));
    end

    % --- Check required columns ---
    REQUIRED = ["MolWt", "LogP", "NumHDonors", "NumHAcceptors"];
    missing  = REQUIRED(~ismember(REQUIRED, tbl.Properties.VariableNames));
    if ~isempty(missing)
        error("emk:filter:lipinski:missingColumns", ...
            "Table is missing required column(s): [%s]. " + ...
            "Use emk.mol.toTable to generate a compatible table.", ...
            strjoin(missing, ", "));
    end

    logDebug("lipinski: applying Ro5 to %d row(s) (MaxViolations=%d)", ...
        height(tbl), maxViol);

    % --- Compute per-row violation counts ---
    % Each criterion contributes 1 to the count when violated.
    % Ro5 thresholds (Lipinski et al. 1997):
    %   MW    > 500   violates
    %   LogP  > 5     violates
    %   HBD   > 5     violates
    %   HBA   > 10    violates
    nRows = height(tbl);
    violations = zeros(nRows, 1);
    violations = violations + double(tbl.MolWt         > 500);
    violations = violations + double(tbl.LogP          > 5  );
    violations = violations + double(tbl.NumHDonors    > 5  );
    violations = violations + double(tbl.NumHAcceptors > 10 );

    pass = violations <= maxViol;

    tbl.Violations_Ro5 = violations;
    tbl.Pass_Ro5       = pass;

    logInfo("lipinski: %d / %d row(s) pass Ro5 (MaxViolations=%d)", ...
        sum(pass), nRows, maxViol);
end
