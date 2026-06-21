function tbl = veber(tbl)
% veber  Annotate a descriptor table with Veber oral bioavailability criteria.
%
%   tbl = emk.filter.veber(tbl)
%
%   Applies Veber's two rules for oral bioavailability (Veber et al. 2002)
%   and appends two annotation columns:
%     Pass_Veber          (logical) - true when both criteria pass
%     Violations_Veber    (double)  - number of violated criteria (0, 1, or 2)
%
%   Veber criteria (both must pass for oral bioavailability):
%     NumRotatableBonds  <= 10
%     TPSA               <= 140 Angstrom^2
%
%   This function is pure MATLAB and does not require Python or RDKit.
%   Required columns are provided by emk.mol.toTable with the default
%   Properties set (NumRotatableBonds and TPSA are in the default 10
%   descriptors).
%
%   Arguments:
%     tbl - (table) Descriptor table with columns NumRotatableBonds and TPSA.
%
%   Returns:
%     tbl - (table) Input table with two new columns appended:
%           Violations_Veber  (double)  - count of violated criteria (0-2)
%           Pass_Veber        (logical) - true when violations == 0
%           All original rows are preserved.
%           Use tbl(tbl.Pass_Veber, :) to extract passing rows.
%
%   NaN handling:
%     MATLAB evaluates NaN > threshold as false, so NaN descriptor values
%     are treated as non-violations.  Remove missing values first with
%     rmmissing(tbl) if NaN molecules should be excluded.
%
%   Errors:
%     emk:filter:veber:invalidInput    - tbl is not a table
%     emk:filter:veber:missingColumns  - required column(s) absent
%
%   Example:
%     mols = {emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O"), ...
%             emk.mol.fromSmiles("CCO")};
%     tbl  = emk.mol.toTable(mols);
%     tbl  = emk.filter.veber(tbl);
%     passing = tbl(tbl.Pass_Veber, :);
%
%   References:
%     Veber, D.F. et al. (2002). Molecular Properties That Influence the
%       Oral Bioavailability of Drug Candidates. J. Med. Chem.
%       45(12):2615-2623. DOI: 10.1021/jm020017n
%
%   See also: emk.filter.lipinski, emk.filter.reos, emk.mol.toTable,
%             docs/algorithm_guide.md

    % --- Validate input table ---
    if ~istable(tbl)
        error("emk:filter:veber:invalidInput", ...
            "Input must be a MATLAB table, got: %s", class(tbl));
    end

    % --- Check required columns ---
    REQUIRED = ["NumRotatableBonds", "TPSA"];
    missing  = REQUIRED(~ismember(REQUIRED, tbl.Properties.VariableNames));
    if ~isempty(missing)
        error("emk:filter:veber:missingColumns", ...
            "Table is missing required column(s): [%s]. " + ...
            "Use emk.mol.toTable to generate a compatible table.", ...
            strjoin(missing, ", "));
    end

    logDebug("veber: applying Veber criteria to %d row(s)", height(tbl));

    % --- Compute per-row violation counts ---
    % Veber (2002) criteria:
    %   NumRotatableBonds > 10  violates
    %   TPSA              > 140 violates
    nRows = height(tbl);
    violations = zeros(nRows, 1);
    violations = violations + double(tbl.NumRotatableBonds > 10);
    violations = violations + double(tbl.TPSA             > 140);

    pass = violations == 0;

    tbl.Violations_Veber = violations;
    tbl.Pass_Veber       = pass;

    logInfo("veber: %d / %d row(s) pass Veber criteria", sum(pass), nRows);
end
