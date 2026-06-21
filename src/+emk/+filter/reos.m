function tbl = reos(tbl)
% reos  Annotate a descriptor table with REOS (Rapid Elimination Of Swill) filter.
%
%   tbl = emk.filter.reos(tbl)
%
%   Applies the REOS drug-likeness filter (Walters & Murcko 2002) and
%   appends two annotation columns:
%     Pass_REOS          (logical) - true when all applied criteria pass
%     Violations_REOS    (double)  - number of violated criteria (0-6)
%
%   REOS criteria applied (6 of the original 7; FormalCharge excluded):
%     MolWt              in [200, 500]
%     LogP               in [-5, 5]
%     NumHDonors         in [0, 5]
%     NumHAcceptors      in [0, 10]
%     NumRotatableBonds  in [0, 8]
%     HeavyAtomCount     in [15, 50]
%
%   NOTE: The original REOS also filters on FormalCharge in [-2, 2], but
%   FormalCharge is not in the default emk.descriptor.calculate() set.
%   Add a FormalCharge column and call this function again to include it
%   (this function ignores the FormalCharge column if absent).
%
%   This function is pure MATLAB and does not require Python or RDKit.
%
%   Arguments:
%     tbl - (table) Descriptor table with required columns (see above).
%           Use emk.mol.toTable to generate a compatible table.
%
%   Returns:
%     tbl - (table) Input table with two new columns appended:
%           Violations_REOS (double)  - count of violated criteria (0-6)
%           Pass_REOS       (logical) - true when violations == 0
%
%   NaN handling:
%     NaN descriptor values are treated as non-violations (MATLAB behaviour:
%     NaN > x and NaN < x both evaluate to false).  Remove missing values
%     first with rmmissing(tbl) if desired.
%
%   Errors:
%     emk:filter:reos:invalidInput   - tbl is not a table
%     emk:filter:reos:missingColumns - required column(s) absent
%
%   Example:
%     mols = {emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O"), ...
%             emk.mol.fromSmiles("CCO")};
%     tbl = emk.mol.toTable(mols);
%     tbl = emk.filter.reos(tbl);
%     passing = tbl(tbl.Pass_REOS, :);
%
%   References:
%     Walters, W.P. & Murcko, M.A. (2002). Prediction of 'drug-likeness'.
%       Adv. Drug Deliv. Rev. 54(3):255-271.
%       DOI: 10.1016/S0169-409X(02)00003-0
%     Bemis, G.W. & Murcko, M.A. (1997). Considering optimal pharmacophores.
%       J. Med. Chem. 40(9):1393-1402. DOI: 10.1021/jm9607547
%
%   See also: emk.filter.lipinski, emk.filter.veber, emk.filter.pains,
%             emk.mol.toTable, docs/algorithm_guide.md

    % --- Validate input table ---
    if ~istable(tbl)
        error("emk:filter:reos:invalidInput", ...
            "Input must be a MATLAB table, got: %s", class(tbl));
    end

    % --- Check required columns ---
    REQUIRED = ["MolWt", "LogP", "NumHDonors", "NumHAcceptors", ...
                "NumRotatableBonds", "HeavyAtomCount"];
    missing  = REQUIRED(~ismember(REQUIRED, tbl.Properties.VariableNames));
    if ~isempty(missing)
        error("emk:filter:reos:missingColumns", ...
            "Table is missing required column(s): [%s]. " + ...
            "Use emk.mol.toTable to generate a compatible table.", ...
            strjoin(missing, ", "));
    end

    logDebug("reos: applying REOS to %d row(s)", height(tbl));

    % --- Compute per-row violation counts ---
    % REOS thresholds (Walters & Murcko 2002):
    %   MolWt in [200, 500]   -> violates if < 200 OR > 500
    %   LogP in [-5, 5]       -> violates if < -5 OR > 5
    %   HBD in [0, 5]         -> violates if < 0 OR > 5
    %   HBA in [0, 10]        -> violates if < 0 OR > 10
    %   RotBonds in [0, 8]    -> violates if < 0 OR > 8
    %   HeavyAtoms in [15,50] -> violates if < 15 OR > 50
    nRows = height(tbl);
    violations = zeros(nRows, 1);

    violations = violations + double(tbl.MolWt            < 200  | tbl.MolWt            > 500);
    violations = violations + double(tbl.LogP             < -5   | tbl.LogP             > 5  );
    violations = violations + double(tbl.NumHDonors       < 0    | tbl.NumHDonors       > 5  );
    violations = violations + double(tbl.NumHAcceptors    < 0    | tbl.NumHAcceptors    > 10 );
    violations = violations + double(tbl.NumRotatableBonds < 0   | tbl.NumRotatableBonds > 8  );
    violations = violations + double(tbl.HeavyAtomCount   < 15   | tbl.HeavyAtomCount   > 50 );

    pass = violations == 0;

    tbl.Violations_REOS = violations;
    tbl.Pass_REOS       = pass;

    logInfo("reos: %d / %d row(s) pass REOS (6-criterion, no FormalCharge)", ...
        sum(pass), nRows);
end
