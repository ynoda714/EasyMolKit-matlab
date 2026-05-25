function result = pubchemFetch(identifier, options)
% pubchemFetch  Fetch extended compound data from PubChem via PubChemPy.
%
%   result = emk.db.pubchemFetch(identifier)
%   result = emk.db.pubchemFetch(identifier, NameSpace="smiles")
%   result = emk.db.pubchemFetch(702, NameSpace="cid")
%
%   Queries PubChem using the PubChemPy Python library and returns a MATLAB
%   struct with an extended property set beyond what searchPubchem() provides:
%   synonyms, InChI, InChIKey, XLogP, TPSA, HBD/HBA counts, and more.
%
%   Requires: emk.setup.installExtra("pubchempy") before first use.
%
%   Arguments:
%     identifier  - string | char | numeric scalar.
%                   For NameSpace="cid", a numeric CID is also accepted.
%     NameSpace   - string (default "name"). Accepted values:
%                   "name"     - compound name (IUPAC, trivial, synonym, etc.)
%                   "smiles"   - SMILES string (canonical or isomeric)
%                   "cid"      - PubChem Compound ID (string or numeric)
%                   "inchi"    - InChI string
%                   "inchikey" - InChIKey (27-character hash)
%                   "formula"  - molecular formula (returns first match)
%     MaxSynonyms - positive integer (default 10).
%                   Maximum number of synonym strings included in result.Synonyms.
%
%   Returns:
%     result  - struct with fields:
%       .CID               double  - PubChem Compound ID
%       .IUPACName         string  - IUPAC systematic name (may be "" if absent)
%       .MolecularFormula  string  - molecular formula (e.g. "C2H6O")
%       .MolecularWeight   double  - average molecular weight (g/mol)
%       .IsomericSMILES    string  - isomeric SMILES
%       .InChI             string  - IUPAC InChI string
%       .InChIKey          string  - InChIKey (hashed InChI)
%       .XLogP             double  - XLogP (NaN when PubChem has no value)
%       .TPSA              double  - topological polar surface area (Ang^2)
%       .HBondDonors       double  - H-bond donor count
%       .HBondAcceptors    double  - H-bond acceptor count
%       .RotatableBonds    double  - rotatable bond count
%       .HeavyAtomCount    double  - heavy atom count
%       .Charge            double  - formal charge
%       .Complexity        double  - Bertz CT molecular complexity (NaN if absent)
%       .Synonyms          string  - row vector of synonyms (up to MaxSynonyms)
%
%   Errors:
%     emk:db:pubchemFetch:invalidInput     - identifier type is not supported
%     emk:db:pubchemFetch:invalidNamespace - unknown NameSpace value
%     emk:db:pubchemFetch:libraryNotFound  - pubchempy is not installed
%     emk:db:pubchemFetch:notFound         - no compound found for identifier
%     emk:db:pubchemFetch:pythonError      - unexpected Python exception
%
%   Example:
%     s = emk.db.pubchemFetch("ethanol");
%     s = emk.db.pubchemFetch("CCO", NameSpace="smiles");
%     s = emk.db.pubchemFetch(702, NameSpace="cid");
%     fprintf("MW = %.3f g/mol\n", s.MolecularWeight);
%     disp(s.Synonyms(1:3));
%
%   See also: emk.db.searchPubchem, emk.mol.fromSmiles,
%             emk.setup.installExtra

    arguments
        identifier
        options.NameSpace    (1,1) string = "name"
        options.MaxSynonyms  (1,1) double {mustBeInteger, mustBePositive} = 10
    end

    % --- Input validation: identifier ---
    if isnumeric(identifier) && isscalar(identifier)
        identifier = string(identifier);          % numeric CID -> string
    elseif ischar(identifier) || isStringScalar(identifier)
        identifier = string(identifier);
    else
        error("emk:db:pubchemFetch:invalidInput", ...
            "identifier must be a string, char, or numeric scalar CID, got: %s", ...
            class(identifier));
    end

    if strlength(strtrim(identifier)) == 0
        error("emk:db:pubchemFetch:invalidInput", ...
            "identifier must not be empty or whitespace-only");
    end

    % --- NameSpace validation ---
    VALID_NS = ["name", "smiles", "cid", "inchi", "inchikey", "formula"];
    ns = lower(options.NameSpace);
    if ~ismember(ns, VALID_NS)
        error("emk:db:pubchemFetch:invalidNamespace", ...
            "NameSpace must be one of: %s. Got: %s", ...
            strjoin(VALID_NS, ", "), options.NameSpace);
    end

    % --- Check pubchempy availability ---
    try
        py.importlib.import_module("pubchempy");
    catch ME
        if contains(ME.message, "ModuleNotFoundError") || ...
                contains(ME.message, "No module named")
            error("emk:db:pubchemFetch:libraryNotFound", ...
                "PubChemPy is not installed in the active Python environment.\n" + ...
                "Install it with: emk.setup.installExtra(""pubchempy"")");
        end
        rethrow(ME);
    end

    % --- Call PubChemPy ---
    logInfo("pubchemFetch: querying PubChem for '%s' (NameSpace=%s) ...", ...
        identifier, ns);
    try
        pcp          = py.importlib.import_module("pubchempy");
        pyCompounds  = pcp.get_compounds(identifier, ns);
    catch ME
        error("emk:db:pubchemFetch:pythonError", ...
            "PubChemPy raised an exception: %s", ME.message);
    end

    compList = cell(pyCompounds);
    if isempty(compList)
        error("emk:db:pubchemFetch:notFound", ...
            "No compound found for identifier: '%s' (NameSpace=%s)", ...
            identifier, ns);
    end

    % Use first match (PubChem returns ranked results)
    c = compList{1};

    % --- Extract properties ---
    result.CID              = pyDbl_(c.cid);
    result.IUPACName        = pyStr_(c.iupac_name);
    result.MolecularFormula = pyStr_(c.molecular_formula);
    result.MolecularWeight  = pyStrDbl_(c.molecular_weight);  % PubChemPy returns str
    % PubChem API (2025+) no longer returns isomeric_smiles via the compound
    % object; it now exposes a unified "SMILES" property via get_properties.
    % Fall back to get_properties("SMILES") when isomeric_smiles is None.
    smilesVal = pyStr_(c.isomeric_smiles);
    if strlength(smilesVal) == 0
        try
            smileProps = cell(pcp.get_properties("SMILES", string(result.CID), "cid"));
            if ~isempty(smileProps)
                smilesVal = string(smileProps{1}.get("SMILES"));
                if isequal(smilesVal, "None"); smilesVal = ""; end
            end
        catch
            smilesVal = "";
        end
    end
    result.IsomericSMILES   = smilesVal;
    result.InChI            = pyStr_(c.inchi);
    result.InChIKey         = pyStr_(c.inchikey);
    result.XLogP            = pyDbl_(c.xlogp);
    result.TPSA             = pyDbl_(c.tpsa);
    result.HBondDonors      = pyDbl_(c.h_bond_donor_count);
    result.HBondAcceptors   = pyDbl_(c.h_bond_acceptor_count);
    result.RotatableBonds   = pyDbl_(c.rotatable_bond_count);
    result.HeavyAtomCount   = pyDbl_(c.heavy_atom_count);
    result.Charge           = pyDbl_(c.charge);
    result.Complexity       = pyDbl_(c.complexity);

    % Synonyms: take up to MaxSynonyms
    try
        synList         = cell(c.synonyms);
        nSyn            = min(numel(synList), options.MaxSynonyms);
        result.Synonyms = string(synList(1:nSyn));
    catch
        result.Synonyms = string.empty(1, 0);
    end

    logInfo("pubchemFetch: fetched CID=%d (%s).", result.CID, result.MolecularFormula);
end

% =========================================================================
% Private helpers
% =========================================================================

function s = pyStr_(val)
% Convert a Python str (or None) to a MATLAB string.
    try
        s = string(val);
        if isequal(s, "None")
            s = "";
        end
    catch
        s = "";
    end
end

function d = pyDbl_(val)
% Convert a Python numeric or None to a MATLAB double (NaN for None).
    try
        d = double(val);
    catch
        d = NaN;
    end
    if ~isfinite(d) && ~isnan(d)
        d = NaN;   % guard against unexpected Inf from Python
    end
end

function d = pyStrDbl_(val)
% Convert a Python str-encoded number (e.g. "46.07") to MATLAB double.
% PubChemPy returns molecular_weight as a string.
    try
        s = string(val);
        if isequal(s, "None") || strlength(s) == 0
            d = NaN;
        else
            d = str2double(s);
        end
    catch
        d = NaN;
    end
end
