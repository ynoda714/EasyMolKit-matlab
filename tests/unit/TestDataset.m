classdef TestDataset < matlab.unittest.TestCase
% TestDataset  Unit tests for emk.dataset.* (L07).
%
%   Covers: esol, freesolv, bbbp, tox21
%   These tests verify argument validation (RDKit-free) and, where
%   network access is available, the schema of downloaded tables.
%
%   Network tests are tagged with "network" and can be excluded:
%     suite = testsuite("tests/unit/TestDataset");
%     % To skip network tests: (no built-in filter in basic runner)
%     results = runner.run(suite);

    properties (Access = private)
        % Temporary cache directory for isolated test downloads
        tmpCacheDir (1,1) string = ""
    end

    methods (TestMethodSetup)
        function setupTmpCache(tc)
            tc.tmpCacheDir = fullfile(tempdir(), ...
                sprintf("emk_test_dataset_%s", char(java.util.UUID.randomUUID())));
            mkdir(tc.tmpCacheDir);
        end
    end

    methods (TestMethodTeardown)
        function cleanupTmpCache(tc)
            if isfolder(tc.tmpCacheDir)
                rmdir(tc.tmpCacheDir, "s");
            end
        end
    end

    % ======================================================================
    % TC-1: Input validation -- emk.dataset.esol (RDKit NOT required)
    % ======================================================================
    methods (Test)
        function test_esol_unknownOption_throwsInvalidInput(tc)
            tc.verifyError( ...
                @() emk.dataset.esol("Bogus", 1), ...
                "emk:dataset:esol:invalidInput");
        end

        function test_esol_oddArgs_throwsInvalidInput(tc)
            tc.verifyError( ...
                @() emk.dataset.esol("CacheDir"), ...
                "emk:dataset:esol:invalidInput");
        end

        function test_freesolv_unknownOption_throwsInvalidInput(tc)
            tc.verifyError( ...
                @() emk.dataset.freesolv("Bogus", 1), ...
                "emk:dataset:freesolv:invalidInput");
        end

        function test_bbbp_unknownOption_throwsInvalidInput(tc)
            tc.verifyError( ...
                @() emk.dataset.bbbp("Bogus", 1), ...
                "emk:dataset:bbbp:invalidInput");
        end

        function test_tox21_unknownOption_throwsInvalidInput(tc)
            tc.verifyError( ...
                @() emk.dataset.tox21("Bogus", 1), ...
                "emk:dataset:tox21:invalidInput");
        end
    end

    % ======================================================================
    % TC-2: Network tests -- ESOL download and schema
    % These tests require internet access; skip gracefully if unavailable.
    % ======================================================================
    methods (Test)
        function test_esol_returnsTable(tc)
            tbl = tc.downloadWithFallback( ...
                @() emk.dataset.esol("CacheDir", tc.tmpCacheDir));
            if isempty(tbl), return; end
            tc.verifyClass(tbl, "table", "esol must return a table");
        end

        function test_esol_hasRequiredColumns(tc)
            tbl = tc.downloadWithFallback( ...
                @() emk.dataset.esol("CacheDir", tc.tmpCacheDir));
            if isempty(tbl), return; end
            required = ["SMILES", "Name", "logS", "logS_Delaney", "MolWt"];
            for col = required
                tc.verifyTrue(ismember(col, tbl.Properties.VariableNames), ...
                    sprintf("esol table must have column '%s'", col));
            end
        end

        function test_esol_smilesIsString(tc)
            tbl = tc.downloadWithFallback( ...
                @() emk.dataset.esol("CacheDir", tc.tmpCacheDir));
            if isempty(tbl), return; end
            tc.verifyClass(tbl.SMILES, "string");
        end

        function test_esol_logSIsDouble(tc)
            tbl = tc.downloadWithFallback( ...
                @() emk.dataset.esol("CacheDir", tc.tmpCacheDir));
            if isempty(tbl), return; end
            tc.verifyClass(tbl.logS, "double");
        end

        function test_esol_rowCount_around1128(tc)
            tbl = tc.downloadWithFallback( ...
                @() emk.dataset.esol("CacheDir", tc.tmpCacheDir));
            if isempty(tbl), return; end
            % Dataset has 1128 molecules; allow small variation across versions
            tc.verifyGreaterThan(height(tbl), 1000, ...
                "ESOL should have > 1000 molecules");
            tc.verifyLessThan(height(tbl), 1300, ...
                "ESOL should have < 1300 molecules");
        end

        function test_esol_cache_prevents_redownload(tc)
            % Call twice; second call should use cache (no error)
            tbl1 = tc.downloadWithFallback( ...
                @() emk.dataset.esol("CacheDir", tc.tmpCacheDir));
            if isempty(tbl1), return; end
            tbl2 = emk.dataset.esol("CacheDir", tc.tmpCacheDir);
            tc.verifyEqual(height(tbl1), height(tbl2), ...
                "Cached and fresh tables must have same row count");
        end

        % --- BBBP schema ---
        function test_bbbp_hasRequiredColumns(tc)
            tbl = tc.downloadWithFallback( ...
                @() emk.dataset.bbbp("CacheDir", tc.tmpCacheDir));
            if isempty(tbl), return; end
            required = ["SMILES", "Name", "BBB"];
            for col = required
                tc.verifyTrue(ismember(col, tbl.Properties.VariableNames), ...
                    sprintf("bbbp table must have column '%s'", col));
            end
        end

        function test_bbbp_BBBIsLogical(tc)
            tbl = tc.downloadWithFallback( ...
                @() emk.dataset.bbbp("CacheDir", tc.tmpCacheDir));
            if isempty(tbl), return; end
            tc.verifyClass(tbl.BBB, "logical");
        end

        function test_bbbp_rowCount_around2039(tc)
            tbl = tc.downloadWithFallback( ...
                @() emk.dataset.bbbp("CacheDir", tc.tmpCacheDir));
            if isempty(tbl), return; end
            tc.verifyGreaterThan(height(tbl), 1800);
            tc.verifyLessThan(height(tbl), 2500);
        end

        % --- Tox21 schema ---
        function test_tox21_hasRequiredColumns(tc)
            tbl = tc.downloadWithFallback( ...
                @() emk.dataset.tox21("CacheDir", tc.tmpCacheDir));
            if isempty(tbl), return; end
            required = ["SMILES", "MolID", "NR_AR", "SR_ARE", "SR_p53"];
            for col = required
                tc.verifyTrue(ismember(col, tbl.Properties.VariableNames), ...
                    sprintf("tox21 table must have column '%s'", col));
            end
        end

        function test_tox21_rowCount_around7831(tc)
            tbl = tc.downloadWithFallback( ...
                @() emk.dataset.tox21("CacheDir", tc.tmpCacheDir));
            if isempty(tbl), return; end
            tc.verifyGreaterThan(height(tbl), 7000);
            tc.verifyLessThan(height(tbl), 9000);
        end

        function test_tox21_endpointValues_inRange(tc)
            % All endpoint values must be 0, 1, or NaN
            tbl = tc.downloadWithFallback( ...
                @() emk.dataset.tox21("CacheDir", tc.tmpCacheDir));
            if isempty(tbl), return; end
            endpoints = ["NR_AR", "SR_ARE", "SR_p53"];
            for col = endpoints
                vals = tbl.(col);
                validMask = ~isnan(vals);
                uniqueVals = unique(vals(validMask));
                tc.verifyTrue(all(ismember(uniqueVals, [0, 1])), ...
                    sprintf("Endpoint %s must contain only 0, 1, or NaN", col));
            end
        end
    end

    % ======================================================================
    % Helper: attempt download, skip test on network error
    % ======================================================================
    methods (Access = private)
        function tbl = downloadWithFallback(tc, fnHandle)
            try
                tbl = fnHandle();
            catch ME
                if contains(ME.identifier, "downloadFailed") || ...
                   contains(ME.message, "network") || ...
                   contains(ME.message, "Unable to")
                    tc.assumeFail("Network unavailable; skipping download test");
                end
                tbl = [];
                rethrow(ME);
            end
        end
    end
end
