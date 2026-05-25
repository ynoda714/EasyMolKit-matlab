classdef TestLoadConfig < matlab.unittest.TestCase
% TestLoadConfig  Unit tests for src/config/loadConfig.m
%
% Run with:
%   addpath(genpath("src"));
%   results = run(TestLoadConfig);

    methods (TestMethodSetup)
        function setupPath(tc) %#ok<MANU>
            addpath(genpath("src"));
        end
    end

    methods (Test)
        % ------------------------------------------------------------------
        function test_defaults_returnStruct(tc)
        % loadConfig returns a struct with expected top-level fields
            cfg = loadConfig();
            tc.verifyTrue(isstruct(cfg), "cfg must be a struct");
            tc.verifyTrue(isfield(cfg, "python"),  "missing field: python");
            tc.verifyTrue(isfield(cfg, "rdkit"),   "missing field: rdkit");
            tc.verifyTrue(isfield(cfg, "runtime"), "missing field: runtime");
            tc.verifyTrue(isfield(cfg, "output"),  "missing field: output");
            tc.verifyTrue(isfield(cfg, "run"),     "missing field: run");
        end

        % ------------------------------------------------------------------
        function test_defaults_pythonVersion(tc)
        % Default Python version is 3.10
            setenv("EMK_PYTHON_VERSION", "");
            cfg = loadConfig();
            tc.verifyEqual(cfg.python.version, "3.10");
        end

        % ------------------------------------------------------------------
        function test_defaults_rdkitVersion(tc)
        % Default RDKit version is set and non-empty string.
        % verifyClass confirms the conversion codepath returns a MATLAB string
        % (not a char or numeric), mirroring the string-priority policy.
            cfg = loadConfig();
            tc.verifyClass(cfg.rdkit.version, "string", ...
                "rdkit.version must be a MATLAB string");
            tc.verifyNotEmpty(cfg.rdkit.version, ...
                "rdkit.version must be non-empty");
        end

        % ------------------------------------------------------------------
        function test_envVar_overridesPythonVersion(tc)
        % EMK_PYTHON_VERSION env var overrides python.version
            setenv("EMK_PYTHON_VERSION", "3.11");
            tc.addTeardown(@() setenv("EMK_PYTHON_VERSION", ""));
            cfg = loadConfig();
            tc.verifyEqual(cfg.python.version, "3.11");
        end

        % ------------------------------------------------------------------
        function test_envVar_overridesEvalMode(tc)
        % EMK_RUNTIME_EVAL_MODE env var overrides runtime.eval_mode
            setenv("EMK_RUNTIME_EVAL_MODE", "ref");
            tc.addTeardown(@() setenv("EMK_RUNTIME_EVAL_MODE", ""));
            cfg = loadConfig();
            tc.verifyEqual(cfg.runtime.eval_mode, "ref");
        end

        % ------------------------------------------------------------------
        function test_envVar_booleanTrue(tc)
        % Boolean env var "true" maps to logical true.
        % verifyClass guards against numeric 1 being returned instead of logical.
            setenv("EMK_RDKIT_AUTO_INIT", "true");
            tc.addTeardown(@() setenv("EMK_RDKIT_AUTO_INIT", ""));
            cfg = loadConfig();
            tc.verifyClass(cfg.rdkit.auto_init, "logical", ...
                "auto_init must be logical (not double 1)");
            tc.verifyTrue(cfg.rdkit.auto_init, ...
                "auto_init must be true for EMK_RDKIT_AUTO_INIT=true");
        end

        % ------------------------------------------------------------------
        function test_envVar_boolean1(tc)
        % Boolean env var "1" maps to logical true.
        % verifyClass guards against "1" being stored as string or double.
            setenv("EMK_RDKIT_AUTO_INIT", "1");
            tc.addTeardown(@() setenv("EMK_RDKIT_AUTO_INIT", ""));
            cfg = loadConfig();
            tc.verifyClass(cfg.rdkit.auto_init, "logical", ...
                "auto_init must be logical (not double 1)");
            tc.verifyTrue(cfg.rdkit.auto_init, ...
                "auto_init must be true for EMK_RDKIT_AUTO_INIT=1");
        end

        % ------------------------------------------------------------------
        function test_envVar_booleanFalse(tc)
        % Boolean env var "false" maps to logical false.
        % verifyClass guards against 0 being returned as double instead of logical.
            setenv("EMK_RDKIT_AUTO_INIT", "false");
            tc.addTeardown(@() setenv("EMK_RDKIT_AUTO_INIT", ""));
            cfg = loadConfig();
            tc.verifyClass(cfg.rdkit.auto_init, "logical", ...
                "auto_init must be logical (not double 0)");
            tc.verifyFalse(cfg.rdkit.auto_init, ...
                "auto_init must be false for EMK_RDKIT_AUTO_INIT=false");
        end

        % ------------------------------------------------------------------
        function test_jsonFile_overridesDefault(tc)
        % settings.json values override built-in defaults
        % Write a minimal temp JSON and point CWD to its parent dir.
            tmpDir = fullfile(tempdir, "emk_test_cfg_" + ...
                char(datetime("now", "Format", "yyyyMMdd_HHmmssSSS")));
            mkdir(tmpDir);
            cfgDir = fullfile(tmpDir, "config");
            mkdir(cfgDir);
            tc.addTeardown(@() rmdir(tmpDir, "s"));

            json = '{"python":{"version":"3.9"},"rdkit":{"version":"2023.03.1","auto_init":true},"runtime":{"eval_mode":"prod","text_only_mode":true},"output":{"root_dir":"result/intermediate"},"run":{"root_dir":"result/runs","publish_latest":true}}';
            fid = fopen(fullfile(cfgDir, "settings.json"), "w");
            fprintf(fid, "%s", json);
            fclose(fid);

            prevDir = pwd;
            tc.addTeardown(@() cd(prevDir));
            cd(tmpDir);
            addpath(genpath(fullfile(prevDir, "src")));

            cfg = loadConfig();
            tc.verifyEqual(cfg.python.version, "3.9");
            tc.verifyEqual(cfg.rdkit.version, "2023.03.1");
        end

        % ------------------------------------------------------------------
        function test_jsonFile_unknownKeysIgnored(tc)
        % JSON keys absent from defaults (e.g., _comment) are silently ignored
            tmpDir = fullfile(tempdir, "emk_test_cfg_" + ...
                char(datetime("now", "Format", "yyyyMMdd_HHmmssSSS")));
            mkdir(tmpDir);
            cfgDir = fullfile(tmpDir, "config");
            mkdir(cfgDir);
            tc.addTeardown(@() rmdir(tmpDir, "s"));

            json = '{"_comment":"ignore me","python":{"version":"3.10"},"rdkit":{"version":"2024.03.6","auto_init":true},"runtime":{"eval_mode":"prod","text_only_mode":true},"output":{"root_dir":"result/intermediate"},"run":{"root_dir":"result/runs","publish_latest":true}}';
            fid = fopen(fullfile(cfgDir, "settings.json"), "w");
            fprintf(fid, "%s", json);
            fclose(fid);

            prevDir = pwd;
            tc.addTeardown(@() cd(prevDir));
            cd(tmpDir);
            addpath(genpath(fullfile(prevDir, "src")));

            cfg = loadConfig();
            tc.verifyFalse(isfield(cfg, "_comment"), "_comment key must not appear in cfg");
            tc.verifyEqual(cfg.python.version, "3.10");
        end

        % ------------------------------------------------------------------
        function test_defaults_extraLibraries_fieldExists(tc)
        % cfg must have a top-level extraLibraries struct (Track 1 versions).
            cfg = loadConfig();
            tc.verifyTrue(isfield(cfg, "extraLibraries"), ...
                "cfg must have field 'extraLibraries'");
            tc.verifyTrue(isstruct(cfg.extraLibraries), ...
                "cfg.extraLibraries must be a struct");
        end

        function test_defaults_extraLibraries_pubchempy_isString(tc)
        % cfg.extraLibraries.pubchempy must be a non-empty string (pinned version).
            cfg = loadConfig();
            tc.verifyTrue(isfield(cfg.extraLibraries, "pubchempy"), ...
                "extraLibraries must have field 'pubchempy'");
            tc.verifyClass(cfg.extraLibraries.pubchempy, "string", ...
                "pubchempy version must be a MATLAB string");
            tc.verifyGreaterThan(strlength(cfg.extraLibraries.pubchempy), 0, ...
                "pubchempy version must be non-empty");
        end

        function test_defaults_extraLibraries_mordred_isString(tc)
        % cfg.extraLibraries.mordred must be a non-empty string (pinned version).
            cfg = loadConfig();
            tc.verifyTrue(isfield(cfg.extraLibraries, "mordred"), ...
                "extraLibraries must have field 'mordred'");
            tc.verifyClass(cfg.extraLibraries.mordred, "string", ...
                "mordred version must be a MATLAB string");
            tc.verifyGreaterThan(strlength(cfg.extraLibraries.mordred), 0, ...
                "mordred version must be non-empty");
        end

        function test_defaults_extraLibraries_mordred_version_notDateBased(tc)
        % mordredcommunity uses 2.x.x versioning (NOT date-based like 2024.x.x).
        % Major version must be a single or double digit, not a 4-digit year.
        % This guards against the 2024.1.1 mistake (version that does not exist
        % on PyPI).
            cfg = loadConfig();
            ver = cfg.extraLibraries.mordred;
            parts = strsplit(ver, ".");
            major = str2double(parts{1});
            tc.verifyLessThan(major, 100, ...
                "mordred version must not use a 4-digit year prefix (e.g. 2024.1.1 " + ...
                "does not exist on PyPI). Use 2.x.x series instead.");
        end

        function test_defaults_extraLibraries_biopython_isString(tc)
        % cfg.extraLibraries.biopython must be a non-empty string (pinned version).
            cfg = loadConfig();
            tc.verifyTrue(isfield(cfg.extraLibraries, "biopython"), ...
                "extraLibraries must have field 'biopython'");
            tc.verifyClass(cfg.extraLibraries.biopython, "string", ...
                "biopython version must be a MATLAB string");
            tc.verifyGreaterThan(strlength(cfg.extraLibraries.biopython), 0, ...
                "biopython version must be non-empty");
        end

    end
end
