classdef TestVerify < matlab.unittest.TestCase
% TestVerify  Unit tests for src/+emk/+setup/verify.m
%
% verify() is a non-throwing diagnostic function.  Tests are grouped by
% what can be verified without assumptions about Python state (structural
% checks), and what requires Python to be loaded or unloaded (state checks).
%
% Coverage:
%   TC1 : Return value always has fields .python, .rdkit, .version
%   TC2 : Field types: .python -> logical, .rdkit -> logical, .version -> string
%   TC3 : verify() never throws in any Python state (non-throwing design)
%   TC4 : Python NOT loaded + python_env absent (env var override) -> .python=false, .rdkit=false  [assumeTrue NotLoaded]
%   TC4b: Same scenario -> .version is empty string                                                [assumeTrue NotLoaded]
%   TC5 : Python loaded -> .python=true                                                            [assumeTrue Loaded]
%   TC5b: Python loaded -> .version is non-empty string (verifyClass confirmed)                    [assumeTrue Loaded]
%
% Note (ADR-005): python_env absence is simulated via EMK_PYTHON_EMBEDDED_DIR env var override,
% NOT by cd() -- initPython() resolves python_env via mfilename('fullpath'), so cd() has no effect.
%
% Tests marked [assumeTrue] are skipped gracefully if the precondition is
% not met in the current MATLAB session.
%
% Run with:
%   addpath(genpath("src"));
%   results = run(TestVerify);

    methods (TestMethodSetup)
        function setupPath(tc) %#ok<MANU>
            addpath(genpath("src"));
        end
    end

    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TC1: Return struct field existence
        % ------------------------------------------------------------------

        function test_returnStruct_hasField_python(tc)
        % verify() must always return a struct with a .python field.
            result = emk.setup.verify();
            tc.verifyTrue(isfield(result, "python"), ...
                "Return struct must have field 'python'");
        end

        function test_returnStruct_hasField_rdkit(tc)
        % verify() must always return a struct with a .rdkit field.
            result = emk.setup.verify();
            tc.verifyTrue(isfield(result, "rdkit"), ...
                "Return struct must have field 'rdkit'");
        end

        function test_returnStruct_hasField_version(tc)
        % verify() must always return a struct with a .version field.
            result = emk.setup.verify();
            tc.verifyTrue(isfield(result, "version"), ...
                "Return struct must have field 'version'");
        end

        % ------------------------------------------------------------------
        % TC2: Field types
        % ------------------------------------------------------------------

        function test_pythonField_isLogical(tc)
        % result.python must be a scalar logical, not numeric or string.
            result = emk.setup.verify();
            tc.verifyClass(result.python, "logical", ...
                "result.python must be logical");
            tc.verifySize(result.python, [1, 1], ...
                "result.python must be scalar");
        end

        function test_rdkitField_isLogical(tc)
        % result.rdkit must be a scalar logical.
            result = emk.setup.verify();
            tc.verifyClass(result.rdkit, "logical", ...
                "result.rdkit must be logical");
            tc.verifySize(result.rdkit, [1, 1], ...
                "result.rdkit must be scalar");
        end

        function test_versionField_isString(tc)
        % result.version must be a scalar string (MATLAB string type).
            result = emk.setup.verify();
            tc.verifyClass(result.version, "string", ...
                "result.version must be a MATLAB string");
            tc.verifySize(result.version, [1, 1], ...
                "result.version must be scalar");
        end

        % ------------------------------------------------------------------
        % TC3: Non-throwing design - the most critical contract of verify()
        % ------------------------------------------------------------------

        function test_neverThrows(tc)
        % verify() must never throw regardless of Python environment state.
        % This is the fundamental non-throwing contract (non-throwing design).
        % No Python assumptions required: passes whether Python is loaded or not.
            ME = [];
            try
                emk.setup.verify();
            catch e
                ME = e;
            end
            tc.verifyEmpty(ME, ...
                "verify() must never throw in any Python state");
        end

        % ------------------------------------------------------------------
        % TC4: Python NOT loaded, no python_env/ -> all-false result
        % [assumeTrue: Python must be NotLoaded]
        % ------------------------------------------------------------------

        function test_initPythonFails_pythonField_isFalse(tc)
        % When Python is not loaded and initPython() fails (python_env absent),
        % verify() must return .python = false.
        % python_env absence is simulated via env var override (ADR-005):
        % initPython() uses mfilename('fullpath') for resolution, so cd() has no effect.
            tc.assumeTrue(strcmp(string(pyenv().Status), "NotLoaded"), ...
                "Skipped: Python already loaded; initPython() will not be invoked");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            setenv("EMK_PYTHON_EMBEDDED_DIR", "nonexistent_python_env_for_test");
            tc.addTeardown(@() setenv("EMK_PYTHON_EMBEDDED_DIR", ""));

            result = emk.setup.verify();

            tc.verifyFalse(result.python, ...
                "result.python must be false when initPython() fails");
        end

        function test_initPythonFails_rdkitField_isFalse(tc)
        % When Python is unavailable, result.rdkit must also be false.
        % python_env absence is simulated via env var override (ADR-005).
            tc.assumeTrue(strcmp(string(pyenv().Status), "NotLoaded"), ...
                "Skipped: Python already loaded; initPython() will not be invoked");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            setenv("EMK_PYTHON_EMBEDDED_DIR", "nonexistent_python_env_for_test");
            tc.addTeardown(@() setenv("EMK_PYTHON_EMBEDDED_DIR", ""));

            result = emk.setup.verify();

            tc.verifyFalse(result.rdkit, ...
                "result.rdkit must be false when Python is unavailable");
        end

        function test_initPythonFails_versionField_isEmpty(tc)
        % When Python is unavailable, result.version must be "" (not erroring out).
        % python_env absence is simulated via env var override (ADR-005).
            tc.assumeTrue(strcmp(string(pyenv().Status), "NotLoaded"), ...
                "Skipped: Python already loaded; initPython() will not be invoked");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            setenv("EMK_PYTHON_EMBEDDED_DIR", "nonexistent_python_env_for_test");
            tc.addTeardown(@() setenv("EMK_PYTHON_EMBEDDED_DIR", ""));

            result = emk.setup.verify();

            tc.verifyEqual(result.version, "", ...
                "result.version must be empty string when Python is unavailable");
            tc.verifyClass(result.version, "string", ...
                "result.version must be a MATLAB string even when Python is unavailable");
        end

        % ------------------------------------------------------------------
        % TC5: Python loaded -> positive checks
        % [assumeTrue: Python must be Loaded]
        % ------------------------------------------------------------------

        function test_pythonLoaded_pythonField_isTrue(tc)
        % When Python is already loaded, .python must be true.
            tc.assumeTrue(~strcmp(string(pyenv().Status), "NotLoaded"), ...
                "Skipped: Python must be loaded for this test");

            result = emk.setup.verify();

            tc.verifyTrue(result.python, ...
                "result.python must be true when Python is loaded");
        end

        function test_pythonLoaded_versionField_nonEmpty(tc)
        % When Python is loaded, .version must report a non-empty string.
        % Also verify type (string, scalar) since this is a state-dependent path.
            tc.assumeTrue(~strcmp(string(pyenv().Status), "NotLoaded"), ...
                "Skipped: Python must be loaded for this test");

            result = emk.setup.verify();

            tc.verifyClass(result.version, "string", ...
                "result.version must be a MATLAB string when Python is loaded");
            tc.verifyGreaterThan(strlength(result.version), 0, ...
                "result.version must be non-empty when Python is loaded");
        end

    end
end
