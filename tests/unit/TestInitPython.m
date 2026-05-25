classdef TestInitPython < matlab.unittest.TestCase
% TestInitPython  Unit tests for src/+emk/+setup/initPython.m
%
% Run with:
%   addpath(genpath("src"));
%   results = run(TestInitPython);
%
% Coverage:
%   - Desktop + python_env absent => notInstalled error (TC1-TC3)
%   - Double-call guard fires BEFORE file check (TC4)
%   - Python already loaded => silent return, no re-configuration (TC5-TC6)
%   - Online mode => file check bypassed, notInstalled never thrown (TC7)
%   - Track 2: external_path set (non-empty) => useExternal() called,
%              not notInstalled (EXT1)
%   - Track 2: external_path empty => falls through to embedded python
%              (Track 1), notInstalled when python_env absent (EXT2)
%   - Track 2: external_path whitespace-only => treated as empty => Track 1,
%              notInstalled when python_env absent (EXT3)
%
% Tests that require Python to already be loaded use assumeTrue to skip
% gracefully in a fresh MATLAB session where Python has not been configured.
%
% Note: initPython() resolves python_env via mfilename('fullpath')
% rather than pwd.  TC1-TC4, TC7, and EXT1-EXT3 simulate an absent python_env
% by setting EMK_PYTHON_EMBEDDED_DIR to a nonexistent path instead of cd().

    methods (TestMethodSetup)
        function setupPath(tc) %#ok<MANU>
            addpath(genpath("src"));
        end
    end

    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TC1-TC3: Desktop + NotLoaded + python_env absent
        % ------------------------------------------------------------------

        function test_desktop_noPythonEnv_throwsNotInstalled(tc)
        % Desktop mode with no python_env/python.exe must throw notInstalled.
        % Python must NOT be loaded for this path to be reached (guard check).
        % python_env absence is simulated via env var override (EMK_PYTHON_EMBEDDED_DIR).
            tc.assumeTrue(strcmp(string(pyenv().Status), "NotLoaded"), ...
                "Skipped: Python is already loaded; guard would fire before file check");

            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            setenv("EMK_PYTHON_EMBEDDED_DIR", "nonexistent_python_env_for_test");
            tc.addTeardown(@() setenv("EMK_PYTHON_EMBEDDED_DIR", ""));

            tc.verifyError(@() emk.setup.initPython(), ...
                "emk:setup:initPython:notInstalled", ...
                "Must throw notInstalled when python_env/python.exe is absent on Desktop");
        end

        % ------------------------------------------------------------------
        function test_desktop_noPythonEnv_errorMessage_containsPath(tc)
        % Error message must include the expected path to aid debugging.
        % python_env absence is simulated via env var override (EMK_PYTHON_EMBEDDED_DIR).
            tc.assumeTrue(strcmp(string(pyenv().Status), "NotLoaded"), ...
                "Skipped: Python is already loaded; guard would fire before file check");

            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            setenv("EMK_PYTHON_EMBEDDED_DIR", "nonexistent_python_env_for_test");
            tc.addTeardown(@() setenv("EMK_PYTHON_EMBEDDED_DIR", ""));

            ME = tc.captureError(@() emk.setup.initPython());
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "nonexistent_python_env_for_test", ...
                "Error message must contain the python_env path for debugging");
        end

        % ------------------------------------------------------------------
        function test_desktop_noPythonEnv_errorMessage_mentionsInstall(tc)
        % Error message must include 'install' to guide the user to the fix.
        % python_env absence is simulated via env var override (EMK_PYTHON_EMBEDDED_DIR).
            tc.assumeTrue(strcmp(string(pyenv().Status), "NotLoaded"), ...
                "Skipped: Python is already loaded; guard would fire before file check");

            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            setenv("EMK_PYTHON_EMBEDDED_DIR", "nonexistent_python_env_for_test");
            tc.addTeardown(@() setenv("EMK_PYTHON_EMBEDDED_DIR", ""));

            ME = tc.captureError(@() emk.setup.initPython());
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(lower(ME.message), "install", ...
                "Error message must mention install() to guide the user");
        end

        % ------------------------------------------------------------------
        % TC4: Guard fires BEFORE file check (critical idempotency property)
        % ------------------------------------------------------------------

        function test_guardFiresBeforeFileCheck_noNotInstalledWhenLoaded(tc)
        % When Python is already loaded, the double-call guard must fire
        % BEFORE the file-existence check.  Even on Desktop with no
        % python_env/python.exe present, no notInstalled error must occur.
        %
        % This test directly verifies the "guard first" implementation order.
        % python_env absence is simulated via env var override (EMK_PYTHON_EMBEDDED_DIR).
            tc.assumeTrue(~strcmp(string(pyenv().Status), "NotLoaded"), ...
                "Skipped: Python must be loaded to verify guard-before-file-check");

            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            setenv("EMK_PYTHON_EMBEDDED_DIR", "nonexistent_python_env_for_test");
            tc.addTeardown(@() setenv("EMK_PYTHON_EMBEDDED_DIR", ""));

            ME = tc.captureError(@() emk.setup.initPython());
            tc.verifyEmpty(ME, ...
                "Guard must fire before file check; no error when Python already loaded");
        end

        % ------------------------------------------------------------------
        % TC5-TC6: Python already loaded => silent return
        % ------------------------------------------------------------------

        function test_alreadyLoaded_desktop_silentReturn(tc)
        % When Python is already loaded (Desktop), initPython must return
        % without error and without reconfiguring pyenv.
            tc.assumeTrue(~strcmp(string(pyenv().Status), "NotLoaded"), ...
                "Skipped: Python must be loaded to test double-call guard");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            statusBefore  = string(pyenv().Status);
            versionBefore = string(pyenv().Version);

            ME = tc.captureError(@() emk.setup.initPython());
            tc.verifyEmpty(ME, ...
                "Must not throw when Python is already loaded (Desktop)");

            % pyenv status and version must not change (no re-initialization)
            tc.verifyEqual(string(pyenv().Status), statusBefore, ...
                "pyenv Status must be unchanged after guard fires");
            tc.verifyClass(string(pyenv().Version), "string", ...
                "pyenv Version must be a MATLAB string after guard fires");
            tc.verifyGreaterThan(strlength(string(pyenv().Version)), 0, ...
                "pyenv Version must be non-empty after guard fires");
            tc.verifyEqual(string(pyenv().Version), versionBefore, ...
                "pyenv Version must be unchanged after guard fires");
        end

        % ------------------------------------------------------------------
        function test_alreadyLoaded_online_silentReturn(tc)
        % When Python is already loaded (Online mode), initPython must return
        % without error.  Verifies that Online mode + guard combination works.
            tc.assumeTrue(~strcmp(string(pyenv().Status), "NotLoaded"), ...
                "Skipped: Python must be loaded to test double-call guard (Online)");
            setenv("MATLAB_ONLINE", "true");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            statusBefore = string(pyenv().Status);

            ME = tc.captureError(@() emk.setup.initPython());
            tc.verifyEmpty(ME, ...
                "Must not throw when Python is already loaded (Online mode)");

            % Guard must not alter pyenv state
            tc.verifyEqual(string(pyenv().Status), statusBefore, ...
                "pyenv Status must be unchanged after guard fires (Online mode)");
        end

        % ------------------------------------------------------------------
        % TC7: Online mode bypasses file check
        % ------------------------------------------------------------------

        function test_online_noPythonEnv_notInstalledNeverThrown(tc)
        % In Online mode, the file-existence check is skipped entirely.
        % Even without python_env/python.exe, notInstalled must NOT be thrown.
        %
        % Verification strategy: capture any error that occurs.  If an error
        % is thrown, its identifier must NOT be notInstalled.  pyenvFailed is
        % acceptable (system Python absent); success is also acceptable.
        % python_env absence is simulated via env var override (EMK_PYTHON_EMBEDDED_DIR).
            setenv("MATLAB_ONLINE", "true");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            setenv("EMK_PYTHON_EMBEDDED_DIR", "nonexistent_python_env_for_test");
            tc.addTeardown(@() setenv("EMK_PYTHON_EMBEDDED_DIR", ""));

            ME = tc.captureError(@() emk.setup.initPython());

            if ~isempty(ME)
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:setup:initPython:notInstalled", ...
                    "Online mode must never throw notInstalled regardless of python_env");
            end
            % No error: trivially passes (Online mode succeeded or guard fired)
        end

        % ------------------------------------------------------------------
        % EXT1-EXT3: Track 2 external_path detection (ADR-007)
        % ------------------------------------------------------------------

        function test_externalPathSet_nonEmpty_desktop_callsUseExternal(tc)
        % EXT1: When external_path is set to a non-empty string (via env var)
        % on Desktop with Python NotLoaded, useExternal() must be called
        % instead of the embedded Python path check.
        %
        % The path does not exist, so useExternal throws fileNotFound.
        % This verifies that we reach useExternal (Track 2) and NOT the
        % embedded python_env check (Track 1, which would throw notInstalled).
            tc.assumeTrue(strcmp(string(pyenv().Status), "NotLoaded"), ...
                "Skipped: Python already loaded; guard fires before external_path check");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            setenv("EMK_PYTHON_EXTERNAL_PATH", "C:\nonexistent_t2_python_for_test\python.exe");
            tc.addTeardown(@() setenv("EMK_PYTHON_EXTERNAL_PATH", ""));
            setenv("EMK_PYTHON_EMBEDDED_DIR", "nonexistent_python_env_for_test");
            tc.addTeardown(@() setenv("EMK_PYTHON_EMBEDDED_DIR", ""));

            ME = tc.captureError(@() emk.setup.initPython());

            tc.assertNotEmpty(ME, "Expected an error since the path does not exist");
            % Must NOT be notInstalled (that belongs to Track 1)
            tc.verifyNotEqual(string(ME.identifier), ...
                "emk:setup:initPython:notInstalled", ...
                "Non-empty external_path must route to useExternal, not Track 1 notInstalled");
            % Must be a useExternal error (Track 2 path was taken)
            tc.verifySubstring(ME.identifier, "useExternal", ...
                "Non-empty external_path: error must originate from useExternal()");
        end

        function test_externalPathEmpty_desktop_fallsBackToEmbedded(tc)
        % EXT2: When external_path is an empty string (env var = ""),
        % initPython must ignore Track 2 and proceed with Track 1 (embedded
        % python_env check).  With python_env absent, notInstalled is thrown.
            tc.assumeTrue(strcmp(string(pyenv().Status), "NotLoaded"), ...
                "Skipped: Python already loaded; guard fires before external_path check");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            setenv("EMK_PYTHON_EXTERNAL_PATH", "");
            tc.addTeardown(@() setenv("EMK_PYTHON_EXTERNAL_PATH", ""));
            setenv("EMK_PYTHON_EMBEDDED_DIR", "nonexistent_python_env_for_test");
            tc.addTeardown(@() setenv("EMK_PYTHON_EMBEDDED_DIR", ""));

            ME = tc.captureError(@() emk.setup.initPython());

            tc.assertNotEmpty(ME, "Expected notInstalled with absent embedded Python");
            tc.verifyEqual(string(ME.identifier), ...
                "emk:setup:initPython:notInstalled", ...
                "Empty external_path must fall back to Track 1 and throw notInstalled");
        end

        function test_externalPathWhitespace_desktop_treatedAsEmpty(tc)
        % EXT3: Whitespace-only external_path must be treated as empty
        % (strtrim check) and fall back to Track 1.
        % With python_env absent, notInstalled must be thrown.
            tc.assumeTrue(strcmp(string(pyenv().Status), "NotLoaded"), ...
                "Skipped: Python already loaded; guard fires before external_path check");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            setenv("EMK_PYTHON_EXTERNAL_PATH", "   ");
            tc.addTeardown(@() setenv("EMK_PYTHON_EXTERNAL_PATH", ""));
            setenv("EMK_PYTHON_EMBEDDED_DIR", "nonexistent_python_env_for_test");
            tc.addTeardown(@() setenv("EMK_PYTHON_EMBEDDED_DIR", ""));

            ME = tc.captureError(@() emk.setup.initPython());

            tc.assertNotEmpty(ME, "Expected notInstalled: whitespace-only treated as empty");
            tc.verifyEqual(string(ME.identifier), ...
                "emk:setup:initPython:notInstalled", ...
                "Whitespace-only external_path must be treated as empty (Track 1)");
        end

    end

    % ======================================================================
    methods (Access = private)

        function ME = captureError(~, fcn)
        % Call fcn and return the MException if thrown, or [] if no error.
            ME = [];
            try
                fcn();
            catch e
                ME = e;
            end
        end

    end
end
