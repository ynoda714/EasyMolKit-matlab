classdef TestInstallTrack2 < matlab.unittest.TestCase
% TestInstallTrack2  Unit tests for src/+emk/+setup/installTrack2.m
%
% These tests cover logic paths that do NOT require network access,
% an actual pip install, or a base Python executable in PATH.
% Integration tests (real venv creation + pip install + import) belong
% in the smoke test suite.
%
% Coverage:
%   TC1 : Online environment -> notDesktop error
%   TC1b: notDesktop message guides user with pip install alternative
%   TC1c: notDesktop message for mdanalysis contains "MDAnalysis" (pip package name)
%   TC1d: notDesktop message for pymol contains "pymol-open-source" (pip package name)
%   TC2 : Unknown library name -> unknownLibrary error (before platform check)
%   TC2b: unknownLibrary message contains the bad library name
%   TC2c: unknownLibrary message lists supported names (mdanalysis, pymol)
%   TC2d: unknownLibrary message mentions openbabel recipe() guidance
%   TC2e: unknownLibrary fires BEFORE notDesktop even on MATLAB Online (ordering)
%   TC3 : "mdanalysis" is recognized (throws notDesktop on Online, not unknownLibrary)
%   TC4 : "pymol" is recognized (throws notDesktop on Online, not unknownLibrary)
%   TC5 : char input is accepted (arguments coercion) and reaches unknownLibrary
%   TC6 : Non-existent BasePython path -> basePythonNotFound error
%   TC6b: basePythonNotFound message contains the given path
%   TC7 : settings.example.json has python.external_path field
%   TC8 : Online guard fires for "pymol" (notDesktop before filesystem check)
%
% Notes:
%   Online mode is simulated via MATLAB_ONLINE env var, consistent with
%   TestInstallExtra and TestInstall.
%
% Run with:
%   addpath(genpath("src")); addpath(genpath("tests"));
%   results = run(TestInstallTrack2);

    methods (TestMethodSetup)
        function setupPath(tc) %#ok<MANU>
            addpath(genpath("src"));
        end
    end

    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TC1 / TC1b: Online guard
        % ------------------------------------------------------------------

        function test_onlineEnv_throwsNotDesktop(tc)
        % installTrack2() must be refused on MATLAB Online.
            setenv("MATLAB_ONLINE", "true");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            tc.verifyError(@() emk.setup.installTrack2("mdanalysis"), ...
                "emk:setup:installTrack2:notDesktop", ...
                "installTrack2() must throw notDesktop on MATLAB Online");
        end

        function test_onlineEnv_errorMessage_guidesPipInstall(tc)
        % Error message must mention pip so the user knows how to install on Online.
            setenv("MATLAB_ONLINE", "true");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            ME = tc.captureError(@() emk.setup.installTrack2("mdanalysis"));
            tc.assertNotEmpty(ME, "Expected a notDesktop error");
            tc.verifySubstring(lower(ME.message), "pip", ...
                "Error must mention pip install alternative for Online users");
        end

        function test_onlineEnv_notDesktopMessage_containsLibPipPackage_mdanalysis(tc)
        % Error message must contain the pip package name for mdanalysis
        % ("MDAnalysis") so users can copy-paste the install command.
            setenv("MATLAB_ONLINE", "true");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            ME = tc.captureError(@() emk.setup.installTrack2("mdanalysis"));
            tc.assertNotEmpty(ME, "Expected a notDesktop error");
            tc.verifySubstring(ME.message, "MDAnalysis", ...
                "notDesktop message must contain the pip package name for mdanalysis");
        end

        function test_onlineEnv_notDesktopMessage_containsLibPipPackage_pymol(tc)
        % Error message must contain the pip package name for pymol
        % ("pymol-open-source") so users can copy-paste the install command.
            setenv("MATLAB_ONLINE", "true");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            ME = tc.captureError(@() emk.setup.installTrack2("pymol"));
            tc.assertNotEmpty(ME, "Expected a notDesktop error");
            tc.verifySubstring(ME.message, "pymol-open-source", ...
                "notDesktop message must contain the pip package name for pymol");
        end

        % ------------------------------------------------------------------
        % TC2: Unknown library name (fires before platform check)
        % ------------------------------------------------------------------

        function test_unknownLibrary_throwsUnknownLibrary(tc)
        % An unrecognized name must throw unknownLibrary.  This must fire
        % before any Online/file-system check so the caller gets a precise error.
            tc.verifyError(@() emk.setup.installTrack2("nosuchlibrary"), ...
                "emk:setup:installTrack2:unknownLibrary", ...
                "Unrecognized library must throw unknownLibrary");
        end

        function test_unknownLibrary_errorMessage_containsName(tc)
        % Error message must echo the invalid name.
            badName = "nosuchlibrary";
            ME = tc.captureError(@() emk.setup.installTrack2(badName));
            tc.assertNotEmpty(ME, "Expected unknownLibrary error");
            tc.verifySubstring(ME.message, char(badName), ...
                "Error message must contain the rejected library name");
        end

        function test_unknownLibrary_errorMessage_listsSupportedNames(tc)
        % Error message must list the supported Track 2 names.
            ME = tc.captureError(@() emk.setup.installTrack2("bad_lib"));
            tc.assertNotEmpty(ME, "Expected unknownLibrary error");
            hasMDA = contains(ME.message, "mdanalysis");
            hasPMOL = contains(ME.message, "pymol");
            tc.verifyTrue(hasMDA || hasPMOL, ...
                "Error message must list at least one supported library name");
        end

        function test_unknownLibrary_errorMessage_mentionsOpenBabelRecipe(tc)
        % Error message must mention recipe() to guide openbabel users.
            ME = tc.captureError(@() emk.setup.installTrack2("openbabel"));
            tc.assertNotEmpty(ME, "Expected unknownLibrary error");
            tc.verifySubstring(lower(ME.message), "recipe", ...
                "Error must mention recipe() for openbabel guidance");
        end

        function test_unknownLibrary_firesBeforeNotDesktop_onOnline(tc)
        % CRITICAL ordering test: RDKit-free validations must fire first.
        % An unrecognized library name must throw unknownLibrary even when
        % MATLAB_ONLINE=true.  This proves unknownLibrary fires in STEP 1
        % (resolveLibInfo_) BEFORE the STEP 2 Online guard.
            setenv("MATLAB_ONLINE", "true");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            ME = tc.captureError(@() emk.setup.installTrack2("nosuchlibrary"));
            tc.assertNotEmpty(ME, "Expected unknownLibrary error");
            tc.verifyEqual(string(ME.identifier), ...
                "emk:setup:installTrack2:unknownLibrary", ...
                "unknownLibrary must fire BEFORE notDesktop even on MATLAB Online");
        end

        % ------------------------------------------------------------------
        % TC3 / TC4: Valid library names are recognized
        % ------------------------------------------------------------------

        function test_mdanalysis_isRecognized_notUnknownLibrary(tc)
        % "mdanalysis" is a valid name -- must reach notDesktop (Online env),
        % not unknownLibrary.
            setenv("MATLAB_ONLINE", "true");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            ME = tc.captureError(@() emk.setup.installTrack2("mdanalysis"));
            tc.assertNotEmpty(ME, "Expected an error");
            id = string(ME.identifier);
            tc.verifyNotEqual(id, ...
                "emk:setup:installTrack2:unknownLibrary", ...
                """mdanalysis"" must not throw unknownLibrary");
            tc.verifyEqual(id, ...
                "emk:setup:installTrack2:notDesktop", ...
                """mdanalysis"" must throw notDesktop on Online");
        end

        function test_pymol_isRecognized_notUnknownLibrary(tc)
        % "pymol" is a valid name -- must reach notDesktop (Online env).
            setenv("MATLAB_ONLINE", "true");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            ME = tc.captureError(@() emk.setup.installTrack2("pymol"));
            tc.assertNotEmpty(ME, "Expected an error");
            id = string(ME.identifier);
            tc.verifyNotEqual(id, ...
                "emk:setup:installTrack2:unknownLibrary", ...
                """pymol"" must not throw unknownLibrary");
            tc.verifyEqual(id, ...
                "emk:setup:installTrack2:notDesktop", ...
                """pymol"" must throw notDesktop on Online");
        end

        % ------------------------------------------------------------------
        % TC5: char input coercion
        % ------------------------------------------------------------------

        function test_charInput_isCoerced_reachesUnknownLibrary(tc)
        % char input must be accepted (coerced to string by arguments block)
        % and reach the unknownLibrary guard for an unrecognized name.
        % Uses single-quote char literal 'invalidcharname', NOT string "...".
            tc.verifyError(@() emk.setup.installTrack2('invalidcharname'), ...
                "emk:setup:installTrack2:unknownLibrary", ...
                "char input must be coerced to string and reach unknownLibrary guard");
        end

        % ------------------------------------------------------------------
        % TC6 / TC6b: Non-existent BasePython path
        % ------------------------------------------------------------------

        function test_nonExistentBasePython_throwsBasePythonNotFound(tc)
        % When BasePython points to a non-existent file, the function must
        % throw basePythonNotFound (Desktop-only, so guard against Online).
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            fakePath = fullfile(tempdir, "nonexistent_python_xyz.exe");
            tc.assumeFalse(isfile(fakePath), ...
                "Pre-condition: fake path must not exist");

            tc.verifyError( ...
                @() emk.setup.installTrack2("mdanalysis", BasePython=fakePath), ...
                "emk:setup:installTrack2:basePythonNotFound", ...
                "Non-existent BasePython must throw basePythonNotFound");
        end

        function test_nonExistentBasePython_errorMessage_containsPath(tc)
        % Error message must echo the given path so the user knows what was checked.
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            fakePath = fullfile(tempdir, "nonexistent_python_xyz.exe");
            tc.assumeFalse(isfile(fakePath), ...
                "Pre-condition: fake path must not exist");

            ME = tc.captureError( ...
                @() emk.setup.installTrack2("mdanalysis", BasePython=fakePath));
            tc.assertNotEmpty(ME, "Expected basePythonNotFound error");
            tc.verifySubstring(ME.message, char(fakePath), ...
                "Error message must contain the rejected path");
        end

        % ------------------------------------------------------------------
        % TC7: settings.example.json has external_path field
        % ------------------------------------------------------------------

        function test_settingsExample_hasPythonExternalPath(tc)
        % config/settings.example.json must declare the python.external_path
        % field so users know the key name (ADR-007).
            projectRoot = fileparts(fileparts(fileparts(fileparts( ...
                which("emk.setup.installTrack2")))));
            settingsPath = fullfile(projectRoot, "config", "settings.example.json");
            tc.assumeTrue(isfile(settingsPath), ...
                "Pre-condition: config/settings.example.json must exist");

            raw = jsondecode(fileread(settingsPath));
            tc.verifyTrue(isfield(raw, "python"), ...
                "settings.example.json must have a python section");
            tc.verifyTrue(isfield(raw.python, "external_path"), ...
                "settings.example.json python section must have external_path field");
        end

        % ------------------------------------------------------------------
        % TC8: Online guard fires for "pymol"
        % ------------------------------------------------------------------

        function test_pymolOnlineEnv_throwsNotDesktop(tc)
        % "pymol" must throw notDesktop on MATLAB Online (valid name check
        % must succeed before platform check fires).
            setenv("MATLAB_ONLINE", "true");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            tc.verifyError(@() emk.setup.installTrack2("pymol"), ...
                "emk:setup:installTrack2:notDesktop", ...
                "installTrack2(""pymol"") must throw notDesktop on Online");
        end

    end

    % ======================================================================
    % Private helpers
    % ======================================================================

    methods (Access = private)
        function ME = captureError(tc, fn) %#ok<INUSL>
        % Run fn and return the MException, or [] if no error was thrown.
            ME = [];
            try
                fn();
            catch e
                ME = e;
            end
        end
    end

end
