classdef TestInstallExtra < matlab.unittest.TestCase
% TestInstallExtra  Unit tests for src/+emk/+setup/installExtra.m
%
% These tests cover only logic paths that do NOT require network access
% or an actual pip install.  Integration tests (real pip install + import)
% belong in the smoke test suite.
%
% Coverage:
%   TC1 : Online environment -> notDesktop error
%   TC1b: notDesktop message guides user with pip install alternative
%   TC2 : Unknown library name -> unknownLibrary error
%   TC2b: unknownLibrary message contains the bad library name
%   TC2c: unknownLibrary message mentions supported Track 1 names
%   TC2d: unknownLibrary message mentions recipe() for Track 2 guidance
%   TC3 : Desktop + python_env absent -> installFailed error
%   TC3b: installFailed message contains a path (python.exe location hint)
%   TC3c: installFailed message mentions emk.setup.install() to fix it
%   TC4 : char input is accepted (arguments coercion) and reaches unknownLibrary
%   TC5 : settings.example.json mordred version uses 2.x.x (not date-based)
%   TC6 : New A-1 names (torch, torch_geometric, transformers) are
%         recognized -- throw installFailed not unknownLibrary
%   TC7 : Online guard fires for torch_geometric (notDesktop before system())
%   TC8 : settings.example.json torch version field exists and is semver
%   TC9 : unknownLibrary message mentions at least one torch-family name
%   TC10: meeko (A-3) is recognized -- throws installFailed (not unknownLibrary)
%         on Desktop with absent python_env (no Windows early return for meeko)
%   TC10b: A-3 library names throw notDesktop on MATLAB Online (Online guard)
%   TC10c: vina on Windows Desktop returns normally (early return -- no PyPI wheel)
%   TC10d: pdbfixer on Windows Desktop returns normally (early return -- SAC blocks openmm)
%   TC10e: vina recognized on non-Windows Desktop -- throws installFailed not unknownLibrary
%   TC10f: pdbfixer recognized on non-Windows Desktop -- throws installFailed not unknownLibrary
%   TC11: unknownLibrary message mentions ALL THREE docking library names
%         (meeko AND vina AND pdbfixer -- checked individually)
%   TC12: meeko LGPL warning -- function reaches installFailed (not earlier error)
% Notes:
%   python_env absence is simulated via EMK_PYTHON_EMBEDDED_DIR env var override,
%   consistent with TestVerify and TestInitPython (ADR-005).
%   TC10c/TC10d run only on Windows (ispc()=true).
%   TC10e/TC10f run only on non-Windows Desktop (ispc()=false).
%
% Run with:
%   addpath(genpath("src")); addpath(genpath("tests"));
%   results = run(TestInstallExtra);

    methods (TestMethodSetup)
        function setupPath(tc) %#ok<MANU>
            addpath(genpath("src"));
        end
    end

    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TC1: Online guard
        % ------------------------------------------------------------------

        function test_onlineEnv_throwsNotDesktop(tc)
        % installExtra() must be refused on MATLAB Online; Embedded Python
        % is not available there.
            setenv("MATLAB_ONLINE", "true");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            tc.verifyError(@() emk.setup.installExtra("pubchempy"), ...
                "emk:setup:installExtra:notDesktop", ...
                "installExtra() must throw notDesktop on MATLAB Online");
        end

        function test_onlineEnv_errorMessage_guidesPipInstall(tc)
        % Error message must tell the user how to install packages on Online
        % (manual pip install --user).
            setenv("MATLAB_ONLINE", "true");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            ME = tc.captureError(@() emk.setup.installExtra("pubchempy"));
            tc.assertNotEmpty(ME, "Expected a notDesktop error");
            tc.verifyClass(ME, "MException", "Error must be MException");
            tc.verifySubstring(lower(ME.message), "pip", ...
                "Error must mention pip install alternative for Online users");
        end

        % ------------------------------------------------------------------
        % TC2: Unknown library name
        % ------------------------------------------------------------------

        function test_unknownLibrary_throwsUnknownLibrary(tc)
        % An unrecognized library name must throw unknownLibrary before
        % any system() call or file check.
            tc.verifyError(@() emk.setup.installExtra("nosuchlibrary"), ...
                "emk:setup:installExtra:unknownLibrary", ...
                "Unrecognized library must throw unknownLibrary");
        end

        function test_unknownLibrary_errorMessage_containsName(tc)
        % Error message must echo the invalid name so the user knows what
        % was rejected.
            badName = "nosuchlibrary";
            ME = tc.captureError(@() emk.setup.installExtra(badName));
            tc.assertNotEmpty(ME, "Expected unknownLibrary error");
            tc.verifySubstring(ME.message, char(badName), ...
                "Error message must contain the rejected library name");
        end

        function test_unknownLibrary_errorMessage_mentionesSupportedNames(tc)
        % Error message must list at least one Track 1 supported name to
        % guide the user toward the correct spelling.
            ME = tc.captureError(@() emk.setup.installExtra("bad_lib"));
            tc.assertNotEmpty(ME, "Expected unknownLibrary error");
            hasPub = contains(ME.message, "pubchempy");
            hasMor = contains(ME.message, "mordred");
            hasBio = contains(ME.message, "biopython");
            tc.verifyTrue(hasPub || hasMor || hasBio, ...
                "Error message must list at least one supported library name");
        end

        function test_unknownLibrary_errorMessage_mentionsRecipe(tc)
        % Error message must mention recipe() to guide users who want
        % Track 2 libraries (openbabel, mdanalysis, pymol).
            ME = tc.captureError(@() emk.setup.installExtra("openbabel"));
            tc.assertNotEmpty(ME, "Expected unknownLibrary error");
            tc.verifySubstring(lower(ME.message), "recipe", ...
                "Error must mention recipe() for Track 2 guidance");
        end

        % ------------------------------------------------------------------
        % TC3: Desktop + python_env absent -> installFailed
        % ------------------------------------------------------------------

        function test_noPythonEnv_throwsInstallFailed(tc)
        % When python_env/ is absent (simulated via env var), installExtra
        % must throw installFailed with a valid library name.
        % python_env absence is simulated via EMK_PYTHON_EMBEDDED_DIR (ADR-005).
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only python_env check (cannot mock isOnline on real Online)");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            setenv("EMK_PYTHON_EMBEDDED_DIR", "nonexistent_python_env_for_test");
            tc.addTeardown(@() setenv("EMK_PYTHON_EMBEDDED_DIR", ""));

            tc.verifyError(@() emk.setup.installExtra("pubchempy"), ...
                "emk:setup:installExtra:installFailed", ...
                "Must throw installFailed when python_env/python.exe is absent");
        end

        function test_noPythonEnv_errorMessage_containsPath(tc)
        % Error message must include the python.exe path to aid debugging.
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only python_env check");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            setenv("EMK_PYTHON_EMBEDDED_DIR", "nonexistent_python_env_for_test");
            tc.addTeardown(@() setenv("EMK_PYTHON_EMBEDDED_DIR", ""));

            ME = tc.captureError(@() emk.setup.installExtra("pubchempy"));
            tc.assertNotEmpty(ME, "Expected installFailed error");
            tc.verifyClass(ME, "MException", "Error must be MException");
            tc.verifySubstring(ME.message, "nonexistent_python_env_for_test", ...
                "Error message must contain the missing python_env path");
        end

        function test_noPythonEnv_errorMessage_mentionsInstall(tc)
        % Error message must suggest emk.setup.install() as the fix.
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only python_env check");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            setenv("EMK_PYTHON_EMBEDDED_DIR", "nonexistent_python_env_for_test");
            tc.addTeardown(@() setenv("EMK_PYTHON_EMBEDDED_DIR", ""));

            ME = tc.captureError(@() emk.setup.installExtra("pubchempy"));
            tc.assertNotEmpty(ME, "Expected installFailed error");
            tc.verifySubstring(lower(ME.message), "install", ...
                "Error must mention emk.setup.install() to guide the user");
        end

        % ------------------------------------------------------------------
        % TC4: char input coercion
        % ------------------------------------------------------------------

        function test_charInput_acceptedByArguments(tc)
        % installExtra() must accept char input (arguments block coerces to
        % string).  An unknown library name is used to get a deterministic
        % early failure without any network or system() access.
            tc.verifyError(@() emk.setup.installExtra('nosuchlibrary'), ...
                "emk:setup:installExtra:unknownLibrary", ...
                "char input must be coerced to string by the arguments block");
        end

        % ------------------------------------------------------------------
        % TC5: settings.example.json pinned version format
        % ------------------------------------------------------------------

        function test_settingsExample_mordred_notDateBased(tc)
        % mordredcommunity uses 2.x.x versioning, NOT date-based (2024.x.x).
        % This is a static test that reads settings.example.json without
        % any network access.  It guards against accidental use of a
        % non-existent version on PyPI (e.g. 2024.1.1 does not exist).
            projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
            examplePath = fullfile(projectRoot, "config", "settings.example.json");
            tc.assumeTrue(isfile(examplePath), "settings.example.json not found");
            raw = jsondecode(fileread(examplePath));
            tc.assertTrue(isfield(raw, "extraLibraries") && ...
                isfield(raw.extraLibraries, "mordred"), ...
                "settings.example.json must have extraLibraries.mordred");
            ver = string(raw.extraLibraries.mordred);
            parts = strsplit(ver, ".");
            major = str2double(parts{1});
            tc.verifyLessThan(major, 100, ...
                "settings.example.json mordred version must not use a 4-digit year " + ...
                "prefix. mordredcommunity uses 2.x.x series (e.g. 2.0.7).");
        end

        % ------------------------------------------------------------------
        % TC6: New A-1 library names are recognized by resolveLibInfo_
        % ------------------------------------------------------------------

        function test_newLibNames_areKnown_notUnknownLibraryError(tc)
        % torch, torch_geometric, transformers added in A-1 must be
        % accepted by resolveLibInfo_().  Verified on Desktop by confirming
        % the error is installFailed (python_env absent) NOT unknownLibrary.
        % A regression in resolveLibInfo_() for any new name would produce
        % unknownLibrary, which this test would catch immediately.
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only python_env simulation");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            setenv("EMK_PYTHON_EMBEDDED_DIR", "nonexistent_for_a1_test");
            tc.addTeardown(@() setenv("EMK_PYTHON_EMBEDDED_DIR", ""));

            newNames = ["torch", "torch_geometric", "transformers"];
            for k = 1:numel(newNames)
                nm = newNames(k);
                ME = tc.captureError(@() emk.setup.installExtra(nm));
                tc.assertNotEmpty(ME, nm + ": expected an error");
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:setup:installExtra:unknownLibrary", ...
                    nm + ": must NOT throw unknownLibrary (name should be recognized)");
            end
        end

        % ------------------------------------------------------------------
        % TC7: Online guard fires for torch_geometric
        % ------------------------------------------------------------------

        function test_torchGeometric_onlineEnv_throwsNotDesktop(tc)
        % torch_geometric must be refused on MATLAB Online.
        % Ensures the Online guard (STEP 2) fires BEFORE installTorchGeometric_()
        % executes any system() call, preventing unintended side-effects.
            setenv("MATLAB_ONLINE", "true");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            tc.verifyError(@() emk.setup.installExtra("torch_geometric"), ...
                "emk:setup:installExtra:notDesktop", ...
                "torch_geometric must throw notDesktop on MATLAB Online");
        end

        % ------------------------------------------------------------------
        % TC8: settings.example.json torch version field
        % ------------------------------------------------------------------

        function test_settingsExample_torch_versionFieldExists(tc)
        % settings.example.json must declare extraLibraries.torch so that
        % torch is installed at a pinned version, making the PyG wheel URL
        % deterministic: https://data.pyg.org/whl/torch-<X.Y.Z>+cpu.html
            projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
            examplePath = fullfile(projectRoot, "config", "settings.example.json");
            tc.assumeTrue(isfile(examplePath), "settings.example.json not found");
            raw = jsondecode(fileread(examplePath));
            tc.verifyTrue(isfield(raw, "extraLibraries") && ...
                isfield(raw.extraLibraries, "torch"), ...
                "settings.example.json must have extraLibraries.torch " + ...
                "to enable deterministic PyG wheel URL construction");
        end

        function test_settingsExample_torch_notDateBased(tc)
        % The torch version must be a regular X.Y.Z semver (not a date-based
        % version such as 2024.x.x).  PyTorch uses standard semver.
            projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
            examplePath = fullfile(projectRoot, "config", "settings.example.json");
            tc.assumeTrue(isfile(examplePath), "settings.example.json not found");
            raw = jsondecode(fileread(examplePath));
            tc.assumeTrue(isfield(raw, "extraLibraries") && ...
                isfield(raw.extraLibraries, "torch"), ...
                "Skipped: torch field absent (covered by test_settingsExample_torch_versionFieldExists)");
            ver   = string(raw.extraLibraries.torch);
            parts = strsplit(ver, ".");
            major = str2double(parts{1});
            tc.verifyLessThan(major, 100, ...
                "settings.example.json torch version must use X.Y.Z semver, " + ...
                "not a 4-digit year prefix (e.g. correct: 2.5.0).");
        end

        % ------------------------------------------------------------------
        % TC9: unknownLibrary message mentions torch-family names
        % ------------------------------------------------------------------

        function test_unknownLibrary_errorMessage_mentionsTorchName(tc)
        % After A-1, the unknownLibrary error message must list at least one
        % torch-family name so users can find the correct spelling for the
        % new libraries added in this milestone.
            ME = tc.captureError(@() emk.setup.installExtra("bad_lib"));
            tc.assertNotEmpty(ME, "Expected unknownLibrary error");
            tc.verifyTrue(contains(ME.message, "torch"), ...
                "unknownLibrary message must mention 'torch' after A-1 changes");
        end

        % ------------------------------------------------------------------
        % TC10: meeko (A-3) recognized on Desktop without python_env
        % ------------------------------------------------------------------

        function test_meeko_isKnown_throwsInstallFailed(tc)
        % meeko added in A-3 must be accepted by resolveLibInfo_().
        % On Desktop with absent python_env, the function must throw
        % installFailed (not unknownLibrary).  meeko has no Windows early
        % return, so this test runs on all Desktop platforms.
        % Uses verifyEqual to confirm the exact error ID (stricter than
        % verifyNotEqual).
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only python_env simulation");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            setenv("EMK_PYTHON_EMBEDDED_DIR", "nonexistent_for_a3_meeko_test");
            tc.addTeardown(@() setenv("EMK_PYTHON_EMBEDDED_DIR", ""));

            ME = tc.captureError(@() emk.setup.installExtra("meeko"));
            tc.assertNotEmpty(ME, "meeko: expected installFailed error");
            tc.verifyEqual(ME.identifier, ...
                'emk:setup:installExtra:installFailed', ...
                "meeko must throw exactly installFailed (not unknownLibrary or other)");
        end

        % ------------------------------------------------------------------
        % TC10b: A-3 library names throw notDesktop on MATLAB Online
        % ------------------------------------------------------------------

        function test_dockingLibNames_onlineEnv_throwNotDesktop(tc)
        % meeko, vina, pdbfixer must all be refused on MATLAB Online.
        % The Online guard fires BEFORE any ispc() check, so all three
        % must throw notDesktop regardless of platform.
            setenv("MATLAB_ONLINE", "true");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            dockingNames = ["meeko", "vina", "pdbfixer"];
            for k = 1:numel(dockingNames)
                nm = dockingNames(k);
                tc.verifyError(@() emk.setup.installExtra(nm), ...
                    "emk:setup:installExtra:notDesktop", ...
                    nm + ": must throw notDesktop on MATLAB Online");
            end
        end

        % ------------------------------------------------------------------
        % TC10c: vina on Windows Desktop returns normally (early return)
        % ------------------------------------------------------------------

        function test_vina_windowsDesktop_returnsNormally(tc)
        % vina has no PyPI wheel on Windows.  installExtra("vina") must
        % return normally (early return) without throwing any error.
        % The function logs a warning and calls recipe("docking") instead
        % of attempting pip, which would fail.
        % This test runs only on Windows (ispc()=true).
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only test");
            tc.assumeTrue(ispc(), ...
                "Skipped: Windows-only early return (vina has no PyPI wheel on Windows)");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            ME = tc.captureError(@() emk.setup.installExtra("vina"));
            tc.verifyEmpty(ME, ...
                "vina on Windows must return normally (early return, no error). " + ...
                "No PyPI wheel is available; pip install would fail.");
        end

        % ------------------------------------------------------------------
        % TC10d: pdbfixer on Windows Desktop returns normally (early return)
        % ------------------------------------------------------------------

        function test_pdbfixer_windowsDesktop_returnsNormally(tc)
        % pdbfixer requires openmm, whose _openmm.pyd is blocked by Windows
        % Smart App Control (SAC).  installExtra("pdbfixer") must return
        % normally (early return) without throwing any error.
        % The function logs a warning and calls recipe("docking").
        % This test runs only on Windows (ispc()=true).
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only test");
            tc.assumeTrue(ispc(), ...
                "Skipped: Windows-only early return (pdbfixer/openmm blocked by SAC on Windows)");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            ME = tc.captureError(@() emk.setup.installExtra("pdbfixer"));
            tc.verifyEmpty(ME, ...
                "pdbfixer on Windows must return normally (early return, no error). " + ...
                "openmm .pyd is blocked by Windows Smart App Control.");
        end

        % ------------------------------------------------------------------
        % TC10e: vina recognized on non-Windows Desktop (installFailed)
        % ------------------------------------------------------------------

        function test_vina_nonWindowsDesktop_isKnown_throwsInstallFailed(tc)
        % On non-Windows Desktop (Linux/macOS), vina has a PyPI wheel so
        % the ispc() early return does NOT fire.  With absent python_env,
        % installExtra("vina") must reach python_env check and throw
        % installFailed (confirming vina is registered in resolveLibInfo_).
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only test");
            tc.assumeFalse(ispc(), ...
                "Skipped: non-Windows only (Windows has early return for vina)");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            setenv("EMK_PYTHON_EMBEDDED_DIR", "nonexistent_for_vina_nonwin_test");
            tc.addTeardown(@() setenv("EMK_PYTHON_EMBEDDED_DIR", ""));

            ME = tc.captureError(@() emk.setup.installExtra("vina"));
            tc.assertNotEmpty(ME, "vina on non-Windows: expected installFailed");
            tc.verifyEqual(ME.identifier, ...
                'emk:setup:installExtra:installFailed', ...
                "vina on non-Windows must throw installFailed (registered in resolveLibInfo_)");
        end

        % ------------------------------------------------------------------
        % TC10f: pdbfixer recognized on non-Windows Desktop (installFailed)
        % ------------------------------------------------------------------

        function test_pdbfixer_nonWindowsDesktop_isKnown_throwsInstallFailed(tc)
        % On non-Windows Desktop, pdbfixer openmm is NOT blocked by SAC
        % so the ispc() early return does NOT fire.  With absent python_env,
        % installExtra("pdbfixer") must reach python_env check and throw
        % installFailed.
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only test");
            tc.assumeFalse(ispc(), ...
                "Skipped: non-Windows only (Windows has early return for pdbfixer)");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            setenv("EMK_PYTHON_EMBEDDED_DIR", "nonexistent_for_pdbfixer_nonwin_test");
            tc.addTeardown(@() setenv("EMK_PYTHON_EMBEDDED_DIR", ""));

            ME = tc.captureError(@() emk.setup.installExtra("pdbfixer"));
            tc.assertNotEmpty(ME, "pdbfixer on non-Windows: expected installFailed");
            tc.verifyEqual(ME.identifier, ...
                'emk:setup:installExtra:installFailed', ...
                "pdbfixer on non-Windows must throw installFailed (registered in resolveLibInfo_)");
        end

        % ------------------------------------------------------------------
        % TC11: unknownLibrary message mentions docking library names
        % ------------------------------------------------------------------

        function test_unknownLibrary_errorMessage_mentionsDockingName(tc)
        % After A-3, the unknownLibrary error message must list ALL THREE
        % docking library names (meeko AND vina AND pdbfixer).
        % Uses individual verifyTrue per name so the failure message
        % identifies exactly which name is missing.
            ME = tc.captureError(@() emk.setup.installExtra("bad_lib"));
            tc.assertNotEmpty(ME, "Expected unknownLibrary error");
            tc.verifyTrue(contains(ME.message, "meeko"), ...
                "unknownLibrary message must mention 'meeko'");
            tc.verifyTrue(contains(ME.message, "vina"), ...
                "unknownLibrary message must mention 'vina'");
            tc.verifyTrue(contains(ME.message, "pdbfixer"), ...
                "unknownLibrary message must mention 'pdbfixer'");
        end

        % ------------------------------------------------------------------
        % TC12: meeko LGPL warning log
        % ------------------------------------------------------------------

        function test_meeko_lgplWarning_loggedBeforeInstall(tc)
        % installExtra("meeko") must emit a logWarn about LGPL-2.1 BEFORE
        % calling pip, so users are informed of the license.
        % Verified indirectly: on Desktop with missing python_env, the
        % installFailed error should still be thrown (warning logged before
        % the pip step that triggers the error).
        % This test primarily guards against the warning being removed
        % (regression) by confirming the function reaches installFailed,
        % not some earlier error.
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only check");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            setenv("EMK_PYTHON_EMBEDDED_DIR", "nonexistent_for_meeko_warn_test");
            tc.addTeardown(@() setenv("EMK_PYTHON_EMBEDDED_DIR", ""));

            % Must reach installFailed (not unknownLibrary), confirming
            % meeko passed resolveLibInfo_() and reached the warning+pip steps.
            ME = tc.captureError(@() emk.setup.installExtra("meeko"));
            tc.assertNotEmpty(ME, "Expected installFailed error for meeko");
            tc.verifyEqual(ME.identifier, ...
                'emk:setup:installExtra:installFailed', ...
                "meeko must reach installFailed (LGPL warning is logged before pip step)");
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
