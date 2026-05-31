classdef TestInstallOnline < matlab.unittest.TestCase
% TestInstallOnline  Unit tests for src/+emk/+setup/installOnline.m
%
% These tests cover only logic paths that can be exercised without
% a real MATLAB Online environment or network access.
% Integration tests (actual pip bootstrap + rdkit install) belong
% in the smoke test suite.
%
% Coverage:
%   TC1 : Desktop environment -> notOnline error (guard fires)
%   TC1b: notOnline message contains "emk.setup.install" (user guidance)
%   TC1c: notOnline message mentions "Online" context
%   TC2 : Online environment -> notOnline NOT thrown (guard passes)
%   TC3 : Config=struct() (empty) still throws notOnline on Desktop
%         (Config param must not bypass the Desktop guard)
%   TC4 : Config with optionalLibraries still throws notOnline on Desktop
%         (Config with data must not bypass the Desktop guard)
%   TC5 : installOnline accepts Config keyword argument without error syntax
%         (arguments block regression: wrong keyword name would produce MATLAB error)
%
% Run with:
%   addpath(genpath("src"));
%   results = run(TestInstallOnline);

    methods (TestMethodSetup)
        function setupPath(tc) %#ok<MANU>
            addpath(genpath("src"));
        end
    end

    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TC1: Desktop guard fires
        % ------------------------------------------------------------------

        function test_desktopEnv_throwsNotOnline(tc)
        % installOnline() must throw notOnline when called on Desktop.
        % Verified by ensuring MATLAB_ONLINE is unset so isOnline() returns
        % false via the heuristic (Windows = not Online).
        % Skip on MATLAB Online: ismatlabonline() overrides MATLAB_ONLINE env
        % var so isOnline() cannot be mocked to return false.
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only guard test (ismatlabonline() overrides env var on MATLAB Online)");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            tc.verifyError(@() emk.setup.installOnline(), ...
                "emk:setup:installOnline:notOnline", ...
                "installOnline() must be refused on Desktop");
        end

        function test_desktopEnv_errorMessage_mentionsEmkInstall(tc)
        % Error message must guide the user to emk.setup.install() with the
        % full namespace to avoid ambiguity.
        % "emk.setup.install" must appear; matching "install()" alone is
        % insufficient because it could match other function names.
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only guard test (ismatlabonline() overrides env var on MATLAB Online)");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            ME = tc.captureError(@() emk.setup.installOnline());
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifyClass(ME, "MException", ...
                "Error must be an MException instance");
            tc.verifySubstring(lower(ME.message), "emk.setup.install", ...
                "Error message must name emk.setup.install() to guide the user");
        end

        function test_desktopEnv_errorMessage_mentionsOnline(tc)
        % Error message must reference "MATLAB Online" so the user understands
        % that installOnline() is only meaningful in that environment.
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only guard test (ismatlabonline() overrides env var on MATLAB Online)");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            ME = tc.captureError(@() emk.setup.installOnline());
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifyClass(ME, "MException", ...
                "Error must be an MException instance");
            tc.verifySubstring(lower(ME.message), "online", ...
                "Error message must mention Online context");
        end

        % ------------------------------------------------------------------
        % TC2: Online env guard passes (notOnline must NOT be thrown)
        % ------------------------------------------------------------------

        function test_onlineEnv_guardPasses_notOnlineNotThrown(tc)
        % When isOnline() returns true, the guard must PASS and notOnline
        % must NOT be thrown.
        %
        % This test runs only on actual MATLAB Online.  On Desktop,
        % overriding MATLAB_ONLINE=true bypasses the guard but lets
        % bootstrapPip_() reach system('python get-pip.py --user'), which
        % executes the Embedded Python and performs a real pip install --
        % an unacceptable side effect in a unit test suite.
        % The Desktop guard path is already covered by TC1-TC3.
        % Use emk.util.isOnline() (not ismatlabonline() directly) so that
        % the check is safe across all supported MATLAB versions.
            tc.assumeTrue(emk.util.isOnline(), ...
                "Skipped: Desktop -- overriding MATLAB_ONLINE=true would " + ...
                "trigger real pip bootstrap via Embedded Python (side effect). " + ...
                "Guard is verified by TC1 on Desktop.");

            ME = tc.captureError(@() emk.setup.installOnline());
            % Guard must not fire regardless of whether another error occurs.
            if ~isempty(ME)
                tc.verifyClass(ME, "MException", ...
                    "Any error thrown must be an MException instance");
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:setup:installOnline:notOnline", ...
                    "notOnline must never fire when isOnline() returns true");
            end
        end

        % ------------------------------------------------------------------
        % TC3: Config=struct() (empty) does not bypass Desktop guard
        % ------------------------------------------------------------------

        function test_emptyConfigStruct_desktopGuardStillFires(tc)
        % Passing Config=struct() must not bypass the Desktop guard.
        % An empty struct triggers the internal emkLoadConfig() fallback, but
        % the notOnline error must still be thrown before any pip access.
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only guard test");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            tc.verifyError(@() emk.setup.installOnline(Config=struct()), ...
                "emk:setup:installOnline:notOnline", ...
                "Config=struct() must not bypass the Desktop guard");
        end

        % ------------------------------------------------------------------
        % TC4: Config with optionalLibraries does not bypass Desktop guard
        % ------------------------------------------------------------------

        function test_configWithOptionalLibs_desktopGuardStillFires(tc)
        % Passing a non-empty Config (with optionalLibraries fields) must
        % not bypass the Desktop guard.  The guard must fire BEFORE any
        % processing of the Config struct contents.
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only guard test");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            libs = struct("pubchempy", true, "pdbfixer", true);
            cfg  = struct("optionalLibraries", libs);

            tc.verifyError(@() emk.setup.installOnline(Config=cfg), ...
                "emk:setup:installOnline:notOnline", ...
                "Config with optionalLibraries must not bypass the Desktop guard");
        end

        % ------------------------------------------------------------------
        % TC5: Config keyword argument accepted by arguments block
        % ------------------------------------------------------------------

        function test_configKeyword_acceptedByArgumentsBlock(tc)
        % installOnline(Config=cfg) must be accepted by the arguments block
        % without a MATLAB syntax error.  A wrong keyword name (e.g. 'cfg'
        % instead of 'Config') would cause MATLAB to throw a different error.
        % Verified by confirming the thrown error is notOnline (the guard),
        % not an arguments-block validation error.
        % Runs only on Desktop to avoid triggering real pip bootstrap on Online.
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only test (avoids real pip bootstrap on Online)");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            ME = tc.captureError(@() emk.setup.installOnline(Config=struct()));
            tc.assertNotEmpty(ME, "Expected notOnline error");
            tc.verifyClass(ME, "MException", ...
                "Error must be an MException (not a MATLAB argument-validation error)");
            tc.verifyEqual(string(ME.identifier), ...
                "emk:setup:installOnline:notOnline", ...
                "Config keyword must be accepted; error must be notOnline not syntax error");
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
