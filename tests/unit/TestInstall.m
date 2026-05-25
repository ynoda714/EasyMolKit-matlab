classdef TestInstall < matlab.unittest.TestCase
% TestInstall  Unit tests for src/+emk/+setup/install.m
%
% These tests cover only logic paths that can be exercised without
% network access.  Integration tests (actual download + install) belong
% in the smoke test suite.
%
% Coverage:
%   TC1 : Online environment -> notDesktop error
%   TC1b: notDesktop message mentions installOnline()
%   TC2 : Install path > 240 chars -> pathTooLong error
%   TC2b: pathTooLong message contains a path-length number
%   TC3 : Install path <= 240 chars -> pathTooLong NOT thrown
%   TC4 : Unsupported Python version -> downloadFailed (before any network call)
%   TC4b: downloadFailed message contains the unsupported version string
%
% Run with:
%   addpath(genpath("src"));
%   results = run(TestInstall);

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
        % install() must throw notDesktop when called on MATLAB Online.
        % Verified by setting MATLAB_ONLINE=true (isOnline() step 2).
            setenv("MATLAB_ONLINE", "true");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            tc.verifyError(@() emk.setup.install(), ...
                "emk:setup:install:notDesktop", ...
                "install() must be refused on MATLAB Online");
        end

        function test_onlineEnv_errorMessage_mentionsInstallOnline(tc)
        % Error message must guide user to emk.setup.installOnline().
            setenv("MATLAB_ONLINE", "true");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            ME = tc.captureError(@() emk.setup.install());
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifyClass(ME, "MException", ...
                "Error must be an MException instance");
            tc.verifySubstring(lower(ME.message), "installonline", ...
                "Error message must mention installOnline() to guide the user");
        end

        % ------------------------------------------------------------------
        % TC2: Path length > 240 chars -> pathTooLong
        % ------------------------------------------------------------------

        function test_pathTooLong_throwsPathTooLong(tc)
        % install() must throw pathTooLong when the computed install path
        % exceeds 240 characters (ADR-001 rev.3, Risk #10).
        % Path length check runs before any download attempt.
        % Skip on MATLAB Online: ismatlabonline() overrides MATLAB_ONLINE env
        % var so isOnline() cannot be mocked to return false.
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only path test (ismatlabonline() cannot be overridden by env var on MATLAB Online)");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            longDir = tc.makeLongPathDir(245);
            origDir = pwd;
            cd(longDir);
            % Teardown: LIFO order -> cd(origDir) runs first, then rmdir
            tc.addTeardown(@() tc.safeRmdir(longDir));  % added 1st, runs LAST
            tc.addTeardown(@() cd(origDir));             % added 2nd, runs FIRST

            tc.verifyError(@() emk.setup.install(), ...
                "emk:setup:install:pathTooLong", ...
                "Must throw pathTooLong when install path exceeds 240 chars");
        end

        function test_pathTooLong_errorMessage_containsLength(tc)
        % Error message must report the path length (3-digit number) for
        % actionable debugging.
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only path test (ismatlabonline() cannot be overridden by env var on MATLAB Online)");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            longDir = tc.makeLongPathDir(245);
            origDir = pwd;
            cd(longDir);
            tc.addTeardown(@() tc.safeRmdir(longDir));
            tc.addTeardown(@() cd(origDir));

            ME = tc.captureError(@() emk.setup.install());
            tc.assertNotEmpty(ME, "Expected a pathTooLong error");
            tc.verifyClass(ME, "MException", ...
                "Error must be an MException instance");
            tc.verifyTrue(~isempty(regexp(ME.message, '\d{3}', 'once')), ...
                "pathTooLong message must contain the path length (3-digit number)");
        end

        % ------------------------------------------------------------------
        % TC3: Path length <= 240 chars -> pathTooLong NOT thrown
        % ------------------------------------------------------------------

        function test_shortPath_doesNotThrowPathTooLong(tc)
        % For a short path, install() must NOT throw pathTooLong.
        % It may throw downloadFailed or another error, but NOT pathTooLong.
        % Calling with an unsupported version ensures we get a deterministic
        % early failure without any network access.
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            tmpDir  = tc.makeTempDir();
            origDir = pwd;
            cd(tmpDir);
            tc.addTeardown(@() tc.safeRmdir(tmpDir));
            tc.addTeardown(@() cd(origDir));

            ME = tc.captureError(@() emk.setup.install(PythonVersion="3.11"));
            if ~isempty(ME)
                tc.verifyClass(ME, "MException", ...
                    "Any error thrown must be an MException instance");
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:setup:install:pathTooLong", ...
                    "A short path must never trigger pathTooLong");
            end
        end

        % ------------------------------------------------------------------
        % TC4: Unsupported Python version -> downloadFailed before network
        % ------------------------------------------------------------------

        function test_unsupportedVersion_throwsDownloadFailed(tc)
        % resolvePatchVersion_() must throw downloadFailed immediately for
        % unsupported versions (e.g. "3.11"), before any network activity.
        % Verified by: the error fires even in a tmpDir with no network.
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only version test (ismatlabonline() cannot be overridden by env var on MATLAB Online)");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            tmpDir  = tc.makeTempDir();
            origDir = pwd;
            cd(tmpDir);
            tc.addTeardown(@() tc.safeRmdir(tmpDir));
            tc.addTeardown(@() cd(origDir));

            tc.verifyError(@() emk.setup.install(PythonVersion="3.11"), ...
                "emk:setup:install:downloadFailed", ...
                "Unsupported Python version must throw downloadFailed");
        end

        function test_unsupportedVersion_errorMessage_containsVersion(tc)
        % Error message must include the unsupported version string so the
        % user knows exactly which version was rejected.
            tc.assumeFalse(emk.util.isOnline(), ...
                "Skipped: Desktop-only version test (ismatlabonline() cannot be overridden by env var on MATLAB Online)");
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));

            tmpDir  = tc.makeTempDir();
            origDir = pwd;
            cd(tmpDir);
            tc.addTeardown(@() tc.safeRmdir(tmpDir));
            tc.addTeardown(@() cd(origDir));

            ME = tc.captureError(@() emk.setup.install(PythonVersion="3.11"));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifyClass(ME, "MException", ...
                "Error must be an MException instance");
            tc.verifySubstring(ME.message, "3.11", ...
                "Error message must contain the rejected version string '3.11'");
        end

    end

    % ======================================================================
    methods (Access = private)

        function longDir = makeLongPathDir(~, targetTotalLen)
        % Create a temporary directory such that:
        %   length(pwd) + length("\python_env\Lib\site-packages\rdkit")
        %   >= targetTotalLen
        % This simulates the MAX_PATH risk condition checked by install().
        %
        % Suffix components: \python_env\Lib\site-packages\rdkit = 35 chars
            suffixLen = 35;
            base = strtrim(tempdir);
            if ~isempty(base) && base(end) == filesep
                base = base(1:end-1);
            end
            needed  = max(targetTotalLen - length(base) - suffixLen, 10);
            dirName = repmat('a', 1, needed);
            longDir = fullfile(base, dirName);
            if ~isfolder(longDir)
                mkdir(longDir);
            end
        end

        function tmpDir = makeTempDir(~)
        % Create a short, unique temporary directory for test isolation.
            tmpDir = fullfile(tempdir, "emk_inst_" + ...
                char(datetime("now", "Format", "yyyyMMdd_HHmmssSSS")));
            if ~isfolder(tmpDir)
                mkdir(tmpDir);
            end
        end

        function safeRmdir(~, d)
        % Remove directory if it still exists (teardown helper).
            if isfolder(d)
                rmdir(d, "s");
            end
        end

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
