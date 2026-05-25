classdef TestMakeRunDir < matlab.unittest.TestCase
% TestMakeRunDir  Unit tests for src/util/makeRunDir.m
%
% Run with:
%   addpath(genpath("src"));
%   results = run(TestMakeRunDir);

    properties
        TempBase  % temporary base directory for test isolation
    end

    methods (TestMethodSetup)
        function setupPath(tc)
            addpath(genpath("src"));
            tc.TempBase = fullfile(tempdir, "emk_test_makeRunDir_" + ...
                char(datetime("now", "Format", "yyyyMMdd_HHmmssSSS")));
            mkdir(tc.TempBase);
        end
    end

    methods (TestMethodTeardown)
        function cleanupTempDir(tc)
            if isfolder(tc.TempBase)
                rmdir(tc.TempBase, "s");
            end
        end
    end

    methods (Test)
        % ------------------------------------------------------------------
        function test_defaultCall_createsDirUnderResultRuns(tc)
        % makeRunDir() creates a directory under the specified BaseDir
            runDir = makeRunDir("BaseDir", tc.TempBase);
            tc.verifyTrue(isfolder(runDir), "run directory must exist after makeRunDir");
            tc.verifySubstring(runDir, tc.TempBase);
        end

        % ------------------------------------------------------------------
        function test_prefix_appendedToTimestamp(tc)
        % Prefix option is appended after timestamp in directory name
            runDir = makeRunDir("Prefix", "smoke_mol", "BaseDir", tc.TempBase);
            [~, dirName] = fileparts(runDir);
            tc.verifySubstring(dirName, "smoke_mol");
        end

        % ------------------------------------------------------------------
        function test_noPrefix_dirNameIsTimestampOnly(tc)
        % Without Prefix, directory name matches yyyyMMdd_HHmmss pattern
            runDir = makeRunDir("BaseDir", tc.TempBase);
            [~, dirName] = fileparts(runDir);
            tc.verifyMatches(dirName, "^\d{8}_\d{6}$");
        end

        % ------------------------------------------------------------------
        function test_returnedPath_isString(tc)
        % makeRunDir must return a MATLAB string scalar.
        % fullfile() with a string baseDir always produces a string,
        % so char is not an expected return type here.
            runDir = makeRunDir("BaseDir", tc.TempBase);
            tc.verifyClass(runDir, "string", ...
                "makeRunDir must return a MATLAB string");
            tc.verifyTrue(isStringScalar(runDir), ...
                "makeRunDir must return a scalar string");
        end

        % ------------------------------------------------------------------
        function test_twoCalls_createDistinctDirs(tc)
        % Two successive calls create two distinct directories
        %   (relies on at least 1-second apart or unique names)
            runDir1 = makeRunDir("Prefix", "a", "BaseDir", tc.TempBase);
            runDir2 = makeRunDir("Prefix", "b", "BaseDir", tc.TempBase);
            tc.verifyNotEqual(runDir1, runDir2);
            tc.verifyTrue(isfolder(runDir1));
            tc.verifyTrue(isfolder(runDir2));
        end
    end
end
