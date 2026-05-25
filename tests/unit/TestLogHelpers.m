classdef TestLogHelpers < matlab.unittest.TestCase
% TestLogHelpers  Unit tests for src/util/log*.m helper functions
%
% Run with:
%   addpath(genpath("src"));
%   results = run(TestLogHelpers);

    methods (TestMethodSetup)
        function setupPath(tc) %#ok<MANU>
            addpath(genpath("src"));
        end
    end

    methods (Test)
        % ------------------------------------------------------------------
        function test_logInfo_producesOutput(tc)
        % logInfo writes a line containing [INFO] to stdout
            captured = evalc("logInfo(""hello world"")");
            tc.verifySubstring(captured, "[INFO]");
            tc.verifySubstring(captured, "hello world");
        end

        % ------------------------------------------------------------------
        function test_logWarn_producesOutput(tc)
        % logWarn writes a line containing [WARN] to stdout
            captured = evalc("logWarn(""something suspicious"")");
            tc.verifySubstring(captured, "[WARN]");
            tc.verifySubstring(captured, "something suspicious");
        end

        % ------------------------------------------------------------------
        function test_logError_producesOutput(tc)
        % logError writes a line containing [ERROR] to stdout
            captured = evalc("logError(""bad input"")");
            tc.verifySubstring(captured, "[ERROR]");
            tc.verifySubstring(captured, "bad input");
        end

        % ------------------------------------------------------------------
        function test_logInfo_formatArgs(tc)
        % logInfo accepts sprintf-style format arguments
            captured = evalc("logInfo(""%d items processed"", 42)");
            tc.verifySubstring(captured, "42 items processed");
        end

        % ------------------------------------------------------------------
        function test_logDebug_suppressedByDefault(tc)
        % logDebug produces no output unless EMK_LOG_VERBOSE=1
            setenv("EMK_LOG_VERBOSE", "");
            tc.addTeardown(@() setenv("EMK_LOG_VERBOSE", ""));
            captured = evalc("logDebug(""verbose detail"")");
            tc.verifyEmpty(strtrim(captured), "logDebug must be silent when EMK_LOG_VERBOSE is not set");
        end

        % ------------------------------------------------------------------
        function test_logDebug_outputWhenVerbose(tc)
        % logDebug writes [DEBUG] when EMK_LOG_VERBOSE=1
            setenv("EMK_LOG_VERBOSE", "1");
            tc.addTeardown(@() setenv("EMK_LOG_VERBOSE", ""));
            captured = evalc("logDebug(""verbose detail"")");
            tc.verifySubstring(captured, "[DEBUG]");
            tc.verifySubstring(captured, "verbose detail");
        end

        % ------------------------------------------------------------------
        function test_logProgress_finalStepHasNewline(tc)
        % logProgress(n, n, label) shows 100% and label
            captured = evalc("logProgress(5, 5, ""done"")");
            tc.verifySubstring(captured, "100%");
            tc.verifySubstring(captured, "done");
        end

        % ------------------------------------------------------------------
        function test_logProgress_partialStep(tc)
        % logProgress shows correct percentage for partial progress
            captured = evalc("logProgress(1, 4, ""step"")");
            tc.verifySubstring(captured, "25%");
        end

        % ------------------------------------------------------------------
        function test_logProgress_n1_shows100(tc)
        % logProgress(1, 1, label) immediately shows 100%
            captured = evalc("logProgress(1, 1, ""single"")");
            tc.verifySubstring(captured, "100%");
        end

        % ------------------------------------------------------------------
        function test_logInfo_timestampFormat(tc)
        % Output line starts with [HH:MM:SS] timestamp pattern
            captured = evalc("logInfo(""ts check"")");
            tc.verifyMatches(strtrim(captured), "^\[\d{2}:\d{2}:\d{2}\]\[INFO\].*");
        end

        % ------------------------------------------------------------------
        function test_logInfo_stringContainingPercent_doesNotThrow(tc)
        % logInfo("%s", str) must not throw when str contains "%" characters.
        % Regression for G15: logInfo(str) where str contains "%" would fail
        % because the internal sprintf(msg, varargin{:}) treats "%" as a
        % format specifier.  Callers must always use logInfo("%s", str).
            captured = evalc("logInfo(""%s"", ""50% done"")");
            tc.verifySubstring(captured, "50% done", ...
                "logInfo with percent sign in string must output correctly");
        end

        % ------------------------------------------------------------------
        function test_logWarn_stringContainingPercent_doesNotThrow(tc)
        % logWarn("%s", str) must not throw when str contains "%" characters.
            captured = evalc("logWarn(""%s"", ""3.5% error rate"")");
            tc.verifySubstring(captured, "3.5% error rate", ...
                "logWarn with percent sign in string must output correctly");
        end
    end
end
