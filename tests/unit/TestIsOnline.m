classdef TestIsOnline < matlab.unittest.TestCase
% TestIsOnline  Unit tests for src/+emk/+util/isOnline.m
%
% Run with:
%   addpath(genpath("src"));
%   results = run(TestIsOnline);
%
% Note: M0-9 verification on actual MATLAB Online is deferred to M1 smoke tests.
% These tests cover the Desktop path and env-var override path.

    methods (TestMethodSetup)
        function setupPath(tc) %#ok<MANU>
            addpath(genpath("src"));
        end
    end

    methods (Test)
        % ------------------------------------------------------------------
        function test_returnType_isLogical(tc)
        % isOnline always returns a logical scalar
            result = emk.util.isOnline();
            tc.verifyClass(result, "logical", "isOnline must return logical");
            tc.verifySize(result, [1 1], "isOnline must return a [1 1] scalar");
        end

        % ------------------------------------------------------------------
        function test_envVar_true_returnsTrue(tc)
        % MATLAB_ONLINE=true env var causes isOnline to return true
        % (simulates the Online environment on Desktop)
            setenv("MATLAB_ONLINE", "true");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            tc.verifyTrue(emk.util.isOnline(), ...
                "isOnline must return true when MATLAB_ONLINE=true");
        end

        % ------------------------------------------------------------------
        function test_envVar_1_returnsTrue(tc)
        % MATLAB_ONLINE=1 env var causes isOnline to return true
            setenv("MATLAB_ONLINE", "1");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            tc.verifyTrue(emk.util.isOnline());
        end

        % ------------------------------------------------------------------
        function test_envVar_false_returnsFalse(tc)
        % MATLAB_ONLINE=false returns false
            setenv("MATLAB_ONLINE", "false");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            tc.verifyFalse(emk.util.isOnline());
        end

        % ------------------------------------------------------------------
        function test_envVar_empty_usesNativeDetection(tc)
        % When MATLAB_ONLINE is unset, detection falls back to native check.
        % On Windows Desktop with ismatlabonline available, result must be false.
            setenv("MATLAB_ONLINE", "");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            result = emk.util.isOnline();
            if ispc()
                % Windows Desktop must never be detected as Online
                tc.verifyFalse(result, ...
                    "Windows Desktop must not be detected as MATLAB Online");
            else
                % Non-Windows: result depends on actual environment
                tc.verifyClass(result, "logical");
            end
        end

        % ------------------------------------------------------------------
        function test_envVar_True_caseInsensitive(tc)
        % MATLAB_ONLINE=TRUE (uppercase) is treated as true
            setenv("MATLAB_ONLINE", "TRUE");
            tc.addTeardown(@() setenv("MATLAB_ONLINE", ""));
            tc.verifyTrue(emk.util.isOnline());
        end
    end
end
