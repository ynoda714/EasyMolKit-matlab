classdef TestUseExternal < matlab.unittest.TestCase
% TestUseExternal  Unit tests for src/+emk/+setup/useExternal.m
%
% These tests cover input validation and idempotency without actually
% invoking pyenv with a real Python binary.  Integration tests (real
% pyenv switch + validate) belong in the smoke test suite.
%
% Coverage:
%   TC1 : Non-string/char input (double) -> invalidInput error
%   TC1b: Cell array input -> invalidInput error
%   TC1c: Empty string "" -> invalidInput error
%   TC1d: Whitespace-only string "   " -> invalidInput error (strtrim guard)
%   TC1e: Empty char '' -> invalidInput error
%   TC2 : Non-existent file path -> fileNotFound error
%   TC2b: fileNotFound message contains the invalid path
%   TC2c: fileNotFound message advises user on the cause
%   TC3 : char path input is accepted (coerced to string) -> fileNotFound
%   TC4 : Python already Loaded -> no exception thrown (idempotency guard)
%
% Notes:
%   TC4 requires Python to be already loaded.  It is skipped gracefully
%   (assumeTrue) when Python is in NotLoaded state.
%
% Run with:
%   addpath(genpath("src")); addpath(genpath("tests"));
%   results = run(TestUseExternal);

    methods (TestMethodSetup)
        function setupPath(tc) %#ok<MANU>
            addpath(genpath("src"));
        end
    end

    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TC1: Input type / content validation
        % ------------------------------------------------------------------

        function test_numericInput_throwsInvalidInput(tc)
        % A numeric value (not a path) must be rejected immediately with
        % invalidInput before any file check or pyenv call.
            tc.verifyError(@() emk.setup.useExternal(42), ...
                "emk:setup:useExternal:invalidInput", ...
                "Numeric input must throw invalidInput");
        end

        function test_cellInput_throwsInvalidInput(tc)
        % A cell array must be rejected; MATLAB arguments do not coerce cell
        % to string automatically for scalar checks.
            tc.verifyError(@() emk.setup.useExternal({"path"}), ...
                "emk:setup:useExternal:invalidInput", ...
                "Cell array input must throw invalidInput");
        end

        function test_emptyString_throwsInvalidInput(tc)
        % An empty MATLAB string "" must be rejected.
        % Uses strlength (not isempty) to detect empty string.
            tc.verifyError(@() emk.setup.useExternal(""), ...
                "emk:setup:useExternal:invalidInput", ...
                "Empty string must throw invalidInput");
        end

        function test_whitespaceOnlyString_throwsInvalidInput(tc)
        % A whitespace-only string ("   ") must be rejected after strtrim.
        % Ensures the guard is not fooled by blank paths.
            tc.verifyError(@() emk.setup.useExternal("   "), ...
                "emk:setup:useExternal:invalidInput", ...
                "Whitespace-only string must throw invalidInput after strtrim");
        end

        function test_emptyChar_throwsInvalidInput(tc)
        % An empty char '' must be rejected.
        % The manual char-check path must handle '' as well as "".
            tc.verifyError(@() emk.setup.useExternal(''), ...
                "emk:setup:useExternal:invalidInput", ...
                "Empty char must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC2: Non-existent file path
        % ------------------------------------------------------------------

        function test_nonexistentPath_throwsFileNotFound(tc)
        % A syntactically valid but non-existent path must throw fileNotFound.
        % The idempotency guard must be bypassed (Python not loaded) for this
        % check; if Python is Loaded the function silently returns first.
        % Test is unconditional: if Python is Loaded the guard fires and the
        % function returns without checking the file (which is also OK).
            fakeExe = fullfile(tempdir, "nonexistent_python_exe_12345.exe");
            if strcmp(string(pyenv().Status), "NotLoaded")
                tc.verifyError(@() emk.setup.useExternal(fakeExe), ...
                    "emk:setup:useExternal:fileNotFound", ...
                    "Non-existent executable must throw fileNotFound (Python NotLoaded)");
            else
                % Python already Loaded: guard fires, no throw expected
                ME = tc.captureError(@() emk.setup.useExternal(fakeExe));
                tc.verifyEmpty(ME, ...
                    "When Python is Loaded, useExternal must not throw (idempotency guard)");
            end
        end

        function test_nonexistentPath_errorMessage_containsPath(tc)
        % fileNotFound error message must contain the invalid path so the
        % user knows which executable was not found.
            tc.assumeTrue(strcmp(string(pyenv().Status), "NotLoaded"), ...
                "Skipped: Python already loaded; fileNotFound guard would not fire");
            fakeExe = fullfile(tempdir, "nonexistent_python_exe_99999.exe");

            ME = tc.captureError(@() emk.setup.useExternal(fakeExe));
            tc.assertNotEmpty(ME, "Expected fileNotFound error");
            tc.verifyClass(ME, "MException", "Error must be MException");
            tc.verifySubstring(ME.message, char(fakeExe), ...
                "fileNotFound message must contain the invalid path");
        end

        function test_nonexistentPath_errorMessage_guidesUser(tc)
        % fileNotFound error message must mention a Python installation
        % so the user understands what kind of path is expected.
            tc.assumeTrue(strcmp(string(pyenv().Status), "NotLoaded"), ...
                "Skipped: Python already loaded; fileNotFound guard would not fire");
            fakeExe = fullfile(tempdir, "nonexistent_python_exe_77777.exe");

            ME = tc.captureError(@() emk.setup.useExternal(fakeExe));
            tc.assertNotEmpty(ME, "Expected fileNotFound error");
            hasInstall = contains(lower(ME.message), "python") || ...
                         contains(lower(ME.message), "install");
            tc.verifyTrue(hasInstall, ...
                "fileNotFound message must mention Python installation context");
        end

        % ------------------------------------------------------------------
        % TC3: char input coercion
        % ------------------------------------------------------------------

        function test_charPath_acceptedAsString(tc)
        % useExternal() must accept char input (manual char check converts
        % it to string internally).  A nonexistent path is used for a
        % deterministic early failure without real pyenv interaction.
        % If Python is Loaded, the idempotency guard fires first (no throw).
            fakeExe = fullfile(tempdir, "nonexistent_char_python.exe");
            ME = tc.captureError(@() emk.setup.useExternal(char(fakeExe)));
            if ~isempty(ME)
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:setup:useExternal:invalidInput", ...
                    "char input must NOT trigger invalidInput; it should be accepted");
            end
        end

        % ------------------------------------------------------------------
        % TC4: Idempotency guard (Python already Loaded)
        % ------------------------------------------------------------------

        function test_alreadyLoaded_noThrow(tc)
        % When Python is already loaded, useExternal() must return silently
        % (with a warning via logWarn) rather than throwing.
        % The pyenv Version/ExecutionMode cannot be changed once Python is
        % loaded, so this silent-skip behaviour is the intended contract.
            tc.assumeTrue(~strcmp(string(pyenv().Status), "NotLoaded"), ...
                "Skipped: Python is not loaded; cannot test loaded-state guard");

            fakeExe = fullfile(tempdir, "some_fake_python.exe");
            ME = tc.captureError(@() emk.setup.useExternal(fakeExe));
            tc.verifyEmpty(ME, ...
                "useExternal() must not throw when Python is already loaded");
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
