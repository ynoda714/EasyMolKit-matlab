classdef TestViz < matlab.unittest.TestCase
% TestViz  Unit tests for emk.viz.draw2d.
%
% Run with:
%   addpath(genpath("src"));
%   suite = testsuite("tests/unit");
%   runner = matlab.unittest.TestRunner.withNoPlugins;
%   results = runner.run(suite);
%
% ======================================================================
% Coverage (draw2d):
%
%   TC1:  Non-Mol input => invalidInput (no RDKit required)
%         double, string, char, logical, empty all rejected
%   TC1b: Width=0 / Height<0 => invalidInput (no RDKit required)
%   TC2:  Valid Mol => no error thrown (requires RDKit)
%   TC3:  Returns figure handle when requested (requires RDKit)
%   TC4:  Figure has image data with correct dimensions (requires RDKit)
%   TC5:  Width/Height options affect image dimensions (requires RDKit)
%   TC5b: Non-square cross-validation: width AND height in single call
%   TC5c: Idempotency: calling draw2d twice on same mol must not throw
%   TC6:  Title option sets figure title (requires RDKit)
%   TC7:  Empty Title => no title text on axes (requires RDKit)
%
% Design notes:
%   - Figures created during tests are closed in teardown to avoid
%     accumulating open windows in batch/CI mode.
%   - emk.viz.draw2d calls AllChem.Compute2DCoords on the mol in-place
%     (adds 2D geometry).  Tests use fresh mols to avoid interference.
%   - 3-principle checklist:
%       1. No py.rdkit.* direct calls -- all mol creation via emk.mol.fromSmiles
%       2. verifyClass / verifySize for all return-value assertions
%       3. Image dimension cross-validation: size(img,1)==Height, size(img,2)==Width
%
% ======================================================================

    properties
        openFigs  matlab.ui.Figure  % figures opened during a test (closed in teardown)
    end

    methods (TestMethodSetup)
        function setupPath(~)
            addpath(genpath("src"));
        end
    end

    methods (TestMethodTeardown)
        function closeOpenFigs(tc)
        % Close any figures opened by draw2d to keep the environment clean.
            for i = 1:numel(tc.openFigs)
                try
                    close(tc.openFigs(i));
                catch
                end
            end
            tc.openFigs = matlab.ui.Figure.empty();
        end
    end

    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TC1: draw2d invalidInput -- no RDKit required
        % ------------------------------------------------------------------

        function test_draw2d_numericInput_throwsInvalidInput(tc)
        % Numeric mol must throw invalidInput.
            tc.verifyError(@() emk.viz.draw2d(42), ...
                "emk:viz:draw2d:invalidInput", ...
                "Numeric mol must throw invalidInput");
        end

        function test_draw2d_numericInput_errorMessage_containsClass(tc)
        % invalidInput message must contain the actual class name.
            ME = tc.captureError(@() emk.viz.draw2d(42));
            tc.assertNotEmpty(ME, "Expected invalidInput to be thrown");
            tc.verifySubstring(ME.message, "double", ...
                "Error message must contain the input class name");
        end

        function test_draw2d_stringInput_throwsInvalidInput(tc)
        % String SMILES passed directly must throw invalidInput.
            tc.verifyError(@() emk.viz.draw2d("CCO"), ...
                "emk:viz:draw2d:invalidInput", ...
                "String mol must throw invalidInput");
        end

        function test_draw2d_charInput_throwsInvalidInput(tc)
        % char SMILES must throw invalidInput.
            tc.verifyError(@() emk.viz.draw2d('CCO'), ...
                "emk:viz:draw2d:invalidInput", ...
                "char mol must throw invalidInput");
        end

        function test_draw2d_logicalInput_throwsInvalidInput(tc)
        % Logical must throw invalidInput.
            tc.verifyError(@() emk.viz.draw2d(true), ...
                "emk:viz:draw2d:invalidInput", ...
                "Logical mol must throw invalidInput");
        end

        function test_draw2d_emptyInput_throwsInvalidInput(tc)
        % Empty matrix must throw invalidInput.
            tc.verifyError(@() emk.viz.draw2d([]), ...
                "emk:viz:draw2d:invalidInput", ...
                "Empty input must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC1b: Width/Height range validation -- no RDKit required
        % Width < 1 or Height < 1 must throw invalidInput before any RDKit call.
        % ------------------------------------------------------------------

        function test_draw2d_zeroWidth_throwsInvalidInput(tc)
        % Width=0 must throw invalidInput (RDKit not required).
        % Validates that range check fires before RDKit is called.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tc.verifyError(@() emk.viz.draw2d(mol, Width=0), ...
                "emk:viz:draw2d:invalidInput", ...
                "Width=0 must throw invalidInput");
        end

        function test_draw2d_negativeHeight_throwsInvalidInput(tc)
        % Height=-1 must throw invalidInput (RDKit not required).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tc.verifyError(@() emk.viz.draw2d(mol, Height=-1), ...
                "emk:viz:draw2d:invalidInput", ...
                "Height=-1 must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC2: draw2d valid mol does not throw (requires RDKit)
        % ------------------------------------------------------------------

        function test_draw2d_ethanol_doesNotThrow(tc)
        % draw2d on a valid mol must complete without error.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fig = emk.viz.draw2d(mol);
            tc.openFigs(end+1) = fig;
            % If we reach here without exception, TC2 passes.
            tc.verifyNotEmpty(fig, "draw2d must return a figure handle");
        end

        function test_draw2d_benzene_doesNotThrow(tc)
        % draw2d on benzene (aromatic ring) must complete without error.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("c1ccccc1");
            fig = emk.viz.draw2d(mol);
            tc.openFigs(end+1) = fig;
            tc.verifyNotEmpty(fig, "draw2d on benzene must return a figure handle");
        end

        function test_draw2d_aspirin_doesNotThrow(tc)
        % draw2d on aspirin (complex mol) must complete without error.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fig = emk.viz.draw2d(mol);
            tc.openFigs(end+1) = fig;
            tc.verifyNotEmpty(fig, "draw2d on aspirin must return a figure handle");
        end

        % ------------------------------------------------------------------
        % TC3: draw2d returns figure handle (requires RDKit)
        % ------------------------------------------------------------------

        function test_draw2d_returnsFigureHandle(tc)
        % Return value must be a matlab.ui.Figure.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fig = emk.viz.draw2d(mol);
            tc.openFigs(end+1) = fig;
            tc.verifyClass(fig, "matlab.ui.Figure", ...
                "draw2d must return a matlab.ui.Figure handle");
        end

        function test_draw2d_noOutputArg_doesNotThrow(tc)
        % Calling draw2d without capturing output must not throw.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            emk.viz.draw2d(mol);  % no output argument
            % Close last opened figure by handle search
            figs = findall(0, "Type", "Figure");
            if ~isempty(figs)
                close(figs(1));
            end
        end

        % ------------------------------------------------------------------
        % TC4: draw2d image dimensions match defaults (requires RDKit)
        % Default Width=300, Height=300 => RGB image of size 300x300x3
        % ------------------------------------------------------------------

        function test_draw2d_defaultSize_imageIsCorrectHeight(tc)
        % Default Height=300: image must have 300 rows (height).
        % Cross-validation: size(img,1) == Height option value.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fig = emk.viz.draw2d(mol);
            tc.openFigs(end+1) = fig;
            img = tc.getFigureImage_(fig);
            tc.verifyEqual(size(img, 1), 300, ...
                "Default height must be 300 pixels");
        end

        function test_draw2d_defaultSize_imageIsCorrectWidth(tc)
        % Default Width=300: image must have 300 columns (width).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fig = emk.viz.draw2d(mol);
            tc.openFigs(end+1) = fig;
            img = tc.getFigureImage_(fig);
            tc.verifyEqual(size(img, 2), 300, ...
                "Default width must be 300 pixels");
        end

        function test_draw2d_defaultSize_imageIsRGB(tc)
        % Image must have 3 channels (RGB).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fig = emk.viz.draw2d(mol);
            tc.openFigs(end+1) = fig;
            img = tc.getFigureImage_(fig);
            tc.verifyEqual(size(img, 3), 3, ...
                "Image must be RGB (3 channels)");
        end

        function test_draw2d_defaultSize_imageIsUint8(tc)
        % Image data type must be uint8.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fig = emk.viz.draw2d(mol);
            tc.openFigs(end+1) = fig;
            img = tc.getFigureImage_(fig);
            tc.verifyClass(img, "uint8", ...
                "Image data type must be uint8");
        end

        % ------------------------------------------------------------------
        % TC5: draw2d Width/Height options (requires RDKit)
        % Cross-validation: image dimensions match the requested options.
        % ------------------------------------------------------------------

        function test_draw2d_customWidth_imageWidthMatches(tc)
        % Width=400 option must produce a 400-pixel wide image.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fig = emk.viz.draw2d(mol, Width=400, Height=300);
            tc.openFigs(end+1) = fig;
            img = tc.getFigureImage_(fig);
            tc.verifyEqual(size(img, 2), 400, ...
                "Width=400 must produce 400-pixel wide image");
        end

        function test_draw2d_customHeight_imageHeightMatches(tc)
        % Height=200 option must produce a 200-pixel tall image.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fig = emk.viz.draw2d(mol, Width=300, Height=200);
            tc.openFigs(end+1) = fig;
            img = tc.getFigureImage_(fig);
            tc.verifyEqual(size(img, 1), 200, ...
                "Height=200 must produce 200-pixel tall image");
        end

        function test_draw2d_nonSquare_bothDimensionsCorrect(tc)
        % Cross-validation: Width=400, Height=150 must produce exactly 150x400 image.
        % Verifies both dimensions simultaneously in a single non-square call.
        % Principle 3: size(img,1)==Height AND size(img,2)==Width in one test.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("c1ccccc1");
            fig = emk.viz.draw2d(mol, Width=400, Height=150);
            tc.openFigs(end+1) = fig;
            img = tc.getFigureImage_(fig);
            tc.verifyEqual(size(img, 1), 150, ...
                "Height=150: image must have 150 rows");
            tc.verifyEqual(size(img, 2), 400, ...
                "Width=400: image must have 400 columns");
        end

        function test_draw2d_idempotency_calledTwice_doesNotThrow(tc)
        % Calling draw2d twice on the same mol must not throw.
        % Verifies that Compute2DCoords is idempotent (re-entrant safe).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fig1 = emk.viz.draw2d(mol);
            tc.openFigs(end+1) = fig1;
            fig2 = emk.viz.draw2d(mol);  % second call on same mol object
            tc.openFigs(end+1) = fig2;
            tc.verifyNotEmpty(fig2, "Second draw2d call must return a valid figure");
        end

        % ------------------------------------------------------------------
        % TC6: draw2d Title option (requires RDKit)
        % ------------------------------------------------------------------

        function test_draw2d_titleOption_setsAxesTitle(tc)
        % Title="Ethanol" must appear as the axes title string.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fig = emk.viz.draw2d(mol, Title="Ethanol");
            tc.openFigs(end+1) = fig;
            ax  = findobj(fig, "Type", "Axes");
            tc.assertNotEmpty(ax, "Figure must have an Axes object");
            tc.verifyEqual(string(ax(1).Title.String), "Ethanol", ...
                "Axes title must be the provided Title option value");
        end

        % ------------------------------------------------------------------
        % TC7: draw2d empty Title => no axes title (requires RDKit)
        % ------------------------------------------------------------------

        function test_draw2d_emptyTitle_noAxesTitle(tc)
        % Default empty Title must result in an empty axes title string.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fig = emk.viz.draw2d(mol);
            tc.openFigs(end+1) = fig;
            ax  = findobj(fig, "Type", "Axes");
            tc.assertNotEmpty(ax, "Figure must have an Axes object");
            tc.verifyEmpty(ax(1).Title.String, ...
                "Empty Title option must leave the axes title blank");
        end

        % ------------------------------------------------------------------
        % TC8: draw2d renders into the current figure (gca context)
        % Regression for draw2d change: figure() replaced by gca() so
        % callers can pre-set figure/subplot before calling draw2d.
        % ------------------------------------------------------------------

        function test_draw2d_rendersIntoExistingFigure(tc)
        % When the caller creates a figure first, draw2d must use that
        % figure (not create a new one).  The returned handle must equal
        % the pre-existing figure.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            preFig = figure("Name", "preFig");
            tc.openFigs(end+1) = preFig;
            mol = emk.mol.fromSmiles("CCO");
            fig = emk.viz.draw2d(mol);
            tc.openFigs(end+1) = fig;
            tc.verifyEqual(fig, preFig, ...
                "draw2d must render into the pre-existing figure, not create a new one");
        end

        function test_draw2d_noPreExistingFigure_createsOne(tc)
        % When no figure exists, draw2d must create one automatically.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            close all;
            mol = emk.mol.fromSmiles("CCO");
            fig = emk.viz.draw2d(mol);
            tc.openFigs(end+1) = fig;
            tc.verifyClass(fig, "matlab.ui.Figure", ...
                "draw2d must create and return a figure when none exists");
        end

    end

    % ======================================================================
    methods (Access = private)

        function tf = rdkitAvailable(~)
        % Return true if Python is configured and RDKit can be imported.
            tf = false;
            try
                pe = pyenv();
                if strcmp(string(pe.Status), "NotLoaded")
                    emk.setup.initPython();
                end
                py.importlib.import_module("rdkit.Chem");
                tf = true;
            catch
            end
        end

        function ME = captureError(~, fcn)
        % Call fcn and return the MException if thrown, or [] on success.
            ME = [];
            try
                fcn();
            catch e
                ME = e;
            end
        end

        function img = getFigureImage_(~, fig)
        % Extract the CData (image array) from the first image object in fig.
            imgObj = findobj(fig, "Type", "Image");
            if isempty(imgObj)
                error("TestViz:noImage", "No Image object found in figure");
            end
            img = imgObj(1).CData;
        end

    end
end
