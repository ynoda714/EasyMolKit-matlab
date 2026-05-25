classdef TestRecipe < matlab.unittest.TestCase
% TestRecipe  Unit tests for src/+emk/+setup/recipe.m
%
% recipe() displays installation instructions in the Command Window and
% returns no value.  Tests verify error handling for unknown names and
% confirm that all supported library names produce no exception.
%
% Coverage:
%   TC1 : Unknown library name -> unknownLibrary error
%   TC1b: unknownLibrary message contains the bad library name
%   TC1c: unknownLibrary message lists at least one supported name
%   TC2 : 'pubchempy' -> no exception (Track 1 recipe)
%   TC3 : 'mordred'   -> no exception (Track 1 recipe)
%   TC4 : 'biopython' -> no exception (Track 1 recipe)
%   TC5 : 'openbabel' -> no exception (Track 2 recipe + GPL note)
%   TC6 : 'mdanalysis'-> no exception (Track 2 recipe + GPL note)
%   TC7 : 'pymol'     -> no exception (Track 2 recipe)
%   TC8 : Canonical name list is complete (all known names are covered)
%   TC9 : char input accepted
%   TC10: 'meeko'     -> no exception (Track 1 recipe, A-3)
%   TC10b: 'meeko' recipe mentions LGPL-2.1 and compliance.md/CL-7
%   TC10c: 'meeko' recipe identifies Track 1 classification
%   TC11: 'vina'      -> no exception (Track 1 recipe, A-3)
%   TC12: 'pdbfixer'  -> no exception (Track 1 recipe, A-3)
%   TC12b: 'pdbfixer' recipe mentions openmm dependency (in printPdbfixer_ body)
%   TC12c: 'pdbfixer' recipe on Windows mentions Smart App Control (SAC)
%   TC12d: 'pdbfixer' recipe mentions MATLAB Online as workaround path
%   TC12e: 'pdbfixer' recipe mentions cfg.optionalLibraries flag for Online install
%   TC13: 'docking'   -> no exception (combined Track 1 recipe, A-3)
%   TC13b: 'docking' recipe mentions all three components (meeko, vina, pdbfixer)
%
% Run with:
%   addpath(genpath("src")); addpath(genpath("tests"));
%   results = run(TestRecipe);

    methods (TestMethodSetup)
        function setupPath(tc) %#ok<MANU>
            addpath(genpath("src"));
        end
    end

    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TC1: Unknown library name
        % ------------------------------------------------------------------

        function test_unknownLibrary_throwsUnknownLibrary(tc)
        % An unrecognized name must throw unknownLibrary immediately.
        % recipe() must not call any external process or file access.
            tc.verifyError(@() emk.setup.recipe("nosuchlibrary"), ...
                "emk:setup:recipe:unknownLibrary", ...
                "Unknown library name must throw unknownLibrary");
        end

        function test_unknownLibrary_errorMessage_containsName(tc)
        % Error message must echo the rejected name for debugging.
            badName = "nosuchlibrary";
            ME = tc.captureError(@() emk.setup.recipe(badName));
            tc.assertNotEmpty(ME, "Expected unknownLibrary error");
            tc.verifyClass(ME, "MException", "Error must be MException");
            tc.verifySubstring(ME.message, char(badName), ...
                "Error message must contain the rejected library name");
        end

        function test_unknownLibrary_errorMessage_listsSupportedNames(tc)
        % Error message must list supported names so the user can correct
        % the spelling.  At least 'pubchempy' must appear.
            ME = tc.captureError(@() emk.setup.recipe("bad_name"));
            tc.assertNotEmpty(ME, "Expected unknownLibrary error");
            hasPub = contains(ME.message, "pubchempy");
            hasMor = contains(ME.message, "mordred");
            hasBio = contains(ME.message, "biopython");
            tc.verifyTrue(hasPub || hasMor || hasBio, ...
                "Error message must list at least one supported Track 1 name");
        end

        % ------------------------------------------------------------------
        % TC2-TC4: Track 1 recipes (no exception)
        % ------------------------------------------------------------------

        function test_pubchempy_noThrow(tc)
        % recipe("pubchempy") must display a Track 1 recipe without throwing.
            ME = tc.captureError(@() emk.setup.recipe("pubchempy"));
            tc.verifyEmpty(ME, ...
                "recipe(""pubchempy"") must not throw any exception");
        end

        function test_mordred_noThrow(tc)
        % recipe("mordred") must display a Track 1 recipe without throwing.
            ME = tc.captureError(@() emk.setup.recipe("mordred"));
            tc.verifyEmpty(ME, ...
                "recipe(""mordred"") must not throw any exception");
        end

        function test_biopython_noThrow(tc)
        % recipe("biopython") must display a Track 1 recipe without throwing.
            ME = tc.captureError(@() emk.setup.recipe("biopython"));
            tc.verifyEmpty(ME, ...
                "recipe(""biopython"") must not throw any exception");
        end

        % ------------------------------------------------------------------
        % TC5-TC7: Track 2 recipes (no exception)
        % ------------------------------------------------------------------

        function test_openbabel_noThrow(tc)
        % recipe("openbabel") must display a Track 2 recipe (including GPL
        % warning) without throwing.
            ME = tc.captureError(@() emk.setup.recipe("openbabel"));
            tc.verifyEmpty(ME, ...
                "recipe(""openbabel"") must not throw any exception");
        end

        function test_mdanalysis_noThrow(tc)
        % recipe("mdanalysis") must display a Track 2 recipe without throwing.
            ME = tc.captureError(@() emk.setup.recipe("mdanalysis"));
            tc.verifyEmpty(ME, ...
                "recipe(""mdanalysis"") must not throw any exception");
        end

        function test_pymol_noThrow(tc)
        % recipe("pymol") must display a Track 2 recipe without throwing.
            ME = tc.captureError(@() emk.setup.recipe("pymol"));
            tc.verifyEmpty(ME, ...
                "recipe(""pymol"") must not throw any exception");
        end

        % ------------------------------------------------------------------
        % TC8: All canonical names are covered (regression guard)
        % ------------------------------------------------------------------

        function test_allCanonicalNames_listedInErrorMessage(tc)
        % The unknownLibrary error message must list every canonical name.
        % This is a regression guard: if a new library is added to recipe()
        % but its name is omitted from the error message, this test fails.
        % Individual recipe output is already tested in TC2-TC7, TC10-TC13,
        % so this test avoids calling recipe() again (Command Window output).
            ME = tc.captureError(@() emk.setup.recipe("__invalid__"));
            tc.assertNotEmpty(ME, "Expected unknownLibrary for unknown name");
            allNames = ["pubchempy", "mordred", "biopython", ...
                        "meeko", "vina", "pdbfixer", "docking", ...
                        "openbabel", "mdanalysis", "pymol"];
            for i = 1:numel(allNames)
                tc.verifyTrue(contains(ME.message, allNames(i)), ...
                    "Error message must list canonical name: " + allNames(i));
            end
        end

        % ------------------------------------------------------------------
        % TC9: char input coercion
        % ------------------------------------------------------------------

        function test_charInput_unknownLibrary_throwsUnknownLibrary(tc)
        % recipe() must accept char input (arguments block coerces char to
        % string).  Use an unknown name for a deterministic early failure.
            tc.verifyError(@() emk.setup.recipe('nosuchlibrary'), ...
                "emk:setup:recipe:unknownLibrary", ...
                "char input must be coerced to string and throw unknownLibrary");
        end

        % ------------------------------------------------------------------
        % TC10-TC13: A-3 docking recipes (no exception)
        % ------------------------------------------------------------------

        function test_meeko_noThrow(tc)
        % recipe("meeko") must display a Track 1 recipe without throwing.
            ME = tc.captureError(@() emk.setup.recipe("meeko"));
            tc.verifyEmpty(ME, ...
                "recipe(""meeko"") must not throw any exception");
        end

        function test_meeko_recipe_mentionsLGPL(tc)
        % meeko is licensed under LGPL-2.1.  The recipe must inform the user
        % about this license so they can make an informed decision.
        % Recipe output is captured by redirecting logInfo via env var.
        % Indirect approach: call recipe(), then verify it did not throw
        % AND that the meeko recipe case exists in the switch (non-throw is
        % already tested; this test verifies the LGPL mention is in the
        % recipe source, not the output, via a static code read).
        %
        % Implementation note: since logInfo writes to the Command Window
        % (not capturable without a custom test plugin), we verify the
        % recipe source string contains 'LGPL'.
            projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
            recipePath = fullfile(projectRoot, "src", "+emk", "+setup", "recipe.m");
            tc.assumeTrue(isfile(recipePath), "recipe.m not found");
            src = fileread(recipePath);
            % Find the meeko case block in the switch statement
            % and verify it mentions LGPL.
            idx = strfind(src, 'case "meeko"');
            tc.assertNotEmpty(idx, "recipe.m must contain a 'meeko' case");
            % Extract a 600-char window after the meeko case
            winEnd = min(idx(1) + 600, numel(src));
            meekoBlock = src(idx(1):winEnd);
            tc.verifyTrue(contains(meekoBlock, "LGPL"), ...
                "meeko recipe block must mention 'LGPL' to inform users of license");
        end

        function test_meeko_recipe_mentionsComplianceDoc(tc)
        % meeko recipe must reference docs/compliance.md so users can find
        % the full LGPL-2.1 usage rationale (CL-7).
            projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
            recipePath = fullfile(projectRoot, "src", "+emk", "+setup", "recipe.m");
            tc.assumeTrue(isfile(recipePath), "recipe.m not found");
            src = fileread(recipePath);
            idx = strfind(src, 'case "meeko"');
            tc.assertNotEmpty(idx, "recipe.m must contain a 'meeko' case");
            winEnd = min(idx(1) + 600, numel(src));
            meekoBlock = src(idx(1):winEnd);
            tc.verifyTrue(contains(meekoBlock, "compliance"), ...
                "meeko recipe must reference docs/compliance.md for full rationale");
        end

        function test_meeko_recipe_mentionsTrack1(tc)
        % meeko is a Track 1 library.  The recipe output must identify it as
        % Track 1 (not Track 2), so users know installExtra() is the right path.
        % Verified via recipe source: printTrack1_ is called (not printTrack2_).
            projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
            recipePath = fullfile(projectRoot, "src", "+emk", "+setup", "recipe.m");
            tc.assumeTrue(isfile(recipePath), "recipe.m not found");
            src = fileread(recipePath);
            idx = strfind(src, 'case "meeko"');
            tc.assertNotEmpty(idx, "recipe.m must contain a 'meeko' case");
            winEnd = min(idx(1) + 400, numel(src));
            meekoBlock = src(idx(1):winEnd);
            tc.verifyTrue(contains(meekoBlock, "printTrack1_"), ...
                "meeko recipe must call printTrack1_ (Track 1 classification)");
        end

        function test_vina_noThrow(tc)
        % recipe("vina") must display a Track 1 recipe without throwing.
            ME = tc.captureError(@() emk.setup.recipe("vina"));
            tc.verifyEmpty(ME, ...
                "recipe(""vina"") must not throw any exception");
        end

        function test_pdbfixer_noThrow(tc)
        % recipe("pdbfixer") must display a Track 1 recipe without throwing.
            ME = tc.captureError(@() emk.setup.recipe("pdbfixer"));
            tc.verifyEmpty(ME, ...
                "recipe(""pdbfixer"") must not throw any exception");
        end

        function test_pdbfixer_recipe_mentionsOpenmm(tc)
        % pdbfixer auto-installs openmm >= 8.2 as a dependency.
        % The recipe must warn about this download so users are not surprised
        % by a 70 MB install.  Verified via printPdbfixer_ function body in
        % recipe.m (the case block itself just calls printPdbfixer_() so the
        % 500-char case window would be too narrow -- search the function body).
            projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
            recipePath = fullfile(projectRoot, "src", "+emk", "+setup", "recipe.m");
            tc.assumeTrue(isfile(recipePath), "recipe.m not found");
            src = fileread(recipePath);
            % Search in the printPdbfixer_ function body (not the case block,
            % which only contains the single-line printPdbfixer_() call).
            fnIdx = strfind(src, "function printPdbfixer_");
            tc.assertNotEmpty(fnIdx, ...
                "recipe.m must define function printPdbfixer_()");
            fnEnd = min(fnIdx(1) + 2000, numel(src));
            pdbBlock = src(fnIdx(1):fnEnd);
            tc.verifyTrue(contains(pdbBlock, "openmm"), ...
                "printPdbfixer_() must mention 'openmm' dependency");
        end

        function test_pdbfixer_recipe_windowsMentionsSAC(tc)
        % On Windows, pdbfixer is blocked by Smart App Control (SAC).
        % The recipe source must mention SAC so Desktop users understand
        % the root cause and are not confused by a silent failure.
            projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
            recipePath = fullfile(projectRoot, "src", "+emk", "+setup", "recipe.m");
            tc.assumeTrue(isfile(recipePath), "recipe.m not found");
            src = fileread(recipePath);
            fnIdx = strfind(src, "function printPdbfixer_");
            tc.assertNotEmpty(fnIdx, "recipe.m must define function printPdbfixer_()");
            fnEnd = min(fnIdx(1) + 2000, numel(src));
            pdbBlock = src(fnIdx(1):fnEnd);
            tc.verifyTrue(contains(pdbBlock, "Smart App Control") || ...
                contains(pdbBlock, "SAC"), ...
                "printPdbfixer_() must mention Smart App Control (SAC) as root cause on Windows");
        end

        function test_pdbfixer_recipe_mentionsOnlineWorkaround(tc)
        % pdbfixer on Windows is blocked; MATLAB Online is the recommended
        % workaround.  The recipe source must mention MATLAB Online so users
        % know an alternative path exists.
            projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
            recipePath = fullfile(projectRoot, "src", "+emk", "+setup", "recipe.m");
            tc.assumeTrue(isfile(recipePath), "recipe.m not found");
            src = fileread(recipePath);
            fnIdx = strfind(src, "function printPdbfixer_");
            tc.assertNotEmpty(fnIdx, "recipe.m must define function printPdbfixer_()");
            fnEnd = min(fnIdx(1) + 2000, numel(src));
            pdbBlock = src(fnIdx(1):fnEnd);
            tc.verifyTrue(contains(lower(pdbBlock), "online"), ...
                "printPdbfixer_() must mention MATLAB Online as workaround path");
        end

        function test_pdbfixer_recipe_mentionsOptionalLibrariesFlag(tc)
        % The pdbfixer recipe must show the cfg.optionalLibraries.pdbfixer
        % flag so Online users know how to trigger the installation from
        % main_rdkit.m without reading docs.
            projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
            recipePath = fullfile(projectRoot, "src", "+emk", "+setup", "recipe.m");
            tc.assumeTrue(isfile(recipePath), "recipe.m not found");
            src = fileread(recipePath);
            fnIdx = strfind(src, "function printPdbfixer_");
            tc.assertNotEmpty(fnIdx, "recipe.m must define function printPdbfixer_()");
            fnEnd = min(fnIdx(1) + 2000, numel(src));
            pdbBlock = src(fnIdx(1):fnEnd);
            tc.verifyTrue(contains(pdbBlock, "optionalLibraries"), ...
                "printPdbfixer_() must mention cfg.optionalLibraries flag for Online install path");
        end

        function test_docking_noThrow(tc)
        % recipe("docking") must display the combined Track 1 pipeline recipe
        % (meeko + vina + pdbfixer) without throwing any exception.
            ME = tc.captureError(@() emk.setup.recipe("docking"));
            tc.verifyEmpty(ME, ...
                "recipe(""docking"") must not throw any exception");
        end

        function test_docking_recipe_mentionsAllThreeComponents(tc)
        % The docking combined recipe must name all three components
        % (meeko, vina, pdbfixer) so users know what will be installed.
        % Verified via recipe source: printDocking_() must reference all three.
            projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
            recipePath = fullfile(projectRoot, "src", "+emk", "+setup", "recipe.m");
            tc.assumeTrue(isfile(recipePath), "recipe.m not found");
            src = fileread(recipePath);
            idx = strfind(src, "printDocking_");
            % Find the first function body (not the call site)
            fnIdx = [];
            for k = 1:numel(idx)
                lineStart = max(1, idx(k) - 10);
                if contains(src(lineStart:idx(k)+15), "function")
                    fnIdx = idx(k);
                    break;
                end
            end
            tc.assumeTrue(~isempty(fnIdx), "printDocking_ function not found in recipe.m");
            fnEnd = min(fnIdx + 1000, numel(src));
            dockBlock = src(fnIdx:fnEnd);
            tc.verifyTrue(contains(dockBlock, "meeko"), ...
                "docking recipe (printDocking_) must name 'meeko'");
            tc.verifyTrue(contains(dockBlock, "vina"), ...
                "docking recipe (printDocking_) must name 'vina'");
            tc.verifyTrue(contains(dockBlock, "pdbfixer"), ...
                "docking recipe (printDocking_) must name 'pdbfixer'");
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
