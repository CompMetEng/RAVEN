classdef tInstallation < RavenTestCase
% tInstallation  Tests for the setup/installation helpers in installation/.
%
%   addRavenToUserPath, removeRavenFromPath and updateDocumentation have
%   destructive / irreversible side effects (rewrite the user's startup.m,
%   remove RAVEN from the live path, regenerate the whole doc tree). They are
%   therefore only checked for existence here, not executed.

    methods (Test)

        function findRAVENrootReturnsInstallFolder(testCase)
            p = findRAVENroot();
            testCase.verifyTrue(isfolder(p));
            testCase.verifyTrue(exist(fullfile(p,'src','installation','findRAVENroot.m'),'file')==2);
        end

        function findRAVENrootSurvivesClearAndRepeatCalls(testCase)
            % The resolved path is cached for the session (findRAVENroot is
            % called on every solver invocation, so re-reading the
            % preference from disk every time is wasteful). Repeated calls
            % must keep returning the same root, and it must still resolve
            % correctly after being cleared, which is what checkRaven
            % does when the RAVEN.ravenPath preference changes.
            p1 = findRAVENroot();
            p2 = findRAVENroot();
            testCase.verifyEqual(p2, p1);
            clear('findRAVENroot');
            p3 = findRAVENroot();
            testCase.verifyEqual(p3, p1);
        end

        function checkRavenReturnsVersion(testCase)
            [~, currVer] = evalc('checkRaven(false, false)');
            testCase.verifyNotEmpty(currVer);
        end

        function checkInstallationWarnsAndForwards(testCase)
            % checkRaven itself issues other warnings during its checks (e.g.
            % SBML id-prefixing), which would overwrite lastwarn before it
            % could be inspected here. verifyWarning tracks every warning
            % issued during the call, not just the most recent one.
            testCase.verifyWarning(@() evalc('checkInstallation(false, false)'), 'RAVEN:deprecated');
            [~, currVer] = evalc('checkInstallation(false, false)');
            testCase.verifyNotEmpty(currVer);
        end

        function checkFunctionUniquenessRuns(testCase)
            [~, status] = evalc('checkFunctionUniqueness()');
            testCase.verifyTrue(islogical(status) || isnumeric(status));
        end

        function fetchToolboxAssetChecksDownloads(testCase)
            % A small published asset is downloaded and matches the
            % checksum of the toolbox-assets manifest; a second call reuses
            % the file.
            d = tempname;
            cleanup = onCleanup(@() rmdir(d, 's'));
            p = fetchToolboxAsset(d, 'transporters-2026-07-01', 'tcdb_substrates.tsv');
            testCase.verifyTrue(isfile(p));
            testCase.verifyFalse(isfile([p '.part']));
            info = dir(p);
            testCase.verifyEqual(info.bytes, 219635);
            testCase.verifyEqual(fetchToolboxAsset(d, 'transporters-2026-07-01', 'tcdb_substrates.tsv'), p);
        end

        function addRavenToUserPathExists(testCase)
            testCase.verifyEqual(exist('addRavenToUserPath','file'), 2);
        end

        function removeRavenFromPathExists(testCase)
            testCase.verifyEqual(exist('removeRavenFromPath','file'), 2);
        end

        function updateDocumentationExists(testCase)
            testCase.verifyEqual(exist('updateDocumentation','file'), 2);
        end

    end
end
