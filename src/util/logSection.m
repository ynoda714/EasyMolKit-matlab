function logSection(scriptId, label, layer)
% logSection  Print a section-start banner to the MATLAB Command Window.
%
%   logSection(scriptId, label, layer)
%
%   Emits an INFO-level banner that identifies the running script, the
%   current section title, and its tutorial layer.  Call this as the first
%   executable line of every %% section in tutorial scripts so that each
%   section leaves a clear trace in the log.
%
%   Arguments:
%     scriptId  (1,1) string  - Script identifier (e.g. "F01", "S01", "R01")
%     label     (1,1) string  - Section header text (e.g. "Section 0: Setup")
%     layer     (1,1) string  - Tutorial layer (e.g. "Foundation L1")
%
%   Output format:
%     [HH:MM:SS][INFO]  --- F01 | Section 0: Setup  [Foundation L1] ---
%
%   Example:
%     logSection("S01", "Section 3: Compute Fingerprints", "Stories L2")
%     % [11:36:40][INFO]  --- S01 | Section 3: Compute Fingerprints  [Stories L2] ---
%
%   See also: logInfo, logWarn, logError, logProgress

    arguments
        scriptId  (1,1) string
        label     (1,1) string
        layer     (1,1) string
    end

    logInfo("--- %s | %s  [%s] ---", scriptId, label, layer);
end
