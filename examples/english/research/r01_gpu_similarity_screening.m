%[text] # R01: Large-scale Similarity Screening Using GPU Acceleration
%[text] EasyMolKit Research -- Layer 4
%[text] ## Story
%[text] Suppose the pharmaceutical informatics team has constructed a virtual compound library of 100,000 entries from a commercial synthesis catalog.
%[text] Before running all compounds through costly experimental assays, quickly list candidates structurally similar to the hit compound "Ibuprofen," which has confirmed COX-2 inhibitory activity.
%[text] Calculating similarity one molecule at a time takes several minutes for 100,000 screenings, but arranging fingerprints as a matrix allows for similarity calculations in a single matrix operation.
%[text] GPUs excel at this matrix operation, completing screenings of a million entries in milliseconds.
%[text] ## Learning Objectives
%[text] 1. Understand why naive loops become a bottleneck as the number of compounds increases.
%[text] 2. Represent fingerprints as single-precision binary matrices.
%[text] 3. Derive the Tanimoto coefficient as a matrix-vector product.
%[text] 4. Implement GPU screening using MATLAB's `gpuArray`.
%[text] 5. Benchmark CPU and GPU throughput for multiple library sizes (1k / 10k / 100k).
%[text] 6. Select the top N hits from the FDA drug library.
%[text] 7. Understand GPU memory constraints and batch strategies. \
%[text] ## Prerequisites
%[text] - Completion of S04 (Virtual Screening) -- Understanding of Tanimoto and rankBy
%[text] - Completion of A02 (Molecular Clustering) -- Understanding of fingerprint matrices
%[text] - Parallel Computing Toolbox (gpuArray, gpuDevice)
%[text] - CUDA-enabled GPU (NVIDIA recommended) \
%[text]   If no GPU is detected, Sections 5-7 will execute the same matrix operations on the CPU. Learning the code and concepts remains possible.
%[text] ## Operating Environment
%[text] GPU acceleration (gpuArray/gpuDevice) requires a CUDA-enabled GPU on the machine running MATLAB. MATLAB Online operates on cloud servers and cannot connect to local GPUs.
%[text:table]{"columnWidths":[-1,258],"ignoreHeader":true}
%[text] | Environment | Operation |
%[text] | --- | --- |
%[text] |  MATLAB Desktop (with CUDA GPU) | Full GPU Acceleration  |
%[text] |  MATLAB Desktop (without CUDA GPU) |  CPU Matrix Fallback (Automatic) |
%[text] | MATLAB Online |  CPU Matrix Fallback (Automatic) |
%[text:table]
%[text] Sections 1-4, 7-9 can be executed unchanged in any environment.
%[text] Sections 5-6 will fallback to CPU matrix operations if no GPU is available.
%[text] ## Data
%[text] - `data/list/fda_drugs.csv` -- 200 FDA-approved drugs (ChEMBL, CC-BY-SA 3.0) \
%[text] ## References
%[text] - Cao Y et al. (2008) ChemFP: fast fingerprint searching. J Chem Inf Model 48:2208-2215. doi:10.1021/ci8001854
%[text] - Willett P (2013) Combination of similarity rankings using data fusion. J Chem Inf Model 53:1-10. doi:10.1021/ci300547g
%[text] - Rogers D & Hahn M (2010) Extended-connectivity fingerprints. J Chem Inf Model 50:742-754. doi:10.1021/ci100050t
%[text] - MATLAB Parallel Computing Toolbox: https://www.mathworks.com/help/parallel-computing/gpuarray.html \
%[text] Execution: Run each section with Ctrl+Enter. Duration: 45-90 minutes
%%
%[text] ## Section 0: Setup
% Resolve project root (works for Desktop, MCP, and MATLAB Online)
sDir = fileparts(mfilename('fullpath'));
if strlength(sDir) > 0
    addpath(genpath(fullfile(sDir, '..', '..', '..', 'src')));
elseif isfolder(fullfile(pwd, 'src'))
    addpath(genpath(fullfile(pwd, 'src')));
elseif ~isempty(which("logInfo"))
    addpath(genpath(fileparts(fileparts(which("logInfo")))));
end
projectRoot = resolveProjectRoot();
addpath(genpath(fullfile(projectRoot, 'src')));
emk.setup.initPython(); %[output:840b68c2]
%[text] Warm-up Python/RDKit process (initialization before execution)
mol_warmup = emk.mol.fromSmiles("C");   % Methane -- lightweight molecule
clear mol_warmup;
%[text] Detect the presence of a GPU.
hasGPU = false;
try %[output:group:59f4d6c2]
    d      = gpuDevice();
    hasGPU = true;
    logInfo("R01: GPU detected -- %s (%.1f GB VRAM)", ... %[output:47d13b61]
        d.Name, d.TotalMemory / 1e9); %[output:47d13b61]
catch
    logWarn("R01: CUDA GPU not detected.");
    logWarn("     GPU section will fall back to CPU matrix operations.");
    logWarn("     A CUDA-enabled GPU is required to experience full speed-up.");
end %[output:group:59f4d6c2]
logInfo("R01: Setup complete"); %[output:56b8d46d]
%%
%[text] ## Section 1: Hit Compounds and Scaling Issues
%[text] Similarity screening against a large library can be approached at three tiers, each with a very different throughput profile.
%[text] **Tier 1 -- Molecule-by-molecule comparison (****`emk.similarity.tanimoto`****)**
%[text] Every call crosses the MATLAB-Python boundary, so inter-process communication dominates. Expect roughly 500 molecules/second -- screening 100,000 takes about 200 s, and 1 million would take 33 minutes.
%[text] **Tier 2 -- RDKit Bulk API (****`emk.similarity.rankBy`****)**
%[text] All fingerprints are scored in a single call, crossing the MATLAB-Python boundary only once. Throughput rises to roughly 200,000 molecules/second.
%[text] **Tier 3 -- GPU Matrix Multiplication (this exercise)**
%[text] Fingerprints are packed into a binary matrix and all Tanimoto scores are computed simultaneously via a single matrix-vector product. With a CUDA GPU, throughput can exceed 5 million molecules/second -- though the effective rate depends on CPU-to-GPU transfer time.
%[text] Section 6 benchmarks both end-to-end time (including transfer) and pre-loaded time (matrix already on GPU), so you can see exactly where the time goes.
%[text] 
%[text] GPU acceleration pays off most when:
%[text] (a) the library exceeds ~100,000 compounds, or
%[text] (b) the same library is queried repeatedly (the matrix stays resident on the GPU).
QUERY_SMILES = "CC(C)Cc1ccc(cc1)C(C)C(=O)O";   % Ibuprofen
QUERY_NAME   = "Ibuprofen (COX-2 inhibitor)";
mol_query = emk.mol.fromSmiles(QUERY_SMILES);
fp_query  = emk.fingerprint.morgan(mol_query);
logInfo("Query: %s", QUERY_NAME); %[output:80238221]
logInfo("Canonical SMILES: %s", emk.mol.toSmiles(mol_query)); %[output:736aa985]
logInfo("FP size: %d bits, on-bits: %d", ... %[output:group:5457eb7f] %[output:5f53679a]
    double(fp_query.GetNumBits()), double(fp_query.GetNumOnBits())); %[output:group:5457eb7f] %[output:5f53679a]
%%
%[text] ## Section 2: Actual Library Construction (FDA Drugs)
%[text] We now load SMILES data from the FDA drug CSV to build a screening library to compare against the query from Section 1. Entries that RDKit cannot parse are automatically skipped; only valid molecules are kept.
%[text] We then compute ECFP4 fingerprints (Extended Connectivity FingerPrint, radius 2, 2048 bits) for all molecules. ECFP4 encodes the local chemical environment around each atom into a 2048-bit string; converting this to a matrix enables the fast calculations in subsequent sections.
DATA_FILE = "data/list/fda_drugs.csv";
rawTbl    = readtable(DATA_FILE, TextType="string");
logInfo("Loaded %d entries from %s", height(rawTbl), DATA_FILE); %[output:9465a5aa]
nLib  = height(rawTbl);
mols  = cell(1, nLib);
valid = false(1, nLib);
for k = 1:nLib
    try
        mols{k} = emk.mol.fromSmiles(rawTbl.SMILES(k));
        valid(k) = true;
    catch
        % Skip entries that cannot be parsed
    end
end
validIdx = find(valid);
mols     = mols(validIdx);
names    = rawTbl.Name(validIdx);
smiles   = rawTbl.SMILES(validIdx);
nValid   = numel(mols);
logInfo("Parsed %d / %d molecules", nValid, nLib); %[output:5dac7e8f]
%[text] Calculate the ECFP4 fingerprints for the actual library.
fps = cell(1, nValid);
for k = 1:nValid
    fps{k} = emk.fingerprint.morgan(mols{k});
end
logInfo("Completed calculation of %d ECFP4 fingerprints (Radius=2, NBits=2048)", nValid); %[output:64d2b1f1]
%%
%[text] ## Section 3: Convert Fingerprints to Bit Matrix
%[text] Recall that each fingerprint is simply a string of 0s and 1s. To score $N$ molecules simultaneously, we stack these strings as rows of a matrix $F$ of size $N \\times B$, where $B$ is the bit length:
%[text] $F(i,j) = 1$ if bit $j$ of molecule $i$ is ON, and 0 otherwise.
%[text] ### Memory Footprint (Single-Precision float32)
%[text:table]{"ignoreHeader":true}
%[text] | Number of Molecules | Memory |
%[text] | --- | --- |
%[text] | 10,000 molecules x 2048 bits | 82 MB |
%[text] | 100,000 molecules x 2048 bits | 819 MB |
%[text] | 1,000,000 molecules x 2048 bits | 8\.2 GB (batch processing required) |
%[text:table]
%[text] We store the matrix as `single` (float32, 4 bytes/element) rather than `logical` (1 byte) because GPU matrix-multiplication kernels are optimized for float32. Converting upfront avoids an implicit type cast at computation time.
NBITS = double(fps{1}.GetNumBits());   % 2048
tic
fpMat = false(nValid, NBITS);          % N x B logical (on host)
for k = 1:nValid
    fpMat(k, :) = emk.fingerprint.toArray(fps{k});
end
t_convert = toc;
F_host  = single(fpMat);              % N x B single (ready for GPU transfer)
q_vec   = single(emk.fingerprint.toArray(fp_query));   % 1 x B single
logInfo("FP Matrix: %d x %d  (%.1f MB as single, constructed in %.2f s)", ... %[output:group:71794e8c] %[output:19813262]
    nValid, NBITS, nValid * NBITS * 4 / 1e6, t_convert); %[output:group:71794e8c] %[output:19813262]
%%
%[text] ## Section 4: Tier 2 Baseline -- RDKit Bulk API
%[text] `emk.similarity.rankBy()` calls RDKit's BulkTanimotoSimilarity. It scores all fingerprints in the library in a single call, requiring only one data exchange between MATLAB and Python.
%[text] This result will serve as the accuracy benchmark (ground truth) for the subsequent GPU implementation.
tic
result_bulk = emk.similarity.rankBy(fp_query, fps);
t_bulk      = toc;
%[text] Rearrange all scores in the original index order (indices may differ from the original order).
T_rdkit = zeros(nValid, 1);
T_rdkit(result_bulk.Indices) = result_bulk.Scores;
logInfo("Tier 2 (RDKit Bulk): Scored %d entries in %.4f s  (%.1f k mol/s)", ... %[output:group:06361ea7] %[output:28f71f40]
    nValid, t_bulk, nValid / t_bulk / 1e3); %[output:group:06361ea7] %[output:28f71f40]
%%
%[text] ## Section 5: Tier 3 -- GPU Tanimoto via Matrix Multiplication
%[text] For binary (0/1) fingerprints, the Tanimoto coefficient between molecule $i$ and query $q$ is:
%[text]{"align":"center"} $T(i) = \\frac{|A(i) \\cap q|}{|A(i) \\cup q|} = \\frac{A(i) \\cdot q}{|A(i)| + |q| - A(i) \\cdot q}$
%[text] where $A(i)$ denotes the bit string of molecule $i$, $q$ the query, and $A(i) \\cdot q$ the count of bits that are ON in both (the inner product).
%[text] Notice that stacking all $N$ fingerprints into a matrix $F$ ($N \\times B$) lets us compute all $N$ inner products in one shot as the matrix-vector product $F q^\\top$ ($N \\times 1$). In code:
%[text]     intersection = F \* q'         % ON bits shared with the query, for all molecules
%[text]     onBitsA      = sum(F, 2)      % ON-bit count per molecule
%[text]     T            = intersection ./ (onBitsA + sum(q) - intersection)
%[text] The GPU executes this product across thousands of parallel cores, so the advantage over the CPU grows with $N$.
tic
if hasGPU
    F_gpu      = gpuArray(F_host);        % Transfer N x B to GPU
    q_gpu      = gpuArray(q_vec);         % 1 x B (on GPU)
    inter_gpu  = F_gpu * q_gpu';          % N x 1  (parallel GEMV)
    onBitsA    = sum(F_gpu, 2);           % N x 1
    onBitsQ    = sum(q_gpu);              % Scalar
    union_gpu  = onBitsA + onBitsQ - inter_gpu;
    T_gpu      = gather(inter_gpu ./ union_gpu);   % Retrieve to host
    T_gpu      = double(T_gpu);
else
    % CPU fallback: same formula without gpuArray
    inter_cpu  = F_host * q_vec';
    onBitsA    = sum(F_host, 2);
    onBitsQ    = sum(q_vec);
    union_cpu  = onBitsA + onBitsQ - inter_cpu;
    T_gpu      = double(inter_cpu ./ union_cpu);
    logWarn("No GPU -- Calculated Tanimoto with CPU single-precision matrix.");
end
t_gpu = toc;
logInfo("Tier 3 (GPU matrix): Scored %d entries in %.4f s  (%.0f k mol/s)", ... %[output:group:1a22afc2] %[output:86732921]
    nValid, t_gpu, nValid / t_gpu / 1e3); %[output:group:1a22afc2] %[output:86732921]
%[text] Compare with RDKit Bulk scores.
maxDiff = max(abs(T_gpu - T_rdkit));
logInfo("Max |GPU - RDKit| difference: %.2e  (expected < 1e-5)", maxDiff); %[output:4e351304]
if maxDiff > 1e-4 %[output:group:35f181e2]
    logWarn("Scores are misaligned -- Check FP conversion in Section 3.");
else
    logInfo("Validation PASSED -- GPU matrix and RDKit Bulk scores match."); %[output:460419f8]
end %[output:group:35f181e2]
%%
%[text] ## Section 6: Scale-up Benchmark -- Synthetic Library
%[text] The real FDA library (~200 molecules) is far too small to reveal the GPU's advantage -- when $N$ is small, the CPU-to-GPU transfer cost alone can exceed computation time.
%[text] To benchmark at realistic scales, we synthesize a larger library from the FDA fingerprints:
%[text] 1. Tile the FDA rows cyclically to reach $N$ entries (preserving ON-bit density).
%[text] 2. Randomly flip 5% of the bits to introduce chemical variation. \
%[text] These synthetic fingerprints are not chemically meaningful, but their bit density matches the real library, making them representative for timing purposes.
%[text] We compare three approaches across library sizes of 1k, 10k, and 100k:
%[text] - **CPU Matrix**: Single-precision matrix-vector product on the CPU
%[text] - **GPU end-to-end**: Total wall time including CPU-\>GPU transfer, computation, and GPU-\>CPU transfer
%[text] - **GPU pre-loaded**: Computation time only, with the matrix already resident on the GPU -- the realistic scenario when the same library is queried many times \
rng(42);   % Reproducible synthetic library
SCALE_SIZES  = [1e3, 1e4, 1e5];
t_cpu_mat    = zeros(1, numel(SCALE_SIZES));
t_gpu_mat    = zeros(1, numel(SCALE_SIZES));
t_gpu_pre    = zeros(1, numel(SCALE_SIZES));   % GPU pre-loaded (compute only)
logInfo("--- Scale-up Benchmark ---"); %[output:0ec05bfb]
for si = 1:numel(SCALE_SIZES) %[output:group:5e13348d]
    N_syn = SCALE_SIZES(si);
    % Construct synthetic FP matrix
    baseRows  = mod((0 : N_syn - 1), nValid) + 1;   % Cycle execution
    F_syn     = F_host(baseRows, :);                  % N_syn x B (single)
    flipMask  = single(rand(N_syn, NBITS, "single") < 0.05);
    F_syn     = single(abs(F_syn - flipMask) >= 0.5); % XOR using abs trick
    % CPU matrix calculation
    tic
    ig_cpu   = F_syn * q_vec';
    oa_cpu   = sum(F_syn, 2);
    oq_cpu   = sum(q_vec);
    T_cpu_s  = double(ig_cpu ./ (oa_cpu + oq_cpu - ig_cpu)); %#ok<NASGU>
    t_cpu_mat(si) = toc;
    % GPU end-to-end (including transfer + computation + gather)
    tic
    if hasGPU
        Fg   = gpuArray(F_syn);
        qg   = gpuArray(q_vec);
        ig   = Fg * qg';
        oa   = sum(Fg, 2);
        oq   = sum(qg);
        T_g  = gather(ig ./ (oa + oq - ig)); %#ok<NASGU>
    else
        % Same as CPU branch -- reported as "GPU (fallback)"
        ig   = F_syn * q_vec';
        oa   = sum(F_syn, 2);
        T_g  = double(ig ./ (oa + sum(q_vec) - ig)); %#ok<NASGU>
    end
    t_gpu_mat(si) = toc;
    % GPU pre-loaded (assuming matrix is already resident on GPU -- only timing computation + gather)
    if hasGPU
        Fg_pre = gpuArray(F_syn);
        qg_pre = gpuArray(q_vec);
        wait(gpuDevice());       % Ensure H->D transfer is complete before starting timing
        tic
        ig_p   = Fg_pre * qg_pre';
        oa_p   = sum(Fg_pre, 2);
        T_gp   = gather(ig_p ./ (oa_p + sum(qg_pre) - ig_p)); %#ok<NASGU>
        t_gpu_pre(si) = toc;
        clear Fg_pre qg_pre;
    else
        t_gpu_pre(si) = t_gpu_mat(si);   % Same as e2e in fallback
    end
    speedup_e2e = t_cpu_mat(si) / t_gpu_mat(si);
    speedup_pre = t_cpu_mat(si) / t_gpu_pre(si);
    logInfo("N=%7.0f | CPU %.4f s | GPU e2e %.4f s (%.1fx) | GPU preloaded %.4f s (%.1fx)", ... %[output:086b41df]
        N_syn, t_cpu_mat(si), t_gpu_mat(si), speedup_e2e, t_gpu_pre(si), speedup_pre); %[output:086b41df]
end %[output:group:5e13348d]
%%
%[text] ## Section 7: Benchmark Summary Plot
%[text] We now display the Section 6 results as a table and two log-scale plots: one comparing raw elapsed time and one showing the GPU speedup factor relative to the CPU baseline.
tblBench = table(SCALE_SIZES(:), t_cpu_mat(:), t_gpu_mat(:), t_gpu_pre(:), ...
    VariableNames=["N", "CPU_s", "GPU_e2e_s", "GPU_pre_s"]);
tblBench.CPU_kMolPerSec     = round(tblBench.N ./ tblBench.CPU_s      / 1e3, 1);
tblBench.GPU_e2e_kMolPerSec = round(tblBench.N ./ tblBench.GPU_e2e_s  / 1e3, 1);
tblBench.GPU_pre_kMolPerSec = round(tblBench.N ./ tblBench.GPU_pre_s  / 1e3, 1);
tblBench.Speedup_e2e        = round(tblBench.CPU_s ./ tblBench.GPU_e2e_s, 2);
tblBench.Speedup_pre        = round(tblBench.CPU_s ./ tblBench.GPU_pre_s, 2);
logInfo("--- Benchmark Summary ---"); %[output:5897ae42]
disp(tblBench); %[output:66892007]
figure("Name", "R01: CPU vs GPU Tanimoto Throughput"); %[output:42f6c091]
subplot(1, 2, 1); %[output:42f6c091]
loglog(SCALE_SIZES, t_cpu_mat, "b-o", LineWidth=1.5, DisplayName="CPU matrix"); %[output:42f6c091]
hold on; %[output:42f6c091]
loglog(SCALE_SIZES, t_gpu_mat, "r-s", LineWidth=1.5, DisplayName="GPU e2e (incl. transfer)"); %[output:42f6c091]
loglog(SCALE_SIZES, t_gpu_pre, "g-d", LineWidth=1.5, DisplayName="GPU pre-loaded (compute only)"); %[output:42f6c091]
xlabel("Library Size (Number of Molecules)"); %[output:42f6c091]
ylabel("Elapsed Time (s)"); %[output:42f6c091]
title("Screening Time vs Library Size"); %[output:42f6c091]
legend(Location="northwest"); %[output:42f6c091]
grid on; %[output:42f6c091]
subplot(1, 2, 2); %[output:42f6c091]
semilogx(SCALE_SIZES, tblBench.Speedup_e2e, "r-s", LineWidth=2, DisplayName="GPU e2e"); %[output:42f6c091]
hold on; %[output:42f6c091]
semilogx(SCALE_SIZES, tblBench.Speedup_pre, "g-d", LineWidth=2, DisplayName="GPU pre-loaded"); %[output:42f6c091]
yline(1, "k--", LineWidth=1.5, HandleVisibility="off"); %[output:42f6c091]
xlabel("Library Size (Number of Molecules)"); %[output:42f6c091]
ylabel("Speedup (CPU / GPU)"); %[output:42f6c091]
if hasGPU
    title("GPU Speedup Factor"); %[output:42f6c091]
else
    title("Speedup (CPU Fallback)");
end
legend(Location="northwest"); %[output:42f6c091]
grid on; %[output:42f6c091]
if hasGPU %[output:group:34699ed3]
    sgtitle("R01: CPU vs GPU Tanimoto Throughput"); %[output:42f6c091]
else
    sgtitle("R01: Benchmark (CPU fallback -- No GPU)");
end %[output:group:34699ed3]
%%
%[text] ## Section 8: Select Top N Hits from Actual FDA Library
%[text] The appropriate Tanimoto cutoff depends on the goal of the search. Willett (2013) summarizes the commonly accepted ranges:
%[text:table]{"ignoreHeader":true}
%[text] | Threshold | Interpretation |
%[text] | --- | --- |
%[text] | $T \> 0.85$ | Highly similar -- almost certainly the same scaffold |
%[text] | $T \> 0.60$ | Similar -- likely the same chemical class |
%[text] | $T \> 0.40$ | Moderately similar -- may share a pharmacophore |
%[text] | $T \> 0.30$ | Weakly similar -- hits occur but often non-specific |
%[text:table]
%[text] 
%[text] For patent novelty searches, analysts typically require $T \> 0.85$ to establish structural similarity. For lead hopping -- finding new actives by jumping across chemical space -- a lower cutoff of $T \> 0.30$ is common, accepting more non-specific hits in exchange for broader coverage.
TOP_N     = 10;
THRESHOLD = 0.30;
[T_sorted, sortIdx] = sort(T_rdkit, "descend");
hitMask   = T_sorted >= THRESHOLD;
hitIdx    = sortIdx(hitMask);
hitScores = T_sorted(hitMask);
if isempty(hitIdx) %[output:group:34b064d5]
    logWarn("No hits with a threshold of %.2f or higher -- consider lowering the THRESHOLD.", THRESHOLD);
else
    nHits  = min(TOP_N, numel(hitIdx));
    hitTbl = table( ...
        (1:nHits)', ...
        names(hitIdx(1:nHits)), ...
        smiles(hitIdx(1:nHits)), ...
        round(hitScores(1:nHits), 4), ...
        VariableNames=["Rank", "Name", "SMILES", "Tanimoto"]);
    logInfo("Top-%d hits (Tanimoto >= %.2f):", nHits, THRESHOLD); %[output:40ff28ed]
    disp(hitTbl); %[output:91bdcc14]
end %[output:group:34b064d5]
%[text] Notice that amphetamines appear in the hit list ($T \\approx 0.33$--$0.35$). They share the phenylethylamine scaffold (benzene ring + short alkyl chain + polar group) with ibuprofen, which ECFP4 picks up as a structural similarity.
%[text] ECFP4 encodes local chemical environments, not pharmacological activity. At low thresholds, structural coincidences like this are expected -- they are a well-known limitation of fingerprint-based similarity.
%[text] In practice, always combine a similarity search with activity-based filters (Lipinski Ro5, PAINS alerts, known target-class annotations) before selecting leads.
figure("Name", "R01: Tanimoto Score Distribution"); %[output:14cbbcd5]
histogram(T_rdkit, 40, FaceColor=[0.3 0.6 0.9], EdgeColor="w"); %[output:14cbbcd5]
xline(THRESHOLD, "r--", LineWidth=2, ... %[output:14cbbcd5]
    Label=sprintf("Threshold = %.2f", THRESHOLD), LabelVerticalAlignment="bottom"); %[output:14cbbcd5]
xlabel("Tanimoto Score (ECFP4 vs Ibuprofen)"); %[output:14cbbcd5]
ylabel("Count"); %[output:14cbbcd5]
title(sprintf("Tanimoto Distribution: %d FDA Drugs vs Ibuprofen", nValid)); %[output:14cbbcd5]
grid on; %[output:14cbbcd5]
%%
%[text] ## Section 9: Batch Strategy for Memory-Constrained GPUs
%[text] The fingerprint matrix requires $M = N \\times B \\times 4$ bytes in float32, where $N$ is the number of molecules and $B$ the bit length (2048 by default).
%[text] At 1 million compounds this is 8.2 GB; at 10 million it rises to 82 GB -- well beyond most consumer GPUs.
%[text] The standard solution is to stream the library through the GPU in batches, transferring one chunk at a time:
%[text]     BATCH = floor(gpu\_free  *0.8 / (B*  4));   % 80% safety margin
%[text]     T\_all = zeros(N\_total, 1);
%[text]     for start = 1 : BATCH : N\_total
%[text]         batchIdx = start : min(start + BATCH - 1, N\_total);
%[text]         Fb       = gpuArray(F\_big(batchIdx, :));
%[text]         inter    = Fb \* q\_gpu';
%[text]         onA      = sum(Fb, 2);
%[text]         T\_all(batchIdx) = gather(inter ./ (onA + onBitsQ - inter));
%[text]     end
%[text] 
%[text] Here, `gpu_free` is read from `gpuDevice().AvailableMemory`.
%[text] The 80% safety margin reserves headroom for MATLAB's internal GPU buffers and any other concurrent GPU workloads.
if hasGPU %[output:group:0f580fba]
    d        = gpuDevice();
    freeGB   = d.AvailableMemory / 1e9;
    BATCH    = floor(d.AvailableMemory * 0.8 / (NBITS * 4));
    logInfo("GPU free memory: %.1f GB", freeGB); %[output:0f35e493]
    logInfo("Safe batch size: %d molecules (%.1f MB/batch)", ... %[output:5035f545]
        BATCH, BATCH * NBITS * 4 / 1e6); %[output:5035f545]
    logInfo("Number of batches needed for screening 1 million: %.0f", ceil(1e6 / BATCH)); %[output:5f3193cf]
    logInfo("Number of batches needed for screening 10 million: %.0f", ceil(10e6 / BATCH)); %[output:56d9677d]
else
    logInfo("GPU batch size calculation: N/A (No GPU)");
    logInfo("For a GPU with 8 GB VRAM: Safe batch size ~ 980,000 molecules/pass");
    logInfo("Screening 10 million: Requires about 11 batches");
end %[output:group:0f580fba]
logInfo("R01: Complete -- GPU-accelerated similarity screening finished"); %[output:76ce6c50]
%%
%[text] ## Section 10: Summary -- GPU Screening Workflow
%[text] In this exercise, we built a pipeline combining fingerprint matrixing with GPU matrix multiplication to screen 100,000 compounds in milliseconds.
%[text] The table below maps each step to the algorithm and its MATLAB implementation.
%[text:table]{"ignoreHeader":true}
%[text] | Step | Method | MATLAB Implementation |
%[text] | --- | --- | --- |
%[text] | 1\. Encoding | SMILES -\> ECFP4 (2048 bits) | `emk.mol.fromSmiles()` + `emk.fingerprint.morgan()` |
%[text] | 2\. Matrixing | $M \\times 2048$, single precision float32 | `toArray()` + `cat()` / `reshape()` |
%[text] | 3\. GPU Tanimoto | Batch calculation via matrix-vector multiplication | `gpuArray(F) * q_gpu'` + element-wise division |
%[text] | 4\. Benchmark | Tier 2 vs with GPU transfer vs without GPU transfer | `tic`/`toc` x library size 1k/10k/100k |
%[text] | 5\. Hit Selection | Tanimoto score descending + threshold filter | `sort()` + `gather()` |
%[text] | 6\. Batch Strategy | Free VRAM x 0.8 / (B x 4) | `gpuDevice().AvailableMemory` |
%[text:table]
%[text] 
%[text] **Key Takeaways**
%[text] 1. Tier 1 (Python IPC per molecule) yields only ~500 mol/s -- 200 s for 100k compounds. Matrixing fingerprints steps performance up to Tier 2 (Bulk API, ~200k mol/s) and then Tier 3 (GPU, 5M+ mol/s).
%[text] 2. Storing fingerprints as `single` ($M \\times 2048$) uses half the memory of `double` and maps directly to float32 GPU kernels. Be aware that 1 million molecules already requires 8.2 GB -- near the limit of a typical consumer GPU.
%[text] 3. The Tanimoto numerator $F\_i \\cdot q$ is the dot product of row $i$ with the query. Stacking all rows into $F$ lets you compute all $N$ dot products at once as $F q^\\top$; combined with `sum(F,2)` for ON-bit counts, the full Tanimoto formula becomes GPU-parallel.
%[text] 4. Transfer to VRAM with `gpuArray()` and retrieve with `gather()`. Keeping $F$ resident (pre-loaded) eliminates transfer overhead and maximizes throughput when the same library is queried repeatedly.
%[text] 5. ECFP4 encodes topology, not bioactivity. At low thresholds ($T \> 0.30$), structural coincidences produce false positives -- as the amphetamine hits in Section 8 illustrate. Always pair similarity filtering with activity-aware filters (Lipinski Ro5, PAINS alerts).
%[text] 6. Safe batch size $= \\text{Free VRAM} \\times 0.8 / (B \\times 4)$. An 8 GB GPU handles roughly 980k molecules per pass, completing 10 million in about 11 batches. \
%[text] **Next Steps**
%[text] - R02: Pharmacokinetics Simulation -- Compare PK profiles of the hit candidates identified here
%[text] - Scale to a ChEMBL subset (1M+ compounds) to experience the batch strategy in practice
%[text] - Assess hit diversity with GPU k-means clustering (builds on A02)
%[text] - Explore ECFP radius and bit-count sensitivity (radius 3, or 4096 bits) and compare throughput \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:840b68c2]
%   data: {"dataType":"text","outputData":{"text":"[23:22:18][INFO]  initPython: Python 3.10 already active (Status: Loaded) -- skipping.\n","truncated":false}}
%---
%[output:47d13b61]
%   data: {"dataType":"text","outputData":{"text":"[23:22:18][INFO]  R01: GPU detected -- NVIDIA GeForce RTX 3060 Laptop GPU (6.4 GB VRAM)\n","truncated":false}}
%---
%[output:56b8d46d]
%   data: {"dataType":"text","outputData":{"text":"[23:22:18][INFO]  R01: Setup complete\n","truncated":false}}
%---
%[output:80238221]
%   data: {"dataType":"text","outputData":{"text":"[23:22:18][INFO]  Query: Ibuprofen (COX-2 inhibitor)\n","truncated":false}}
%---
%[output:736aa985]
%   data: {"dataType":"text","outputData":{"text":"[23:22:18][INFO]  Canonical SMILES: CC(C)Cc1ccc(C(C)C(=O)O)cc1\n","truncated":false}}
%---
%[output:5f53679a]
%   data: {"dataType":"text","outputData":{"text":"[23:22:18][INFO]  FP size: 2048 bits, on-bits: 25\n","truncated":false}}
%---
%[output:9465a5aa]
%   data: {"dataType":"text","outputData":{"text":"[23:22:18][INFO]  Loaded 200 entries from data\/list\/fda_drugs.csv\n","truncated":false}}
%---
%[output:5dac7e8f]
%   data: {"dataType":"text","outputData":{"text":"[23:22:32][INFO]  Parsed 200 \/ 200 molecules\n","truncated":false}}
%---
%[output:64d2b1f1]
%   data: {"dataType":"text","outputData":{"text":"[23:23:03][INFO]  Completed calculation of 200 ECFP4 fingerprints (Radius=2, NBits=2048)\n","truncated":false}}
%---
%[output:19813262]
%   data: {"dataType":"text","outputData":{"text":"[23:23:09][INFO]  FP Matrix: 200 x 2048  (1.6 MB as single, constructed in 5.48 s)\n","truncated":false}}
%---
%[output:28f71f40]
%   data: {"dataType":"text","outputData":{"text":"[23:23:12][INFO]  Tier 2 (RDKit Bulk): Scored 200 entries in 2.6634 s  (0.1 k mol\/s)\n","truncated":false}}
%---
%[output:86732921]
%   data: {"dataType":"text","outputData":{"text":"[23:23:12][INFO]  Tier 3 (GPU matrix): Scored 200 entries in 0.1063 s  (2 k mol\/s)\n","truncated":false}}
%---
%[output:4e351304]
%   data: {"dataType":"text","outputData":{"text":"[23:23:12][INFO]  Max |GPU - RDKit| difference: 1.41e-08  (expected < 1e-5)\n","truncated":false}}
%---
%[output:460419f8]
%   data: {"dataType":"text","outputData":{"text":"[23:23:12][INFO]  Validation PASSED -- GPU matrix and RDKit Bulk scores match.\n","truncated":false}}
%---
%[output:0ec05bfb]
%   data: {"dataType":"text","outputData":{"text":"[23:23:12][INFO]  --- Scale-up Benchmark ---\n","truncated":false}}
%---
%[output:086b41df]
%   data: {"dataType":"text","outputData":{"text":"[23:23:12][INFO]  N=   1000 | CPU 0.0045 s | GPU e2e 0.0133 s (0.3x) | GPU preloaded 0.0025 s (1.8x)\n[23:23:12][INFO]  N=  10000 | CPU 0.0045 s | GPU e2e 0.0725 s (0.1x) | GPU preloaded 0.0124 s (0.4x)\n[23:23:14][INFO]  N= 100000 | CPU 0.0291 s | GPU e2e 0.2807 s (0.1x) | GPU preloaded 0.0059 s (4.9x)\n","truncated":false}}
%---
%[output:5897ae42]
%   data: {"dataType":"text","outputData":{"text":"[23:23:14][INFO]  --- Benchmark Summary ---\n","truncated":false}}
%---
%[output:66892007]
%   data: {"dataType":"text","outputData":{"text":"      N        CPU_s      GPU_e2e_s    GPU_pre_s    CPU_kMolPerSec    GPU_e2e_kMolPerSec    GPU_pre_kMolPerSec    Speedup_e2e    Speedup_pre\n    _____    _________    _________    _________    ______________    __________________    __________________    ___________    ___________\n\n     1000    0.0044601    0.013265     0.0025326         224.2               75.4                 394.9              0.34           1.76    \n    10000    0.0045062     0.07249      0.012409        2219.2              137.9                 805.8              0.06           0.36    \n    1e+05     0.029056     0.28066     0.0058941        3441.6              356.3                 16966               0.1           4.93    \n\n","truncated":false}}
%---
%[output:42f6c091]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAeUAAAElCAYAAADa2PrWAAAQAElEQVR4Aey9a6wl13Xf+adzA11iPBnyQ4\/USTdYzVCAiQTIAxG6kfBxWqYMJZrEyQiJTNpSF00pkZAPrRiiJUeSb3VExW0Jk4gfBrZgy30oOpTzcD4ECayEcrpoEmETkwD5ktgBaXU5fTNNhEaoiRjkCr4Gp35VZ53ed\/euOlXnnvfZdft\/1lr\/tfZr1ana9TrV3\/d2XGIGYgZiBmIGYgZiBlYiA9+nuMQMxAzEDMQMxAzEDKxEBjZzUl6J1MZOxAzEDMQMxAzEDPTLQJyU++UrRscMxAzEDMQMxAzMLQNxUp5bamdecawwZiBmIGYgZmDDMxAn5Q1fwXF4MQMxAzEDMQPrk4E4Ka\/PutrMnsZRxQzEDMQMxAyMMxAn5XEqohIzEDMQMxAzEDOw3AzESXm5+Y+tb2YG4qhiBmIGYgamykCclKdKWywUMxAzEDMQMxAzMPsMxEl59jmNNcYMbGYG4qhiBmIG5p6BOCkHUnzvvfdqEgLFKsotVxETPoifENLbTZ0+mirx40K2X5YYnwvZXeNCZefJ0S8fTe35cSHbL0uMz4XsSXH4uyJU\/6w4+jCrulalHsbUFdZn4k1fZ7kp41jnddDW9zgpN2Tn29\/+ttrgf7GxgVsGGzQ0UU38Tb5peNoCbh9Mhwehei2mSTaVC9W1yhzjAKFxwoNQ\/0PxLtdULlRXH85tw3TKm+5K+IjuGXByN97OKd3E44uIGVhEBuKkPGWW2Xj9nTGcW51vm49ywOxZSOqjPRCqr4kPxfocZanf59fJpv+MA4T63cSHYn2OstTv85tiM75ljWWT87qsnK5Cu3G9Nq+FOCk356azhy9Y044LHr9bGZzB5eet06bfl3m3uU71x\/ys09qKfW3NQHSubQbipLziq67LJEoME8qKD2Vp3Yv5WVrqY8MxAzEDPTMQJ+WeCVt0+Kwn21nXd5x8MFm2lXf96CG0lZ\/Gt0r5mab\/lOmSJ2KIBegu4FzgC9nwLizG5dCNNwnnw3wmzY\/t6tgG411pvnlKtz3T\/fbgjUM3hDjfZzFIfMgQQj44F5QzG92F8a50\/ej4TKIb4FzAu7aruz50gB8J0DtiK8LipDzlaubLtAk78C7Dn9dYyR91h\/oAjx+f6dg+8BGzTNAH+rXMPljb1hf64wLeYlwJ78ahw7kxIZ0YYl3AAZdDh7M60OF8wFsM0vy+jg2ItxhXwuOfF6jfbc90+FCb8BaDJMbn4AE8\/mlBeepx4XJuvS7vx7tx6KFYOHx9YW1RztWxI+oMxEm5zsNtn3zp2sAX6rZCa0q0jRPfKo91Fn1jjC781ej6Qvos+uC3OY1N3\/r0pW+826emdpp4yra1Rzn8xE0CccSH4uDxh3yz4Ki\/az30w48PcVYfsfjNnqekra71h2LhFtXXrv1c2ziv43FS9hJiJl+6EPDDI9cRoQ2J8bRhnuOkXb9P2PDWLjqcwfhZSeo3hOo0X5MMlVkGR\/9C7ZK3EN8UH4qNXP8MxPz2z1ksIcVJuee3gA2taSfXs6oY3iMD5N1A\/g09qtiKUMuLK8nbVgw+DjJmYAMysCGT8gasiWMMgZ0uO+FjVLFWRRmvocu4ie0St1ZJCHSWMTJWH4HQSG1gBljvfAdcwG3gUDd6SHFSnmL18kXni29Ffdt4JHH40ZeNVeqLmwvyQ9\/gkNjoBjjT5ykX1c48xxDr3o4MhL6rcGw7LpadDfq07D6sW\/txUp7hGvO\/gL7dt6k+8WyItAdC5eCJCfn6ctRDfSBUFh4QF\/I3cZRp8h2Xpy\/UD0J1wRMT8vXlqIf6QKgsPCAu5J81R1uzrnOa+hhvU1\/g8XeplzjiQ7Hw+EO+VeDoG30M9QUef8jnck1xlMXnwy07jU59fjk42vN51+4S48ZHvc5AnJTrPPT+5AvJl84KYgM4AzawmGkkdXUtR1uAMj7gu9bTJY76gN8ONjzoUo\/FWLxJ45Fw1BsCPmK6gFhw3Hpm0Rb96FJPnxjqXMTY+vTJj+3bRz\/e6vN5Gze8xayqpI\/WX1fC+32Gc2PQ4fw4bPPhdwGPf1pQF3W4gHPrw3b96HBujKvjI8bg+rZdj5Ny4BvAFyZA30aF4uAMtxUIEMQG6DE1yT8OdBTK+HDcR1TijhA9Dcr7aK+i2Us9TV58ITTFt\/Ft9eBzy\/q26+uiU95Hl3JNMdTV5IPH78N4pIEY033p+ybZVt6Pa+KJ82GxIenGun6XN93199Wpo61Mm9\/3+bZfL34ffozZTXHwFsPk5trGh2RbXJMP3kVTvaEYuEnxIf+2cnFS3tY1H8cdMxAzsFEZYGL2BxTi\/Jhor1YG4qS8Wusj9mZNMxC7HTOwzAxwNgqYhF3AgWX2LbbdLwNxUu6XrxgdMxAzEDOwshlgAnZx3I5S13HriOX7ZSBOyv3yFaNjBrYoA3GoMQMxA4vOQJyUF53x2F7MQMxAzEDMQMxAQwbipNyQmEjHDMQMbGYG4qhiBlY5A3FSXuW1E\/sWMxAzEDMQM7BVGYiT8lat7jjYmIGYgc3MQBzVpmQgTsqbsibjOGIGYgZiBmIG1j4DcVJe+1UYBxAzEDMQM7CZGdjGUcVJeRvXehxzzEDMQMxAzMBKZiBOyiu5WmKnYgZiBmIGYgY2MwPto4qTcnt+ojdmIGYgZiBmIGZgYRmIk\/LCUh0bihmIGYgZiBmIGWjPwLpOyu2jit6YgZiBmIGYgZiBNcxAnJTXcKXFLscMxAzEDMQMbGYG4qS8Sus19iVmIGYgZiBmYKszECflrV79cfAxAzEDMQMxA6uUgTgpr9La2My+xFHFDMQMxAzEDHTMQJyUOyYqhsUMxAzEDMQMxAzMOwMbPSnfe++98jHvhPatn\/71LdMUT11NsDL4TV9VOamPrt\/VFzqeQGP0xYcbhs+1o94tA+TNh1\/S95vtxsG5tqs3+eB9uOUWodP+rNuhzibMuq1YX78MbOykzBfu29\/+tnzA90vR+kS7Y6XXvg0XMZ8M8L1y8206\/Hxa3I5ayZ\/l0pXwfgZcv+mhOL9ck01Zq8eV8E1l1ol3x+TqxxnDpuTmODk4btmNnZSPm5hFlWdjWFRbtLPo9mhz1ljXMaxBv2e9qo5VHzv4ppzB4z9WA7FwzMAKZmDrJmU2Znc9sGG7MB8cOhK4OjaAM2C7MB4JbxIdYIOQDmcgxmCcK83XVVLWYk1HApcP2XDA4nzZ5DMe6cIv39WmDj8WzmA+bHQkcHVsAGcwG2kwnyvxuXab7n7frBzSh9XRxJs\/ylsZcHN7i52\/5rbL+qJFpAEbmG0SzmCcSeORxiGxXbRx5kMa3LJ9dCtv0i9rvEn86CZNd204AGcwGwmM32a5sZMyGw4r2RBayfiIcwFnsejmc3WXI7bNhx\/4MdjwPuCtfiQ2MUhsAxw68jhw60UH1AvQAboLuFCbxPg+bOORLvCF6unLUU9Tva7P1S0ezm0P23xI7JDf5dAtlngA1wRiXVgc5VweHc78UU6ZgSmLWf5ZB6CpGnzEGrCB2SbhqANpnEm4kM94fF1AvNWJxG4qh88HsXCUdQGHD6C7PnTj8GMDdOOxDXD4DNjmM26b5cZOyqxUW9FIVjyAbwOx5nd146aVXevqGjdtP\/xyfnu+7cfP0p5VW231uD5XbxpHl5i2spQHfNdAU6zxxBBvdpTtGSBfLtxolzf9uLmlvMHqdNtEx4+cB\/rW3SeeWB+MAQ7ZB9OUsfqPU9bq2CS50ZOyu6JY8YANy+X76JR14ZZ1eXTXt+4643ExzXgs98etZ5q23TJu++iub5Y64wVtbeAjxm0XzoXri7qOPLjp54Nc+vBjjmNb3ayfLvUQ58It4\/Lorq+jPtcw+uRiro3Fyo9kYGMnZb5QR0ba0WgrZxulK61alzPdfLOQ1EnfDNizqLdLHbTlo6kccfQRPxIbHaC7wA+\/KNCe2z76pLaJoRxxSGx0H\/h8rs0mPlQXnI+2eqJvfhlgHR2ndn89Ylt96D7Mt2zJuPv2jTLL7nehQufLP635srGT8iLWS9sXsc03Td+oz91QpqljVmXoS9+6pinTtY1p6562XNd+NcXRLuuyye\/yxLr2NunkqGn8TXxbfprqoy58bWXbfJRv85uvLa7J18Q31Tkp3sr1kfOos7H9Yzgu6ZLy8m\/dJ+aNnZTZyPgy+YDXaEFv84\/CKtEW2+arCs\/ow+\/rjKptrWaasVkZpFWO7vcfzvy+9GPN9uOwqcf8SGx4H\/D4Ddh+TMgmjjLIkB8OHzE+4PGrKHT24EDK80o++of\/sFzAf\/uXfum2l92My1eVbN8H4\/dzig0\/TTYoR3kXcKG64N040+HdeGzzIbEBugs4yiFdHh0u5DMeH8Am3oANb8A2HxLbfF0lZShrwHbLYpvPJBwxSON823hi8M0SQ9V\/Khcm5kyZ1nXZ2EmZFcLK9wHvoskP78ahw7mAM7g8uvFI33Y51+fqxAA4vsxIH\/DEhEDsJN6PabLhXYTq9TniQxy8wfebbf6QJAYeCUxHGuABNtIFnAEeHQlcHbsvKO9jXEc5KX\/j5k3p\/HkhfcADJma3jnH5FVD4voUw7665+TCdNtGRwNWx20Csiz6xlAvFwxvMb7ZJ45HGmYQzGIeEM4kOsA1mIw3mQxrnyzYfsfgNZiMN5jNpPNLnzDZJjAHO9Gkll60f1+NHittZ8xFyeUavljd6Uu6ViRUO5osb2hnCr3C3N6JrlveZ5XowkPb2bsdgsBb5Ig8+1qLjsZMbmQEm5PM6Px7bQLe2I3+i1poscVJekxXl7wix16Tra91N8gxmNoiHH5ayTMoyKcukLJOyTIKfWSOxokVlYKbfjWN2epX6csyhdC7OGTETMwVSpbpa\/tnEDO9O2FqTZW0m5Xnkc39\/XxExB4v4Dtj3d\/+\/\/JfqO6fhUMoyo9dG2pUDZFOnF5HP2Ebcbuu7yOV2VH4Rd\/Z39IX9L5Saymn5qmxizpVr3e4vb+2kzEb95JNP6qGHHoqIOZj7d+DRxx4Tyyu\/8ivaf\/e7pccfly5dkvIcem3A2ZghNDHH7SruTxaxT\/2zj\/1Z2eVpJuQTT56otuFHH320OujdU3mLaLRVcTbN5DwyV15s9aT8yiuv6Mtf\/rKee+65mYN6L168OKHe29vtWm5SXJO\/K+\/H+TZj49sNv035Y6yMmfGjGyZxxJOvD771ls7xFHZp7O\/sSEVRauvxj8nY7Sm2PzEzKW\/qdhVax13Xf1OcW6d9R+AsftaSummnb71dy7XF9fWF4uF++OIP6+ZzN8dfxVSp\/unFfyrGxXdP5TLQQFfKP40Wm8BH5kqLrZ2Uba2cPXtW586dmzk+8IEP6C\/9pb+kP\/Nn\/kyvuruWmxTX5O\/K+3G+\/cM\/\/MP6yEc+IuQ25Y+x+rlo4\/7PP\/2n9Wf++T\/Xuc98xr5yUpLo2vvfr4dOn9b+I4\/c4jdI28TtKrTez432HSFfF86NYVua5zZFX932sLuia7m2uL6+UDzc2x956gvB\/gAAEABJREFUW4enDqut5dThKf3t3\/vb1X72rQ++pYOzBxXPB5N1okQs3F\/OlGkdlq2flNdhJa1iH9\/1rnfpwoULwa5Fss7Azv6+Tjz5pHa++MWasM\/BQKdOndLFN9\/UqV\/8RSnLpBdeMO\/KSv+sGJuz5ZXt8Jp1LG5T3VbYxe9c1KlyMiZ6f2df13aviUn3c6c+p5vfuKkHTz0oLldzzxle5cLkHCflMhHb\/G93d1dnzpzRzk55ibJHIrqWmxTX5O\/K+3G+3WNIU4VO217XcpPi2vwh3xGuKKQs005573i3vEVyWwKGw2oyvvid76i6r7wm95aZgJmIDdi3jW3OxJE892ira7m2uL6+ULzP+XaPIU0VOm17Xcu1xfX1heKNe3bn2fH4Hzv5mL64c+vAl4kY2CVrJmT3Uva44Ioq8Ux5RVdM7NaqZaBDf4pCKifj8mhM1WRLkSSR+F3y229LV6+Oce3yZT168qT2ny13Lg7PpW2KrSqYiA2r2sfYr83PwEAD7ZV\/Gi1DDWULk\/AlXTJTF8o\/4sfEiitxUl7xFRS7tyYZyHPxRq7xZEy307SehJmosQcDaTCQBoPqd8mvlFdTKn0wkAYDaTCQkoTIiJiBmIEJGciUaVD+KbBwpgzNfWXi0NcFcVKe05o6ODjQ9evXdXhYP5DQtZmu5SbFNfm78n6cb3cdz7Rx07bXtdykuDb\/EV9RqPp50\/nzUlHUw02SajI++Lmf0\/XyDNm+A0fK1ZFL\/1y3Dkybw67l2uL6+kLxPufb814f07bXtVxbXF9fKN7nfu3g12T3ly13NiEnSrROl601WuKkPErErIXd+4j3lKfL7Krmj9FUfbvjDu089ZSqS9XDIbSUJPWl6vJgTIOBqjjnuQLfrgvFzz4ZmDaHXcu1xfX1heJ9zrf75GKa2Gnb61quLa6vLxTvczzkxcNefi7WdUJmHHFSJgsRMQN9MjAc6rZL1Xt71dmx7FJ1n\/pi7IwzEKvblgy4947dMV\/QBTVd2taKL3FSbllB+\/vxVXZNOXj99dcFmvwby7\/vfdr\/2Meqtwbx8o8Kzz6r\/Y9+VJXe8TtjX7tNy5ONq01u2phnNR62JzCr+ta5nrbvj\/m4TJ0rl7+s431kdwxxUnazUerc\/wN8oeNrOB+qXl0Xem3ee9\/7Xv3oj\/6okCH\/xnK\/\/dvVSz948ccYP\/3TjXlqysNjj9Wv3UQ2xawjz3jKzajxX9yupt+m1vH7MG2f7XWZjV+k0uGeJXO5+uVrL2v32u74HdhlyFr+i5Oyt9pu3LghwM6DV7bxWre+r6SL8be\/PjTmZDtywqsOvU3qiBm3q+34Hhxne+c7xL7Xf6iLL5Jxrx2+JvdnUHap+uRjJ7XuS5yUvTV48uRJAaPn9bpA\/\/V273rXOf3Wb53TV796Tj\/1U7X85jdn\/\/pPv91oxxzP8jvA9mLbTpskbpbtxrrW9Xt8e7\/5brR9d\/D5l625hwy\/CYiTsrcWd3d3BTx6rmae188N8R8HDYdSnkvDoaqfvJ45I+X5XJuPlccMbF4GhkNVP1Vjo2rDcLh5Y9+QEbEf9t+KCHfHmTuOvMFroPpPG7LESXnJKzLP6wm5KOqOJIk0GNQ6n0Wh6kHfPMeKiBmIGeiUAd4lPhxKw6E0HErDoTQcSsOhNBxKw6E0HGod3jneabxbFJQrFw95abTsaW+krY44Tk\/ipHyc7B2zbFHUB\/NUkyS3fuLKWxffHr2VER\/gYB85T9g7jU1aW2a70vWZbpI4048r2+pq8x233Vh+QzKQJNJgIA0G0mAgDQbSYLDQwfE9dWGNu5zprs90k8SY3kUSD7rErlPMM3pm3N1EiQblnzZoiZPyEldmUUhFUXeAn7n6P3EdDKQrV2p\/UUjDYa3P45ON195pbBLO2jLOpOuzmEVL+rLoNmN7a5aBJKl\/P86RrovBYCEDYTvhe+oCzhp3eXTXZzHTSOqhPoA+TR3LLmMPdfFrGOtLrvpPo+XTB5+e6s2Jo+IrKVZ3Ul7JdElMnLOCe\/ZbFFKo3qK4lYhnygPEUMw0XFHcqncRGjuGENy2zW8cNnpXSSywePSINckArymdJfK8Hnieq7r\/49ddFO1+P76rXRR1vS2fTJQt7t4uvu+ga0FiQdf4ZcVx\/9i\/p+z\/DOrjux+XH7Os\/s6q3Tgp98wk\/8verFAUtxpvq9Oi8lzVw19tsV19VmebdHcebMQuXF9bHa6PMgAOCagTG4kN0OHQXYmOz3hsgA2PjsRGj1ijDOS5lOdSnkt5LuW5lOdSnkt5LuW5lOdSnkt5LuW5lOdSnkt5LuW5lOdSnkt5LhXFrcHnuZTnUp5LeS7luVQUtb8opDyX8lzKcynPpTyX8lzKcynPpTyX8lzKcynPpTyX8lzKcynPpTyX8lzK87rOhk++mwYLMdtk3+8u5SgD0K1ebNNdSQw+gO76Vl3nPjLnydbPgQbaxCVOyotdq2vVmrvRshG7mPVAqJv2QFvdxIX88JRFhvyRixlYdgb4bgK3H9guXF9Xne89CMXDU7\/rgwMutw66f5a8t2EPeNk6iJOyZaKj5AGsWcHuF9M0t7tC9box6KGYabgkodXVATsJdh5gml5ZeeQ05WOZJWeAhypmicGgHlCS1E9Q+nUnSbvfj+9qJ0ld7wI\/2WYMbrNsC\/Auhw5nwF4HcJbsvyyEh7zWoe99+xgn5b4Zm2F8mkpJUlfILSvuDdeWVBQStt13ThIpTc07e8lGykbsAm5SS8Tce++96lvOr9etx\/UZ73K+TtvEwSOx0SPWKAN82WeJJKkHXxSqfvbET6Rc1F5pMJBm2a7VO5L2feQ7aRi5WkWoHJxfCM7qRfp+OGD8pHiLWwXpPujlXrambxd0ASE3piI24CNOykteiZz9Whe4H3zmjKrnUpDY+JLk1lPY2PMCG6wLawfO9JDE76Ipxnhim3R8wPxIs03CAbNNwgHfhovY4gzkuYoi1\/m9XMpzKc+lolhYQvg++qBxOGQT8LvoEufGNJVt4t2yq6Rzluxeuj57cFYPHD6wSl2caV\/ipDzTdPavbDCQuHRtJYtCynOzpCSpr74NBre4qC0kA7GRdc7Aww9LaSqlqZSmurQn5eU2dP5qOag0ldJUSlOJuJKK\/1YvA\/b09Us7L4mJWaPl8u5l7ezsVJbFmF2Ra\/4RJ+UVWIGDgcR94b09CQwGUprWZ8fXr0tpugKdjF2IGVinDKRpvQGVl6KGVx7WsDTpPhNzduWesU9pCh2xwhnY9JeF+KmPk7KfkSXaWSZlWX3mfOWKlKZL7ExsejMzsGWj4gzrcT1+ZNRcCs2VKy6rnwHWE7Ce7qk8czFjQ2WclJe9YotCynMpz6U8l\/JcynMpz6U8l\/JcyvNl93Kl23cfZFnpjs6gc9s01uOmiwn5vM6PqxloIFv8idr4KFcnAzzE9ZmDz4w7lCjRjxz8yJE3eBFzvbyc6L71S2u+xEl52SuwKFQ92cXj123I87n3lB2+C2vQ5Ux3faabJMb0LpJ40CXWj6EcD64gfV8Xe1K5Sf4ubUyKoQ0wKQ7\/ccZK+Q1FcFicETMx40yV6mr5ZxMzvDthKy4rl4HXd1\/XK7uvjPvFuvPvIfv2OHiNlTgpr8rKGwxU3VDmprKLwWAhPWRSYIfvAs4ad3l012cx00jqoT6A3qcO4ilHGZPo6wQbQ5\/+E0u5dRrnovvKb1oB7XKGdUXl\/aDScCdmLotmKu8Xab4L68qFteZyprs+000SY\/o85XHb6VO+LZaDKhsn63BPm3\/pmvFu\/KTcttJJwMrg4YelLJOyTMoyKcukLNMynw5l5z\/L\/LAuQNc6iQVd4i0OafDLGY\/0fSHb4pokZfABdIBuwHZhPBLelegGfAY4dJPoEe0Z4CzYLk+zM7cJ2UrtOTt3dvxMzuabtWS9sR25gLN2XB7d9cmCtlS664WzZNblNqRioyfluXzB+QHxLGFvB+F\/mwjVaz9WJi7kn5YrisbvN3kzWJDZJtmBmK+LpBxlALqVwTbdlcTgA+iub5JOGeCWQ4czYHephxjKIAHlzDYdGx0\/wAYuhw5nsDhf4icW3oANb\/Ymyjs0u78zOjNOERM0l6nd2rHHAaWC7fqPo9NeWWXrv1muS74bBmsUG92k6a4N1wWUMbjxx+GoJ1QevgkXdKHJtXH8Sk\/KrLhQxuFByOdys\/zyj+stCqkopKKQikIqCqkopKKQikIqCqkopKKQikIqCqkopKKQikIqCqkopKKQikIqirrqopCKQioKqSikopCKovbxWRRSUUhFIRWFVBRSUUhFIRWFVBRSUUhFIRWFVBRSUUhFIRWFVBRSUUhFIRUFNTaCvAE3ANuF6+uqs85AKB6e+l0fHHC5Wel+W13rdcuh0z8wqXyXWOoBfl2U9blor18GWLcG673ZJqdZ15QB1OHWC4cNjw7Q4QC6D3gDPsoYsPEh2zhiQChOUvX2PytP3CRwlgyI8x\/s8m1i1h0rOSmzMkEoufC2QtGJQbqAmxvSVEpTKU2lNJXSVEpTKU2lNJXSVEpTKU2lNJXSVEpTKU2lNJXSVEpTKU2lNJUGg7qrg4GUplKaSmkqpamUplKS1P40ldJUSlMpTaU0ldJUSlMpTaU0ldJUSlMpTaU0ldJUSlMpTaU0ldJUSlMpTes6F\/xp6w3pNs268zn8cAbsWYI2j1sfdfTpn8VSLtS2+ZEh\/yw5+tCEWbYT67qVAdYruMVI2C5c33F06nTL27p2OWJ8uP4mnTKh+kLxXeNCZY3b03bcS7bxruSkzEoH1kmTrGCXRzcO3WDxc5H8gHiW4KEuOso95VC9F0aXbZAh\/7RcktDqGOSOXLoYO1uUUDk4vwhcW93ms3KT4omzGPRJsFhrB9svg8\/niAvxxJmvyU+MgRgD5Yw3CWd+pPEhiZ\/4kG8SR1lA+SbgB5PqmqWf+76zgr8T58lrt27Xnygpd\/l7cv3H0ZOyPq3g4q7r43aP74bVN6kui0NOig35yedAA9niP23t2xa3cDnDBldyUp7h+CZWxe\/bQphYcNYBvCw\/y6QsU55nUpZJWabqZfqzbquhPjYcH4TCIZuA30WXODemqWwT75e1nQQ8ZZDA1c2GA9igScdnsBiTxiPhDGYjATwSoBuwDXCujg1cznSkO1bsPrCybv2h8vgB8SH\/PDgmzlkhU6a98k+jhQeGrO6BBnLfEHVBF0S8+Y8r5S2WR3Jp8EKCZqgcXCjY6m3yw1sMMlRHiPPLYROHpB6ADSZxk2KpIwR3PYb8m8ht\/KTMl6VtxQ2HQ924ceM2vP76623FZu\/Lc+nSJeUvXNL5QS1V2srz2be1YTVOWsebNNzjjLVv2bZ4drJgVXObKdOg\/FO58PAVD3OVqnjaGhudCZg49HmCPPqgPThkE\/C7mBRnfsqYbhLOYFyTJM586AbjkMYhsQG6ARuYjcQ2YBuM8yVnyawjn990e4Um5eWk+pFHHtHJkydvw1133bWYDmKAZb8AABAASURBVA0G0tWrFYqrV\/R4\/XPKWl6t+cqfJIvpT2xl4zPAZNqGTUnAVV0VO3aVS65c\/ExqqKFY4LlMjR6xmhngf4Py39blP9jl26s5kn692vpJOUkScV8ihH6pPEb0YKBikOjxwTMqkroeJGfMGgykwUBKktoRP2MGjpkBO0NpkkzYbU3gp2xbTOiWEFxbmXn43Il3qM2bkCetB63x8tGdj+r06dPVCPjuAP43KJdz7SpwAz7WalLmC8gOwfKODmf2LCQrHsyiLuroCnYYuXIl5d\/V8ghf5YK9iMtrZVPxX8zAOANsU2xbY2IKZTgc3nZLiNtEi74txCXsPR19epf7yPBTDCsWWWAG3nHzHcHvEN+jEL75zW8usHfza2qtJmXSYDsMdhrocLOErexZ1jmprly5uNfFhHxFV8QOwyZmePyKy8pmgO\/iynZuxh1jrP52B+c3M+m20P7+vhaFj+5\/VI\/sP6Kd\/R0hsRfVdmyn\/3q279KJEyduu60YutVo3P33329F11qu9KTsb\/yWaXhg9iylreBZ1tlWFw+dcK+LGPcInonZjvDNT8w8wc7VhbXlcrV+r7mqFwGMjZFCzEidqzhuO33KN8XCz+u72Cd59KNPfN9Y6gehscLhc+tMkvBtofvuu09nz57VY489poceemhh+O2HflunHzot5CLbjW31X8d8N\/iO8F0J3VZs4u655x73K7i2+kpPysvIqq3wRbRtEzKSSdi\/VI0Nj9+eHp1Xv9ipsnN1AWftuTy667OYbZPkgFxs07gZs4FxoyPJg+nYTTh16pS+\/OUv67nnnpsKv\/RLv6S\/9\/f+nr7+9a\/3Kt+1XFtcX18o3udc29Wnzc+kctO20bVcW1yT72889zd087mbY\/zMcz9TfUdCD3H5nG83fe\/WiY+Tsre2uJ8MPHpsDjXbP7s0\/bAeDtYMr3IhjjPmWbVeVjnxHzvaiUEdA9hhG6wINrpJ010brgsoY3Dj+3DEumXR4QD6JBBncGONQxqPHkLI73JtOvW5ftc23pX4Dcb7tvFI+z4gAdw0YGI+d+6cpsWf+BN\/YqqyXcu1xfX1heJ9zrVdfdr8TCo3bRtdy7XFhXzfPPdNHZw7qPCuc+\/Sh859SHxHpvlubUKZOCl7a3HSPWUmxlmB+8XWPHqoXniLYUIOxUzDcfZt9frSdsxI86G7mGanTBlAPW69cNjw6AAdDqD7gDfgo4wBGx+yjSMGhOJ8HrsNbh20iU08EtuADQ+Mc\/WQ3+WI9UE9cCaJRwfo+FzAgSaOcsCPId7nseEXBa5knTlzRjs7O72a7FquLa6vLxTvc77da1BTBE\/bXtdybXEhH\/shTjhsKHbLDjsU73O+Tbl1R5yUvTW46HvKXvMrYbKjBW5nsF24vuPo1OmWZyIALkeMD9ffpFOGukBTjPHEALOPK2n7uHUcpzxjAdLttVjfkAai0JERMQOLygAnFNZWokTb+LIQeUuclL2EcOQFPHpsXtds\/7hnTOV8If2aeQLb\/Ejffxyb9mh31cDEYDhu35iUutZlccjjtkt52kYuC4zDsKw+xHZjBtoy4J8ls49ri98WX5yUvTXN\/WTg0WOTyWyW4CdQVM4XlMvTbt1c1gFwxCFnBdp0wQ6cicSF62\/SQ+XgQvFWd5Mf3mKQoTpCnF8Omzgk9QBsMImbFEsdBqsL23TKA+w2Ht8kUA8I1eWXJYZYeNOxAVwI+AyUCcWEOGIphzS\/bxs\/a3lwcCDe8tS2jYba7FquLa6vLxTvc74d6vssuWnb61quLc73sb+zsZ06PKUfuvlDcterH0+sz\/k2MeuOOCl7a3DSPWUv\/NgmkyxnxFTE\/WMmYXQkNn6bkOHnCXayPmgPDtkE\/C4mxZmfMqabhDMY1ySJMx+6wTikcUhsgG7ABmYjsQ3YBuNciY8JCQ7dgG0wDulyXXS3DPHYBmwX8GajG4wzSX\/NZxIOPzYSNOn4XFDWjXV9s9a5ihXvKU+f1VXJHych7N9sJPft3KdHTz565FmBUF99zretvnWWcVL21t4y7ilz2cYecOAeC19YJF1zf7uMHbF6GVjUhLQqI2cSdvuybeN3xz5bfXtq46TDHa3t\/1xuW\/U4KXtrniMv4NFzN93fJJ\/RGTExM1nDawOWuOPutxJXOV\/0DfQbUYyOGbiVgWf0zNhgPwfGxJYrcVJeoS8Al6mtO3bZ2uwoYwZiBqRp7yF2LdcW19cXivc53573Op62va7lQnE2JvN96\/Bbcs+UOUs2X7ynLM1tUuYSVxO0Rssi3127s7+jZ\/efFfIL+1+o5CLbj231f0\/vOuaMzc\/vd4jzY7rY1BMRM9CWAfdeMicf8Sz5aLZmPinbRMzlrSZYzNGurIbFkRo4derUUt7R+9MP\/bR4Ry8yvje3\/3tzY84m54wtzc9TiPNjuti8t5i65gVuLcUHvabP7rLzd8eZO\/TSzkvjAXCWjBHqVxcuFEN988f8WpjppMxkaxNxW5cthvi2uGX47OlrJmV7Ry\/v2vXfKQsHXN638cEBdINvw8MBdINvw8MBdINvw8MBdINvw8MB9Od+5mf03Dveoedu3qyB\/jf+RvWOYWJAFffcc+I9tj\/1Uz915B3E+IHFILEBugEbmI3EBugGbGA2EhugG7CB2UhsgG7ABmYjsQG6ARuYjcQG6AZsYDYSG6AbsIHZyI985CPV19u+Y3CbgIsXL1bjih8xA6EM+GfJ8WUht2dpppMyk+3tTTQzfeOba5qdx336momZ98i+\/\/3vF9IFHPA510YnBqAbfBseDqAbfBseDqAbfBseDqAbfBseDpx717t07jOf0bnf+i2dOzhQZf\/ar+ncpz9djf39ZQ4AZQwPPvhg5TMbPzAbiQ3QDdjAbCQ2QDdgA7OR2ADdgA3MRmIDdAM2MBuJDdAN2MBsJDZAN2ADs5HYAN2ADcxGDgaD6kvL\/4aDvSlgPNXA4kfMgJcBHl51f5s8UL0NKC5HMjDTSflIzY7hnhGjA0lOxOqou7u7Am6PfBsfHEA3uLY9uLCzs9OrPivHJXS3Po0WODApjvCb5Rkv9aADv1zFvf66di9fls6ckYoCStrb08Fv\/qauv\/32kR\/z+\/W9ozyTrgvUn9QPaqv+xAa1VX9ig9qqP7FBbdUP9NAe+TMOSQxAN2ADbMsL5YyDB9gAfVIcMe5LKigH4ENlu3AWwzqxuqhvWXC3Q3Qfy+pXW7tuDtvifF\/Xcm1xfX2heJ\/zbb\/fs7anba9rubY4d0LmXvIFXRgPL1SuCxeKGVe6psrcJ2U2dDsjNh0bfU1ztjndLgrp8celS5eqMR2W99GZkJVllR0\/NjcDbH9sh+4IsQGcSfSImIHjZoCzZP\/SdTxTDmd17pNyuNkNZkdD40xoWQ+k0IWm9sf8Sy\/VZ8d5TriUJNp58UUpyyp7HFee6UP4Ntw8MW17XctNimvzh3xduFDMPHPYVHdoQnZjmZCJcblV0afNYddybXF9faF4n\/Pteed52va6lmuKc38CxRj3tIcYI1SuCxeKGVe6pkqclL0Vx6VF4NK+jQ8OoBt8Gx4OoBt8Gx4OoBt8Gx4OoBt8Gx4OoBsquyikLJPOnzdah5\/9rHT9upQkFVfFVdqtDzhwi9GRS9vw+AG6ARuYjcQG6AZsYDYSG6AbsIHZSGyAbsAGZiOxAboBG5iNxAboBmxgNhIboBuwgdlIbIAeETOwjRlwXxaSKFE8S1bjEidlLzX29LVLw7k2OhxAN\/g2PBxAN\/g2PBxAN\/g2PBxAN\/g2PBxAN9x8+WVVk\/HocjWT8JsXL+rGE09YSCX9cpIEB6qA0cfr5f3okVoJ\/KAyRh\/YYGRWAhtUxugDG4zMSmCDyhh9YIORWQlsUBmjD2wwMiuBDSpj9IENRmYlsEFljD6wwcisBDaojNEHNhiZlcAGlTH68HM3opcmODteWuOx4Y3OAGfJwAbpviTJuChvZWDukzIbO5fCADpNuzr2KsF9+tr6BWe6SThgNtK17QGEEydOyOX9OGxADLBynFlh43MBBybF3XXXXfre975362w2y3T6oYekoqirSxLp6lXdefly1T+\/PtemwG31lSRjK8X4H\/0CY6JUsEGpjv9hgzFRKtigVKt\/tE\/\/+7ZBOR7QopxbH5ViA\/RJcYyXelgPxFMOoIfKduEshjqpn7oi+mXAzWGfkl3LtcX19YXifc63+4xpmthp2+taLhTn3ks+dXhKDxw+cFvXQ+W6cKGY2ypfM+L7FtFfJmNgbbm6casiuUcB3P74Nj44gG7wbXg4gG7wbXg4gG7wbXg4gG7wbfg777xT1dPRRSGVE7I9zIWPs+PDV1+VkqR6Mtwv79uUGdeHMUJV\/0hHUA6gG7CB2UhsgG7ABmYjqd\/nsAF+AzYwG4kN0A3YwGyVCjYo1fE\/bMY7JkoFDpTq+B82GBOlgg1KdfwPG4yJUvHtkor\/YgY2LgM84OWeJX\/wux\/cuDHOekAznZQ5A+7Twb7xfepediw73aU\/6HXHHdp53\/s0npDLSZinq+\/+ylfET4bcHPn9nWS7Zeeh++13baNruUlxbf6QrwsXiuk6rlnGcVDctu3hI2aWbc6qrmlz2LVcW1xfXyje53x7Vnlqqmfa9rqW8+Pcs2TuJX\/l7tv3PfTVL9eVC5Wj7DpjppMyGzIbNGhLCn5AfFtc9B0jA5wdu789ZkIuL1cL\/hjVxqIrn4FOHWTbYxsMAV+nSmJQzEBLBgoVcn+bHB\/uUqdlppMyLbJBg9DGbhx+QHzEjDNQFBITr\/MwF2fHcp6unnGLsbo1zQDbYAhrOpzY7RXLgDsh07ULuoCImJCBmU\/K1l5oYzfOYjZZTvsAQtdywbiikPt0dehlIMFy5Yrw+Ul2WWSu\/\/z2uzbWtdykuDZ\/yNeFC8V0HdfWx40SMG0Ou5Zri+vrC8X7nG+Phjk3MW17XctZ3GuHr8m9dM1Z8rmDc+W5wXXxoKM\/QCvn+rpwoRi\/7nWz5zYpr1sirL98KYDZSN82zud9e6FxnB07l6sPzp3T\/rPP6vDHfoxuVJh3\/6gfVI2NPrDByKwENqiM0Qc2GJmVwAaVMfrABiOzEtigMkYf2GBkVgIbVMboAxuMzEpgg8oYfWCDkVkJbFAZow9sMDIrgQ0qI37EDGxBBtyHuxjuno6+LAQuIpyBOCl7eeH3pMClfRsfHEA3uLY9gMD7m12eWN82Dt7K8SAWNj4XcGAct7+v737qU\/If5rr53HN6+557jjzQdaTc6E1dcGBcX8n7Nu2\/8cYbVV30Cxv4v7WlHMBnwAZmI7EBugEbmE1\/aIv8GYckBqAbsAE25XjAjnLGwQNsgD4pjvFSD30gnnIAPVS2C2cx1OnnjnojJmfAzeHk6FsRXcu1xfX1heI9rvr1g\/s9u9Xj+Wh++11b6VrO4v7Bzj8YV80DXpwpm4\/v\/9g5UkK+LlwoZlTl2oo4KXurjt+iApf2bXxwAN3g2\/BwAN3g2\/BwAN3g2\/Cw4HxcAAAQAElEQVRwAF1FUV2uvvvppyuTnzjpyhWpPGsmBtSO+tO3YeEAusG34eEAuoHfApuOxA\/QDdjAbCQ2QDdgA7OR2ADdgA3MRmIDdAM2MBuJDdAN2MBsJDZAN2ADs5HYAN2ADcxGYgN0w7J+p8xzHdaHKGMG5pEBzpKB1b0Xz5ItFZ1knJS9NHHkBVzat\/HBAXSDb8PDAXSDb8PDAXSDb8PDAZUT75H\/2WkwqF+VORgQVh2BV3GVVX\/4NiwcQDf4NjwcQDfwO2LTkfgBugEbmI3EBugGbGA2EhugG7CB2UhsgG7ABmYjsQG6ARuYjcQG6AZsYDYSG6AbsIHZSGyAbvBt4+ctea6DiRnMu61Y\/3Zm4FJ5N9lGzllyqlRx6Z6BhUzK7AAMdA0ducmY9gGEieWKQmJC9p+u5udOTkKb6unK+3G+7TQ1F3Xa9rqWmxTX5g\/5unChmLkkb0KlTMyA7RBMCF8p97Q57FquLa6vLxTvc74972RP217Xcjzg5Z4lX9CtJ67b6gj5unChmHnncN71z31SZqNnBwBsMOjwZkfZMQNFUV2utvvHPF19+PzzEpN0xypiWMyAZYDtEMRt0TIS5XEz8MWdL46r4Cw5U6YNWRY2jLlPygsbyYo1xOXJaR7gCJYrCinLdORy9d6edm7c0M599wVHHqynjOzK+3G+XVY113\/Ttte13KS4Nn\/I14Uj5o47zuiXf3lHX\/3qOd28+Zw+\/\/lT1aqdazInVG4T8zpMzuRwZttVIC9t9ff1heJ9zrcDXZopNW17XcpNellIWx0hXxcuFDPThC2hsjgpLyHpvZosCunxx2Vnx0oSVS8DybJe1cTg5Wcgz1X9jJzVORxKBwfnNBzWq5Zfs+X58vrIxAyYmMHyehJbXtcM5Dr6BXYvXSsunTOw0Em5c69iYJ2BPK\/34nle20ki\/mcnZVltx8+1yUCeq5qQi6Lu8s7OvnZ3r9VG+VkUtT\/PS2OJ\/5iYQdPEDO9iiV2NTa9QBjhL9l8WMtBAcemfgblPyraBsyHTPSSAx95UTPsAAuVuvPiiDj\/3uXovXRR1isrL1XJelUmc+18L1kG3Ppv8XXk\/zrdvtTQfbdr2upabFNfmD\/nauNdeO6wudpCpU6cOqwsdL764r5MnH9ONG\/vVcRY+wFk0ctkIbZ+23eIzwC2yr6E8d2m\/a7m2uL6+ULzP+XaXsRwnZtr2JpXjLJmJ2fq2p9tfFtJWR8jXhQvFWB\/WVc59Um5KzKI35qZ++DxvXgIu79v44ICEVcO3YeEAusG34eHAzn65s37sMe18cfTARJIcuVxNDKCMwbfh4QC6wbfh4QC6wbfh4QC6IWRH7vC21wgWhVQUddYuXnxTn\/vcYW2MPgcDVT8vxywKaThEi4gZWJ8MPKNnxp09dXhKDxw+MLaj0i8D39cvvH80k68dVfuyf23zL8Fbm4Dbkm\/jgwPoBte2BxCa3ixlZUxS9s2vfEU77363mJjhebq6Oo3KMswKxAGrn7fjYFdO54M3UuEDRhPnloOHAy7v28SF6vPfSkU5QLwBG5iNxAboBmxgNv2h7+TPOCQxAN2ADbApx4NAlDMOHmAD9ElxjJd66APxlAPoobJtnE3CSSJ96lMnqjej+blLUylJqF164YVarton22\/XPnFwNg+wPk6fPl11o0\/9Xcu1xfX1heJ9zrf7jGma2GnbayuXq\/7TaNkbnSX7\/WurI+Trwrkxo+bXXsx9Ul63DPHmJeD227fxwQF0g2\/DwwF0g29zGnXyq1\/ViSeftJDq7Pjw1VelJLnFlRplQamO\/\/k2DjiAbvBteDiAbsA23SQcMBsZ3+h1srwEfZJUjEGOwJgolT\/4B3fKTylJKsHq1rVr76+evq6Z+jNJalkUtVzVTw60DU0T9XA4LC\/N34i40T0H++VVMg7WOPhbJ3zm4DPjrypnyT944weXst6\/+c1vjvuxzsrcJ2U2WjbgdUkSZzzA7a9v44MD6AbfhocD6IYjdlGIJ4B2f\/Zna3eSVBOyskxH4mpvxfm8bxMKB9ANvg0PB9ANvg0PB9AN8Y1eu9X6sHwgyRFAN9x3Xz0p57lUFGJ162d\/drd6+pqfRBFXFFJRoElJUst5fbI9tmFSu2zTBuoJxT\/yyCPVAQsHKBH1wdukPLzzne8UB7qT4lbJ\/72T39Mru6+MvwKP7Cxvvd9\/\/\/3jfqyzMvdJ2ZLDxuvDfJsoOz2AUE68R357nCQ6+LVf0\/ULF8Sln7a8TKq\/yd+V9+N8u61vs\/BN2157uVs9mxTX5g\/52rg\/9+du3UN+8MHD8eS7s7OvJ57YrzqV51JRVKoefriW8\/q0CbVJsp2G2m7iQ7FJ+V3m4GTWoC1uTXDZsk\/dXcu1xfX1heJ9zrXvvPNOcaDbZ1x9Y932+pRtKvezu6OTiTIgUaJPH3y6vAJ0s7pF49evcmlad6XrtnJdODfmnnvuwVx7zH1SZkNu2vjXPnvTDqAoJCZk\/1WZztPV01Ydy61eBtJUKueoqmP7+zuVfP\/7r+n06Ycqna+CPXVNXJpW9NI+2F7ZbpfWgdjwWmSgUKFh+afRMtBATMyKy7EyMPdJ+Vi9W+PCHCW6DwqNh1IUqq5f+hMye+YyqLFc6XP\/TYpr8nfl\/TjfdvsyD33a9rqWmxTX5g\/5mjje4PXUUzvjs2DL1W\/9Vv1GrwcfPHXkvTD8J18W01UuKs4mayZsA9yi2qedUJ7hJ6Fruba4vr5QvM\/59qRxHNc\/bXuhcrnKyzu6tVzQhep2TnC\/V4aF6ijp6l\/I14ULxVQVrvHH3CdlNlo24DXO0ey6zsTLq5uKoq4zSeqXgaRpbcfPjcpAUaj6fbJ7\/GUDLApV95TNThKJn6IPBsaspmR7drGavYy9mncGOEu+pEvjZgaq\/xSXY2dg7pOyTchIH8fu\/bpUUBQSE7K7d2YPHC9Xr8sa7N3PolB1QSTP66JJUh9\/vf12Pfl+9KP75VnFNaVp\/RtlvgppWsfGTzIQscoZ4CyZidn6uKfbXxZivij7ZWDuk7J7VO3r\/bq6XtHjB39ee03V3tmdkLlGySQdGNK43OGth4MCYeVZ1oHiG71uz8y880eLoTZcjlXrXhCx4693vateZ\/xumQe8eKPXF76wX03M1Lsq4OCZbXVV+uP2w82zy0\/Su5Zri+vrC8X7nG9PGsdx\/dO255dzXxaSKNGg\/FO5+HElNf7X1xeK9znfHje2xsrcJ+U1zs2xu37X009XLwNRUdR1DQb16dJgUNvxc6MywENcTMh2\/MUrNUNv8Fr2oJl027CqE\/Ky8zaL9jehDs6SgY1lL54lWypmIucyKdsGTw9ND0n8G4mi0O7ly7q7nJSr8SVJfc3y6lUpSSqq6aPrgwuT4pr8XXk\/zreb+j8rftr2upabFNfmD\/lef31Xv\/7rp\/XFL9ZPVyeJ9IUv7OgrX7m7+nkIeQmVg180mHTbsOj+9Glv2hx2LdcW19cXivc53+6Ti2lip23PLefeS+YsOVUqW9w440z29YXifc63ra11lnOZlG2DJzGmhyT+jUM5Id92uZrJOMs2bqhxQHUGmlZ5mtb++BkzsCkZ4D6ye5Y80ECS4jLDDMxlUp5h\/9anKvbMWaYjLwOxm4lJsj7jiD3tlYE8P7rK01TlvX4pSXpVs\/Bg\/8rVwjsQG1zLDPhnyXvx0vXM1+NcJmU2+Jn3dEEV8iYt4Dbn2\/jgAHp1z5gbiQAiSXTw6U9Xb+biQYRxHL4Svl1S1Ru84Im3B7iw8bmAA5Pi3nrrLb366qtVvVbeLwcPB\/z6XJu4pvrwGagHmI3EBugGbGA2EhugA9qn\/0hsAzHAbCQ2QCee\/CGNgwfYAB1\/WxzjxW\/xSOCW5YGt8+dhanD\/+Od+7mCcc78Ns62eutTiP9k+\/StXcIvvSf8Wp81h13JtcX19oXif8+3+GelXYtr2KPfijRc16WUhxLnbjZylry8U73O+7TR3bHVZFcxlUl7WYGbRrr0I3q0LzrXR4YDyXNXl6uEQWkoSyblcTQyonfWnb8PCAXSDb8PDAXSDb8Pby+3RDaE4OGAxSN+GC9XHy\/PxGSgHzEZiA3QDNjAbiQ3QDdTvc9jAYpDYAN2ADcxGYgN0AzYwG4nNeNENcACbB7qefvquI\/ePP\/3pA128+J3xi\/iJM1AOmI1kbMiImIF1yYD7jmv6fEEXEBEzzkCclL2E2sveXRrOtdFPfu974n92qiZkLl1DOper7QEE\/ps5v7xvUxQOWDne7YuNzwUcmBSXJIne8573jB8yog6\/nHE+79vEherj5fn4DJQDZiOxAboBG5iNxAbogPHRf\/KHbSAGmI3EBuiU441ClDMOHmAD9ElxjJd6WA\/EUw6wqt\/3vh09\/fTd0ONjsMuXd8s7F2fktuu3YTZ13nXXXVX5+NEvA24O+5TsWq4trq8vFO9zvt1nTNPETtve67uv61+c+BfjJgeq\/+QtbfX39YXifc63ve6spTm3SZnLYffee6\/a5CpmjJUM3L75tso98+4nPqHdwP\/s5JZDpyxAN\/g2PBxAN\/g2PBxAN\/g2PBxAN\/g2PBxAN\/g2PBxAN\/DyfNOR+AG6ARuYjcQG6AZsYDYSG6AbsIHZSGyAbsAGZiOxAboBG5iNxAboBmyesOZydbn6K7o87rnt\/jFxoAoYfWCDkVkJ367I+BEzsKIZ8B\/wiveS57ei5jYp+\/esQvb8hjXHmvNc5SmRlOd1I+yZncvVNRk\/Ny0DWVavdpuQ7aLIpo0zjidmIJQB\/wGvgQaKy3wyMLdJeT7dXWKt7I2zTDp\/\/lYnbM+cJLe4kTbtAwhdy02Ka\/J35f043x4Nc25i2va6lpsUZ\/7XXjtUlunIfxrB\/eMLF66PH+giCRbvPsTlc75NuWXBv4JFP0Ic\/Cph2hx2LdcW19cXivc53553rqdpL1f9p9Gyp+ZXarbV39cXivc53x51ca1FnJS7rD4mZK5bOk9XV\/97QJZ1KR1j1jQDPNDF\/WN3tW\/CRZHQVasQt6arLXZ7xhnwz5Ldl4XMuKlYXZmBuUzKbOBl3ZvxL8uOXrdMEqnDnpl7hu6DQl2T0bXcpLgmfyO\/u1telT8zfjDMj\/PtruOZNm7a9rqWmxTH\/eOPfOQ0jw9UQ0iSerUniRQq24ULxVSVx4\/OGZg2h13LtcX19YXifc63OydiysC+7fn3kiddtm6rv68vFO9zvj1lWlaq2Fwm5ZUa4bSd4ew4yzS+bkk9LZercUdsRgay7Ohx2KaudrtcbWvNt42Pcnsz4J8l77Vcut7eLM125HFSDuWTCTlerg5l5gh3bffaEXsFjV5dYrVnmcbHYUlSv7I8y3pVsxbBTMBc0QLWYXSAz7gotzcDnCVPelnI9mZnfiOPk7Kf2yw7epqUJPV1yyzzI1vtaR9A6FpuUlyTvyvvx\/k2E\/JjJx\/TSzsvteZhWqffXtd6upbz45iQuXcMaIv\/4elv\/s03xFu7sF34ZfF14UIxlF00fA44HwAAEABJREFUmHSZfJvaxUdMk3+Z\/LQ57FquLa6vLxTvc74979z2aY\/Hu9z+\/NDNHzrycKPrM72t\/r6+ULzP+bb1Y53l9kzK7HXzXMpzKc916rXXdPbgQCr1MbLs9tOk69elJFnndTzzvnME\/bger+o1WRlr+sFXgwsjw2E9gCSRnn\/+UB\/84Fs1sWqfsT8xA3POANv4M3pm3MrZg7M6d3BubEdlfhmYy6TMUXYXzG9YgZptz8vet8SpD39Y37h5U8jqrVwld+S65ZUrUpYFKupGTfsAQtdyk+Ka\/F15P87s\/Z19MRGz0ZIJ5Hk5L4CGnAGsPd5+1ae6ruUs7qWXdnTmjMTXg3bSVNULQX7gB44++IbPYGXdvnXhQjFWZ5TdMjBtDruWa4vr6wvF+5xvd8vC9FFd22O7zpXLlsu7l8vt5NaDoMb7sq3+vr5QvM\/5tt+fdbTnMilz+csFiQnZ8AvHYFDdKNz\/6Ef19F13aX9n52gXBoP6cvVgcJSPVpUB7jGxsZ46PKXnbj5XcdiZpj+A0ZKWLFN1PGbN80AXx2JmR7nQDMTGVigD\/gNeA8X9oRa0zGVSdvvOGTMTssthw7vcwvSHH5ayTLpwoWry1OFhJccf7JmTZGxG5VYGmHzZWBMl+oXDX6guZz1\/+LxY4PGjrzo4K84yHbkwwmrPslXv+Wz6N2n7Y9skZjatxVrWLQNsx8D6vRefuLZULETOfVJeyCimaOTUt76li9\/5Tl0ySaQ0rfUZfU77AELXcpPimvxdeT\/utcPX9OHDD1fZuaALGpR\/KhekbbRc1i6pmfzz2+9a6aRyTMj+nYrQz87b6gn5unChmK7jmnUcky6Tbwj4Zt3erOqbNoddy7XF9fVV8devH3k4yud8e1Z5aqqnS3scYFv5RIl4WUiXciqXtri+vlC8z\/l22YW1\/7e1k\/K1979f13Z3pTStbyLec8\/ar8x5DYD7Sx\/b+Vh5qX9fg\/LPv1SNDU\/cPO4va0aLTchIqjx79kCvvnqoJMHaPjD5hrB9mYgjtgywDbtnyWzX5otyMRmY+6TMRr+qR+M\/eeKE9r\/whblketoHELqWmxTX5O\/Ku3FDDcWGmijRlfJP5cIDX0+eeFJP7TxVWirZK5UkLtPxrwO77avH0lQuyyT3gS4uV1+7tiv3YS23maZ6iAn5unChGOpbFNgO+7TVN75P3dPGTpvDruXa4vr6QvE+59vT5qVruUntsa1bXWzvXBXD9spBBdEW19cXivc53w52as3IuU\/K5IOJ2Qf8snHbQ17L7tAKts+Rs13OssvUTLrv3nm3fvX7f1X47tAd4gzZjqr5KQXlVmE4nBVnmbb2\/rG7DtgGmWiBy\/s6fkC874v25maAbZbt2UbIpGzbtHFRzj8DC5mU2cANDAkduRS88IKUZTr1ta\/p4ptv6tQv\/qJU2oJfSodWu1E2zD3tVZ28pEti8kVWhPPBBs1ZMhQ6R9xI7GWBCbnL\/eNl9W8Z7TLRArbBJuAHy+hfbHN5GbDt13pg273ZGy9XZIBzn5TZ8NnAgY0ZHd7seUnaMIzbyHNx2sRkfPE736l0bOX5OGQWyrQPIHQtNymuyd+Vd+N40IPJmUkWkB9+EvUbN35Drx6+KvzyFibuM6r\/MmWycuq4uO13LFKFWTn+y0UmZCZmHEki8UBXkmBJFuf+V4u1p\/5s84d8XbhQTN3a4j\/ZBpuw+N50b3HaHHYt1xbX1xeK9znf7p6J6SKb2mP75AqX1TpQ\/afR0lRu5B6Ltri+vlC8z\/n2uCNrrMx9Ul5WbpiMbadT9WEwqPfK7JlLXLt8WY+ePKn9Z589witJqvD4UWeAjZWnqpE1I7G5MhkzMTNZ36P6ITn0Pe0JqdFCOX+C9o\/IR6EzE08\/fZfe\/e4d2YTM\/ePr8cVsM8tvrGjzMsB26m6XbMebN8r1GNH3zbabs62NiTVUIzwI+YxjQjZ9LAeDckYZ4eGH9QpPXw8Gt7jBQEqScfhxlGkfQOhablJck78rT9wzZ54R945tY7XJFhuQHySTLj4eAuOs+Lrqvz2FJ2gugXMOTSzlqccH7Z85M\/kNQm45JuHLl3f19NN3V3SSSEzI3J2oCOdjUv1t\/pCvCxeKcboU1Q4ZmDaHXcu1xfX1heJ9zrc7pOBYIU3tsQ1bxWzLHHibjWwqh89FW1xfXyje53zb7cu66is5KTPhglBS4ZlwAToxSBdwBnhizY5ycgY4ambCtA2VjXSvnGCZapHU8IndT4gnsPmpFPYFXZC7IVOGOigDKOf6aYP6u0zQ6rAwIfOfSQDCk6R5QsYfETMQM1BngG3RPThmW6098XMZGZj7pMyEyMQIGCASwGOHgA\/4Pr8cMcahG6yc+cyeVm5TuaGG1Z1gJkyVC5OrnQGXpphomVzZkB86\/ZCQ2PBqWKgD\/1Vd1fXyj42eMhZOHbQ37QRdFKpelzkc1jUmSX1HIk1rO362Z4DtxEdbCTe2LS761iMD3J6ynrKthp4RMX+U88\/A3CdlG4I\/YRo\/L8mOg7qRAD0EHvSZB9566y29+uqr1QNFfervWm5SXJO\/ieeNXZ9885NyN9C\/+bt\/U7958Jt64PABueV4xaZGC\/eVsbuOkfjPHX5OvJ6T+9KfPfxsVf+oOrkT9OnD02Iy50Chqf7hUEd+f\/ypT72lf\/WvXtW73nWgpjLw7niwfbT5Q74unBujFVnYNmzbdCV8qIvwXeJCZWfBTftgT9dybXF9faF4n\/PtWeSorQ6\/Pba3XLlsGWig0OKXC8XAtcX19YXifc636cO6Y+6Tsm3EJMp0Nmp0uHmBNlw0tTMcDnXjxo2ZY39\/X6+\/\/nrveruWmxTX5L\/F3xrzyzdfFmepT9\/9dJUmJs6Pv\/5xfeg3PzTuv1vuD9z4A3p2\/1m96+Bd+tIbXxL2NDmk3BM3ntDXb3xdv3HjN3TxzYv64FsfrPrAx\/7OvjiD5kCBe9A\/cvAj+tR3PyVr61Of+q4ef5xI6dSpQ128+KaeeKJb3vcnrJ82f8jXhXNjvvnNb9YdX6NPtlm2qS5d9g9yon3YepC4rPwMNZQtbPc\/evijK9nPLvmxcay7nPukvOoJeuSRR3Ty5MmZI0kSvec979Hp06d71d213KS4Jr\/Lf+\/k9\/TVk18Vl6H3ywmQdcWR8lVd1d+\/6+8f6b9bjnx94H\/5gL78u1\/WX\/xf\/6Kwj4v3nHiPLt95Wb+888viDHpP9R99AvSPl5Vw4HDvmXv16MlH9fRd9UFEmepyMt7R5ct3yu9nU78mxbX5Q74unBtz\/\/33M6yNxXA4n4PdN954Q7yF7ebNm+ODMztIa5Ndy7XF9fWF4n3OtfcnHCi2ja+rz22Pg3EOeu1L+EcO\/4j+6I0\/qlBdbrmQ37i2uL6+ULzPufY6Huha7l259ZMyO8rd3V0Z2OCB2UjfNs7nfXuV417ffV1\/fvfP62d3f7b6PiRKqmmQCRmdvrtgbMC4O++8U2f\/59lx3uDxA3QDNjAbiQ3QDdgA+76d+5SVf1zm\/r3D3xN92tOe3OWV3Vek7JL09h3S9TNSllV9oTygLoBuwAZmI7EBugEbmI3EBugGbGA2EhugG7CB2ch7VuRd65z5cgbs5hYb3uXQfS4Ud1hesbj53E3N62C36SBrE\/h3vvOdOnHihBY1lv9w4j+wWsfglbmLanse7WzKge5aTcrsFNgR2LcIHc7sWUg74nPrgnNtdDiAbvBteDiAbvBteDiAbvBteDiAbvBteDiAbsAuVIgJj8vB6CoXLlsx8cGXZnWkjHRBWeByr7\/+umtW5fwYbOAGYoOuHGfv9O168bZ05rqU7blFxTg44r9Dd4hxEUv9QM6CDRxqbn2e1IafOzd+kTrbEO0hDSEbzgWxoW3va6e+poNzB\/r8A58\/cpDEgUjErYP\/UC440H3HO96xsLz9w91\/OF6lbGOP7DyysLZD4z8utyoHuuOkTqnMfVJmw2UDBuj009Wx+4A6KA\/Q+5TtEmtHcG4snGujwwF0g2sfHBzo+vXr4yNfi0G6cdgADli5w8PD6ogZnwtiwKS4u+66S9\/73veq+0NW\/nsnvyfuzzKBwXFGvKe96nIx94fpr7Xr1k9sqD6O6vEZ6BcwG4kN0A3YwGwkNkAHtP+9sv\/WRnkirDPlCbGKRLqUaS97W\/Ykt5ylUCHGxyVuLsszQWu0UD\/ApH7GS\/3GwQNsxouffBgHjx4q24WzGOqkfurqhvlFsQ1Ngt9607Y3VP2ncsmVy829ZrS4OexTZddybXF9faF4n\/PtPmOaJtba+9bht8Q60mjZ09ED3RE9FlaO7+6YDChtcX19oXif8+1Al9aOmvuk7GaEjRm4XJvOziLkhwch33E5O1pz64FzbXQ4gG7wbXg4gG7wbXg4gG7wbXg4gG7wbXg76kYH7CDv371ftiEyIV\/VVcH75X1b5eLXV1LiqB5poBwwG4kN0A3YwGwkNkA3UP\/rr+8qy6TQ748ZA\/1\/W28rNEHv7+xXE7SdQV\/evazXd4+e3dMmsDaR2IwX3QAHzEZiA3QDNjAbiQ3QDb5t\/KIl22MTQn0hNrTtFSrEAZ+chYMj+745dFRXIAOsG+sG2xFnymZHudwMzH1Sto2YDdkF\/HKHHm6dI0Hgen0bHxxAN7g2O13eSIXP5UO2ccRZOe5BYuNzAQcmxVGeh8yYmJi4bCNkA9zTnpjE0KkL+PW5Nu1bfUhsQDmkARuYjcQG6AZsYDYSG6AD2r\/jjjN63\/t2jkzIzz9\/qM997pCQMSjHJXjGaRM0P7VifBbEpEEOuLwNmKDvOHNH5aZ8pYw+sBkn6w8JDQfQ6Rs+9D6clbM6Kb9suNuk6fQJHenCtlmkC3J7XufHobvXdse6P1GPHZKm0Xd3d8srJv3e9EY7Xcu1xfX1heJ9zrfp6zxBe3zvX9p5adwML\/4ZGw0K5fjOT\/rutsX19YXifc63G7q\/VvTcJ+W1ykbZWe45glId\/\/NtHHAA3eDb8HAA3eDb8HAA3eDb8HAA3eDb8HBfefMrYgJiMoJj4rqiK2Ly0mghDozMSvg2JBxAN7z++tGzTvzA\/EhsgG7ABmYjsQE6KArpwQcPVRRYUpLULwT5A3+g\/ilXzdaflAO1JTEZ81Orf339X4uDDw5C4DRamETICbkB\/MzKPaOjLjAKrwQ2qIzRBzYYmZXABpUx+sAGI7MSfu4qckU+mJCZdP3uwIdALskp8T9w7Qd08rGTum\/\/PkzBuxO24rL0DLC+rBNsF+7+wPgol5eBOCl7uee+IXBp38YHB9ANvg0PB9ANvg0PB9ANvg0PB9ANvs2OkJ86PXniSQsRExM\/NfIvU1EWjANLxbdLqrq\/7fPcj8VnwA\/MRmIDdAM2MBuJDdCzTOXZkLS\/v4NZvb+6vD2vJFGwH5QDVfDoAxvYTofJmfGTB7hRmLiSwM+smDiYoNlB\/aeT\/6lqx2KQ1AXQDdjAbCQ2QDdgA7ORd911F2LtUd9FHlbjIK8c9IQOBiIAABAASURBVGFc1VXZdy1XLvKqrVhWe5DsG1hn1ktbR2ZHufwMzH1S5sjaP+rGhl\/+8G\/vAZdDgOvxbXxwAN3g2vYAApd7XJ5Y3zYO3spxWRQbnws40BTHRscE4\/\/UiZ2iX4564YBfn2sTB\/htKP1CB9zzRRqoB5iNxAboBmxgNhLbv3\/MC0H8\/1CCOEAZAzbAtn6H8n7fTv1TKyZo3lTGy0rcnRK54yyCn4px\/503nPFADPVSP0APtdGFsxhySF0HZw+obm1BvuzytDsh24A4ADKdvDI5mz2tdHPYp46u5dri+vpC8T7n233GNE2sOyGzzi7oQqdquvazLa6vLxTvc77daTArHjT3SZkJmBwgDSEbLmL6DLCDZOLlbA+dms4enBWvs4THXmUUhap7x\/ZAFxPyl770xm33j2c1BnZIF79zscoPk\/Se9uRP0JxBv2\/nfdUtAHI4i0lFo+Xg3IFufuOmXjv12ohZjrBt0peTDpr5jnHwZ71m5275e+uDb+nDpz4sbDtzVrnYBF6q8d8SMsA64+DImmYbYB2Z3Sajb3EZmPukzMY9CYsb7uSWOIsBbqRv44MD6AbX5kyIByPwuXzINo44K8eZHjY+F3DAjWNjYwdpGxwbGw868TAHZ4dW3i8HDwfc+nybOPrDg2NIbEAc0oANzEZiA3QDNjC7KHTbfyjx4os7+qt\/9YSFjCXlwJgoFWxQqrJxoBuHDrABuhvHvXYm3au6KiZocsf7vokD5JfckmMOenhI7IUzL+Aa\/+TMrS\/UBpzFcMncJiiTVWVL+GjaNid1hTMu8kIc3zfyh3546lBvfPkNcQBDrvABlQvxFleaU\/2zHLrfwy4VdS3XFtfXF4r3Od\/uMpZpY1gnblkORF27Te\/az7a4vr5QvM\/5dtsY1sU390l5XRJh\/eSBHGA20reN83nfnnec7eTY+aHTHjtAJhcedPL749vEwwF0g2\/DwwF0g\/+wEn5gfiQ2QDdgA+zhUNX946LA0pH7x8SA2lN\/YoPaqj+xQW3Vn9igtupPbFBb9Sc2qC2J\/JE73sfNBM2Oy\/1fc8jzJV0Skyl5Z\/J2JxrqAnIWbABFecoizWayR18nMGZyRZ8ZS65cLJ8\/9XlEBXjGhoQgPlOmuCwnA8\/omXHDrIt4ljxOx0opC5mU\/UtjZq9UJkad4YEcMDIr4duQcADd4NvwcADd4NvwcADd4NvwcIAdHTt3Jgh4NjImECYSdGIAPoNvw8MBdINvw8MBdMNxH\/TKMo3\/Q4kkqSfkLLPa+z\/UdavkbMomSsQkckVXZA+JkWONFs54ed0n6+AO3SEmIB6w8\/OEDSjGGSYT2M7+jk4+ehJK2LSjBS22\/U2Sk7pDXiyGsTOOj1\/7uHav7Rp9RLq5O+KIxtwzwLoB1lBcF5aJ1ZNzn5TZ8O0SGcNHdyX6KmF3d7e6BOr2Cc610eEAusG1D0Zv9NrZ2elVn5U7PDy8rZzKhTa4fMpZmm1kTB5XVb8IpAwZ\/7t58+b48iokZd36jfN53yYO+PVN+6AXD3Rdvrxb3UOm3iSRLlyQsgyrBv2kPfJXM\/UnfQO1VX9iAyzK8SYuyhkHD7AB+qQ4YqiH9YDObQAmTvA\/D\/6nnrv5nLjMjc\/A+uABuzt37xTrh3d3U4f1BT8TOOvrqf2ntPvKrvjftigPjx993mAbdEF7IRu+DQMNtFf+abRwkHiP7ql+EvXR\/Y+O2FsCP\/m7xfTXbL3ZeulaQ9dybXF9faF4n\/PtruPpG8f3y8pwu+bHDn\/MzE6yaz\/b4vr6QvE+59udBrPiQcedlKcaHjsAJuupCm9xIc6O2anZBsbGxU7Rzo7XITVFoer+sT3QlSTS1atSlq1D72\/18dzBObEu3tbbIv\/+BF2o0Bd3vihe9\/nunXeLWCYllQsPRT2sh0tNcic282uBC9sh26PbJDa8yzXpjIsx4GfMPOBlOtIH310OWBZ1AOK3v40268XN9we\/e+u\/R93GfKz6mJcyKa96UmbRP87IeNCLs6Q+9TWV47InOzN2atTH2dazO8+KnaICS1M9XXk\/zrcDTU6k8lxH7h8PBpL9\/tgvPG17XctNimvz+z7WxVM7T8km6D3tHRkOO0XWG9Icv6PfMVVZ+cfEhp\/LwFqz5aquihyoXHianAe9vnXqW6UlwYfywTgZt3oufu67Fu9ari2ury8U73O+3XU8feL47lk86+Mrd39Fs9ovWb0m28bT1xeK9znftn6ss4yTsrf2uCx2eHj7axy9sOqyMLEu79v44AC6wbfh4QC6AZsdNTsv9yyKszLOztiRE0sc0gUc8DnXRicGoBt8Gx4OoBtCdhOXZarOkK3sZz97WJ0hm005YDYSG6AbsIHZSGyAbsAGZiOxAboBG5iNxAboBmxgNhIboBu4esFl67f1tlhH\/oREHDvJHzn3I7rxGzd0eKr+rtn9Wc5oWN\/ErROs\/\/SZn0QhmQDgGQ+5sO8rPkAemJwLFYrLfDJAbjmgt9r9dWB8lKuTgblPyu6lMNO5NIa+Omm41ROekgW3GFX\/vZ9roxMD0A2+DQ8H0A2+DQ8H0A0v33xZ7LTYecGxk+OFFzwdjG3wy8HDAXSDb8PDAXSDb8PDAXRDl6evX375pj71qe8euX988eKbeuKJG1ZNJakbVMboAxuMzEpgg8oYfWCDkVkJbFAZow9sMDIrgQ0qY\/SBDUZmJbBBZYw+sMHIrAQ2wGB9uU9us0OEw2fY2a\/fWgbPBAbPU7LsTNHnDbZDtkcf8H3aZmz+AQiX6eGph\/Fd1VX5MRyEcAWIiZu4iNlmgPy6NbJOXDvqq5eBuU\/KDNndwNEB\/Cri5MmTAm7ffBsfHEA3uLY9gMATyi5PrG8bB0+5F2+8KM62+K8HbefMzo2d2uU7L1f9I46HiA7Ls3rKUYcLXuPIf32If8TfVg6essCvz7WJC9XH2PAZqAeYXRTSY4+d1NNP311RSVLfP758+c6qLxU5+qAcGJmiffo\/qQ3iKQfQKUdeKGccPMAG6JPiGC\/1WP4oB5rKhuqDe\/v62\/r0wacpJtYl7+PmSW4ehNp95dZTyvjs4IsdJ5OYFrSwPfqYpmkmVr6nlOVgBBvdBRxnzf74GDsHoG5sSCen7noJxYS4ruXa4vr6QvE+59uhvk\/Lud8p6mDd8CxEzB\/ZWF3MZVL2j7qb7FVMC\/cogNs338YHB9ANvg0PB9ANvg0PB\/iZzWMnH6seEoJn58XZxVXV9+2IAfgMvg3Pfz0YejoanwvKAp9zbfRQfaH6ra6iUHW52t5fnST1hJwkqp4qtzjqBtgA3UD9PocNLAaJDdAN2MBsJDZAN2ADs5HYjBfdAAfMRmIDdAM2MBvJ\/WZ2iOwkHz35qFinT+w\/oRNPnhAL\/ON6XEjimLi0wMXdPmkWGzkNfubaz2j32q6+sP+FxuKMn4mZ77UbxFkdZ81Il4\/6dBkgj3ynrLSfb+OjXK0MzGVS9o+6m+zVSsVse7O7u6s+D3qx8bAz5kldJmZ6w87rajkZw2O7mFR\/k78r78f5ttsXX88ylWOXiqL28P7q69elJKntLp+d2\/Mq61puUlybP+SbxNmlaX7XzCTtdpt7fuxAWd8W5\/rnqTMB2\/Zp7WDDm91X8r9EdSnD9\/pq+f12Y9kOOGPmIMXlTQ\/l2Xxtsmu5tri+vlC8z\/l22xj6+rgNYmX4bnHAN217Xcu1xfX1heJ9zrdtvOss5zIpd9mgu8TMI7FFIQ2H0le\/ek43bz6nz3\/+lLJsHi11r5MdETshLuFRig1oT3vibAIdblVQFNJwKH3iE7t69NGT+tjHdpRlde+KQsoyHbl\/zIScZbV\/mz9Zj1d1tUoB6\/kFvVDpuXJdKv\/wMyEjtUULEwXfc6Q7bA5U4lmzm5F+Ot8rYKX2dPQXAcZHuXoZmMuk7A9zWROw3488V3VJ9fHHpeFQ5X3Lc5XkN7Nnzkh5LnEPEbhlfRsfHEA3+DY8HEA3uHamTO7Oh50yO294N06jBQ6MzEr4NiQcQDf4NjwcQDf4Njxcnh\/N3yvlPdHhUNUkTP7IK7kkPknqy9Wf+9xhlVM4A3UBs5HYAN2ADcxGYgN0AzYwG4kN0A3YwGwkNkA3YAOzkdgA3YANzEZiA3QDNhhoUO4a98Ty8+d+vnry2l5LeUEXhF\/1slWf9p3noMQdeKFC58u\/TJni0i8DHOhZCfKbKlVc1iMDC5mUVyEVea5qQi6Kujc7O\/vlvc1rtVF+FkXt\/8f\/+I3bnra2p2nLsPE\/ODAmSsW17QEOHqpw+TKsqp8dTqZMtvHwUxqerOa\/FOThIHbifjmVCxyw+pviXnvtNb388stHJkS\/XFld1Ref923iyMv581JRYEn8L05\/8k\/+f7VRfhaFlOelUv5LknpCTpL6yXXqK+nxP2wwJkoFG5Rq9Y\/x0X\/yVxGjD2LAyKwENsCgHGWAcfAAG6BPiiN\/1EF+iaccQA+V7cIRc+H6BfG+bNY\/P4lCDjRQVv5pyxcmjtBZM9uIHbiSQ3e9dE1Z13JtcX19oXif8+2u42mL4zuVK5ctF8oDPtOnba9ruba4vr5QvM\/5to1zneVWTMpFIXEWx4pKElX\/6cGLL+7r5MnHyklp\/8jvZX\/yJ0+U\/ElCx7CnbsdEqcCBUh3\/820ccADd8C9P\/kuxk2FnA5coES8C4clqOYtfDhccQDf4Nvw73\/lOnThRP0iEDUJxcAC\/wbeLQiIv+JOkzt+rrx7qG9+4qV\/4hfp3tvgM1537x9QFzIfEBugGbGA2kv77HDbAb8AGZiOxAboBG5iNxAboBmzyZzYSDqAbsIHZSGyAbsAGZrtnhKx717aYRUm7f2xXs5AAfuZ96FAh+biqq9or\/9xwJhrOmv378W5M1OsM2H4Fi3xmyhSX9cnA1kzKRVGvlNA9zsFAunKl9vO08K\/8yq2fqsDyMAHSBRzwObPx8aDX93\/\/95dn5LsVzY6FDeQTu5+obD72yp2PnR1QBlBuZ2dnXE7Ogt\/QFnf33XfrB37gB6o39xSFVBQS75wGd9xxRoyzKG5xLo8O3Dji6caFC1Kaqir\/8z\/\/zuqeMryL4fCWZX29xXR7+ppy9J\/89SlLOfJCOfSmsvja4sgfftYDdRAPTMfntoFvEmcxvEf72f1ntbO\/IybkROWRjpa3MAH7WF5v6pbZTmy7qJn6k1eX\/viZH5c9DFmzkz8t9zvldtUW3RbX1xeK9znfbutbFx\/7GO7HW+xAA7nLtO11LdcW19cXivc533bHuq76VkzK7n3ONL21qnjQ68MfPlVd1rYYvOjcI50pzhc6U5zX+Ci2KHfE2Z4u3ZHpjjs0Vxx3HFy21mix3Lz73TtHfn\/MwU6S1EEvvFDL+Nmcgfv276t+EoVsjlqMhzNjAy2iI5cNDlau6qr2yj+3L7lycaWJiVvLX1aqB+TG7dAFlUfRLhH1lc\/A3CZlNmwDWTDdJNyikSRHWzw4OKc8l\/JcKopbvqKQikIqCqkopKKyzGsDAAAQAElEQVSQikIqCqkopKKQikIqCqkopKKQikIqCqkopKKQikIqCokj1+JCpuLqGSkpCZVLPpDOX5UuZaWx3v+SRLL\/UCJJ6rEURS3jZ3sGdl+pr6C0R83Xy\/ZoZ8nWEja82cuWmTJx1swkLWe5pEvikrZDbbXKvoacWBIGqv8Ul7XKwFwmZTbqLlhUppKkbinPa2mfu7vXNBhIg4E0GBhbzp2JlKZSmkppKqWplKZSmkppKqWplKZSmkppKqWplKZSmkppKqWp9KM\/eqizH\/qtcsY6L2WXqsoTJUqLPaXPXNXehaS6t80ZpgveCc1rKJEu7+v42+JcP5fmDT\/3cwf68pffqO4FG4f0edceDDReuIfsxnNvOUmkopCKog5Lkloe53PaBzi6lpsU1+YP+bpwoZjj5GibyrLtMDHz3nd33JwZctaMdHlf75r7tri+vlC8z\/m23+8+NjlgYrYye9ozdSw7tzcuUStdy7XF9fWF4n3Ot+vervfnXCblVUvJww\/f6lGe39J50OvZZ+sHvS5cuMUzATLxTIu9K4Xu++Wn9Mqv3C8lRVUxO5WruqorSVbdv84yKcukLJOyTMoyKcskfkJ08eJ3Kj3LpCyTskzKMinLpCyTsmxynFtPmkppKqWplKbSBz\/4ViXTVEpTKU2lNL2dtzg3N\/fdt1PFpmkdXw2u\/MhzqShKpfzn5rs047+YgZllgNfPPnfzuSP1MRFxxsxv\/Y84tsx4Rs+MR5wo0aD8U1zWLgNbMSmnaTk3JvW64f5oltU6n0UhZZnGT2fzU58f+7GjTxTbz2KIN8ABs5HY7CDYOdhlJH7qxNH9dV1XUv6pXIgrxZF\/cMB9cAH7SFBpwIFJcTzQcvr0aSHLYtU\/vxwkHPDrc+00DeePckUhZZmO5C9NqbkGMaC26k9sUFv1JzaorfphMPpvtkligNlIbIBu\/UY3Dh1gA\/RJceSNB7eQxFMOoIfKduEsxuqkroh+GSCHvKqU7WmggdxlqKGazpop565PNSxtcX19oXif8+2Gbk2kOUsGFrin28+S8U3bXtdybXF9faF4n\/NtxtiCtXBtxaTMmuCsFwl4WOnRR89Vb\/R68MFTwoZPEulLX5r+d8qf+u6njuwUEiXiqH5V\/lcn+52tnAUOOJR8G9\/f\/bs3ERXI1\/337+qv\/bX\/XTzwhY2DAxryh26gLmA2EhugG7CB2UhsgG7ABmYjsQG6ARuYjcQG6AZsYDYSG6AbsIHZSGyAbsAGZiOxAbrB\/x+2jF+05DYT948BbSMBPPaqgm3rquo\/t48cFHPWnKk8UtT2LHYSwIjJTapUcVnPDGzNpDwYSDyQZKupKOo3epmdJBKXrf\/KX7lb7m9K8fu2ccazI2An8PTdT+MSG8We9nS9\/HvPiff0rk+jxeofmZWAA5Ux+vBtaDiAbvBteDiAbvBtePLy\/PO3riAUhfTv\/\/3\/hqtCkkhf+MKOiKuI0Qd1gZFZCWxQGaMPbDAyK4ENKmP0gQ1GZiWwQWWMPrDByKwENqiM0Qc2GJmVwAaVMfrABiOzEtigMkYf2GBkVgIbVMbo46677hppURwnAwMNyq3sevk5kLswSTWdNbtxm6Cz\/3HPksnJJoxrW8fQe1Je50QNBtLbb9eT70c\/ui8e9EpTVfd4eeFFmqrkdiu44+QSiWujwwG7ZMZOAJ4J+Yqu6DMHnxFvHuJSJXH4DL4NDwfcBxew8bmAA5PiKHPz5k3ZJVdsv5xxPu\/bxIE\/+kdv6Pd+77A6eOFBsrNnD5Smt+ePWAN1AbOR2ADdgA3MZnz0n\/wZhyQGoBuwATbl2vLeNY66qMfyRzkAH2qjC2cx1Gl1Ud8yYWfFnBm7gF9mv5radnNoMYkSXS3\/9rQnd2GysrPmUDk31vS2uL6+ULzP+bb1o4+0fQ9lyIWfB3jDtO11LdcW19cXivc537ZxrrPcqknZVlSWSU88sV+ewT5Wnt3tK03N012ywXN2zP1jK8XGwNnxQANt8pJl9YNm9kavNN3k0caxrUsG2B5D2x+T1v2796vvC0fWYdzshzgxsL4ONFBS\/ikua5uBrZyUb19b\/Rg2BI7A2dgpyUawVx6ls1PABpwJdXmwhFgXXctNimvyd+X9ON92+zwPfdr2upabFNfmD\/m6cKGYeeRuk+uclEO2xaaz5odOP6RJr+lsq7+vLxTvc77dd925l60pe0EXEI2Ytr2u5dri+vpC8T7n240DXyNHnJR7rCwmYybeMzojdJULR6bsBOBLM\/6LGVirDHDJ2r9UjQ3fNhBi2vzL9rE9ctbMJO32hQNpDqhdbl119kGMx\/rPvgiYHeV6ZiBOyh3XGxsAG7NtBGzsnB0zIaN3rGahYbGxmIFJGbDJFWmgjOlIbAM2MHuVJdslEzPbqdtPzi45sEa6\/Lrp9J\/9kvXbH6fxUa5XBrZ6Uj44ezBxbfGl56ibjRidAmzsTMbw2CFM+wBC13KT4pr8XXk\/zrdDY54lN217XctNimvzh3xduFDMLHM2TV2cEU+CW6\/Futwi9WlyyHb6\/OHzR7rJtsxB9uN6\/AjfVn9fXyje53z7SGcmGNO8LGTa9rqWa4vr6wvF+5xvT0jZWri3dlI+OHegm9+4qddOvda4othw2Wj9s2OOvpmYGwtGxxwzEKte9Qz87u\/+rn7yJ39SP\/7jP94KnkKfBL+Oj3\/841XdH\/vYx26ru62ucwfndLW4qgcOH5CeKDM4wvCJoc48cUY\/\/OM\/XNXn199WJz7rn18O3ueI74LQ2KjPxV9+4i+LM2WNFl5QFKrbLYPu9wnOECof4iw+JN36Q2V9zupwyzVxflnfHqVi7cVWTso22bL2mHSRBlvRHF2f0RnZFz9RoqvlH7zKxeJKdfwPzgx7AAHb5UO2ccRZOX4KhI3PBRyYFEd53oiFtPJ+OXg44Nfn2sRRT6g+fAbqAWYjsQG6ARuYjcQG6ID2aQ\/dBTGgiaPcGf5brDLgOHGMl3qQZVWiLoAeaqMLZzFWJ3UtC1yCBm772AaX76v\/\/u\/\/vv7dv\/t3+p3f+Z1W8EKVSfDrIP6\/\/bf\/pv\/8n\/\/zbXXja8Ibb7yhe96+R1+\/8XX9qf\/8p6T\/Wo5qhOK\/Fvpnb\/wz5f9vXr04x63frY86WHc3b96s4vBZ\/9DdcvA+h+3X4dr7+\/vipTKUnYQXbr5QDqD+x1sDf\/DGD477RDsGvx54v58Wg68Jbj8tPiQpb\/Wju3DrMN7qwLZyTRwxfh2u\/c1vfrNOyJp\/bt2kXKgQEzFS5YI8r\/OlVv97+ebL4s1cl3TrP5G4+OZFXS\/\/mJg1WviCgJFZCd+GhAPoBt+GhwPoBt+GhwPoBt+GhwPoBt+GhwPoBt+GhwPoBnYgpiPxA3QDNjAbiQ3QDdjAbCQ2QDdgA7OR2ADdgA3MRmIDdAM2MBuJDdAN2MBsJDZAN2ADs5HYAN3g5874RUgmXv8ytMvhw562L3feeac+9KEP6QMf+EAreKHKJEyqw\/VPqsv8P\/bnf0xPvfcp3f+D90s\/WI5yhOvvv65\/83\/8G535wJlxv61Mk3Tbn6Q31WH8O9\/5Tp04cUJ\/4S\/8hXH7oTr\/2Af+mL7zQ98pO17\/e2TnEVkdvgyVb+L8sk12U3mfbyrv8n6ZNtstF9Lvv79cn3VK1vpz6ybloYbi7Hdnf0cnHz0pFuxMmfDxswn3zVxXdEWX77wsfzl58mS1Ibg8nGujwwF0g2\/DwwF0g2\/DwwF0g2\/DwwF0g2\/DwwF0g2\/DwwF0AzsQ05H4AboBG5iNxAboBmxgNhIboBuwgdlIbIBuwAZmI7EBugEbmI3EBugGbGA2EhugG7CB2UhsgG7Y5Dd6\/aE\/9If01\/\/6X9dP\/MRPtIIrB5MwqQ7XP6ku81Pms3\/rs\/qPn\/yP2vvknvRJjXHwtw70tZ\/4mv77T\/x3EWdlmiQxXdFUh\/EczLzjHe\/QJz\/5yartpnq\/+xPflT6pakmUaK\/8szp82VRHiPfLNtmhsiGuqbzLh8o1cW65kH7PPfdoE5atmpSZfDkD5ov81P5T2n1lV8\/uP1utR3jOoCuj\/OCLztkxPzHgC1BSR\/7BAZd07YODg417oxeX7ewSLuNmB4I0MH5gNhIboBuwgdlIbIAOyB\/t7ezsYI5BDBgTpYINSlWU401clDMOHmAD9ElxxFCPjZdyAD5UtgtnMdRpdVHfssFZMWfHs+\/H7Gt0c9in9qZyHIzbdu7Wx\/7gjM6oKP80WprqwB3ydeFCMdTXBPrDyYP52T+xPzN7kuzbntXXtVxbXF9fKN7nfNv6u85yayZlvsyP6\/FqXV3QBT1c\/mHwhQboAH1Pe2JjVVxiBmIGNj4DbPM8L8J27w6WfQYT8yrtCzixcPvIvsy1o77+GdiKSZmNiwkZyZElG9nhqUO9efFNPXjqQcFrtLCB4h+ZUwvOhNwHhbpW1LXcpLgmf1fej\/PtruOZNm7a9rqWmxTX5g\/5unChmGnzs+xyyzqznjaHXcqx3XPWzD7AzS9nzed1XlYHV2FcP3rI14ULxVBfCOynntEzYxf7MjAmOih92nOr61quLa6vLxTvc77t9nld9a2YlLncwxEmGxv3iFlZ3zr1LX3n4ndQBW9HycSxcSouMQMbmgEmVC5ZG7BtqCHOfNsg2RcwMdv+wMbMfoGzZqRxi5ZMym77fh8X3Z\/Y3nwysPGTMl9kjnRJH19iNjr09197v3av7QqbS1dMxDZhczRKOcUlZmBDM8BEbHCHGOJc\/7bo7A\/YL7jjZZ\/AGTNX3Vx+Ubrtx2iP\/Vbfs2TKNSN6ViUDGz8p8+VlMibh7pca+8RPntCL+y+KGDY483OfBk7HWKZ9AKFruUlxTf6uvB\/n28dITaei07bXtdykuDZ\/yNeFC8V0SkYMGmdg2hx2LefGMelx1owcd6BUhhrKP2t2y2m0dOFCMaPiRwRnyMBI26eZ3VV2bc+vr2u5tri+vlC8z\/m23+91tDd+UmalcNTLhsXEy5EuHNjZr5\/shefoF0kc8fgjYgZiBrY7Axycc8YM3Eywr2Bfsqh9hZ0w0Af6lCpVXCZnYB0jtmJSZsXYpWmONv0NiSNfeL7sFkeZ42B3d1fxQa\/pM7iq+WNEob514UIx1BfRPQPT5rBruaa4gQa6Xv4h5SxMlpw1X9u9dtv2HqrL53zbqXqscgDA\/skIvw\/Gd5Fd2gvV07VcW1xfXyje53w71Pd147ZmUmbCtaNdNqQX9EK1rviyX9Il4WdCRiouMQMxAzEDXgbYN7AP8S8dM2nO86yZ\/ZN1hT747Zsvys3IwORJeTPGWY2CI0z7Qv\/8uZ8XP4v6\/KnPVz7uI+OvjPgRMxAzEDPQkAGutLWdNTNJNxTtTVMXV\/KsIPsoJmazo9y8DGzVpMzqY4Pii82X\/cZv3BASG14zXKZ9AKFruUlxTf6uvB\/n2zNMVbCqadvrWm5SXJs\/5OvChWKCg49kYwamzWHXcm1xro+JsemsmcvZ7E\/ceBuQz\/m2Vbd\/egAAEABJREFUxZnkSp7pSE4ekNNiUntN9XYt1xbX1xeK9znfbur\/OvFbNymzcq7oCqICGxeXrSuj\/Dg8PKz+R6BSHf+DGxsjBQ6MzEr4NiQcQDf4NjwcQDf4NjwcQDf4NjwcQDf4NjwcQDf4NjwcQDeE7MiFvz9+XiyHUa53Bph8OWtmP+KOhEvO79t5n0v11jlh4OeZVpCTB2B2lJuZga2clNmAeOc1T18zIWPb6uV\/8wFmI33bOJ93bXsAgfc3u7yVRbogBli5nZ2d6r9ic2PQiQGT4uy\/NKMeygG\/nHE+79vEherz\/6cjygHiDdjAbCQ2QDdgA7MZH30nf8YhiQHoBmyATTkesKOccfAAG6BPimO81EMfiKccQA+V7cJZDHX6uaPeiMkZcHM4OfpWRNdybXFNPvYfTMx2a8xafWnnJb33zHuFNM6vw7ctDsmk7J4p+\/UT0xdt7bXV1bVcW1xfXyje53y7bQzr4tvKSZmVc9\/+fTrx5AkhsQ38bz7AbKRvG+fzvr3wOBocgb6AkVkJ34aEA+gG34aHA+iGEydOmFpJ\/KAyRh\/YYGRWAhtUxugDG4zMSmCDyhh9YIORWQlsUBmjD2wwMiuBDSpj9IENRmYlsEFljD6wwcisBDaojNEHNhiZlcAGlTH6uOuuu0ZaFJuSAc6auaTtjoeJlYfAgMuH9P2dfT168tGxi7NtM5j441myZWOz5dZOyqzW3Vd2EUfAkRdwSd\/GBwfQDb4NDwfQDb4NDwfQDb4NDwfQDb4NDwfQDb4NDwfQDb4NDwfQDfP8X6Jog\/YAugEbmI3EBugGbGA2EhugG7CB2UhsgG7ABmYjsQG6ARuYjcQG6AbfNj7K9c4AEydnzf7viDnj5V4zsmmEX9z5ol7ZfUVc9mYyd2NncZbc1G7kVysDWz0pz3NVTPsAQtdyk+Ka\/F15P86355k76nbaw+yMruUmxbX5Q74uXCim88BiYJWBaXPYtVxbXFcfZ7XcFnv+8Pmqz\/bBRMsZ8yff\/OT4uRWr8xcPf1H2lDWTMRO4Rgv1+ZP8yNVbWHt9n3HoWq4trq8vFO9zvt07IStYIE7KK7hSYpdiBmIG1j8DDxw+oN+48RtCuqN5+u6n9e6dd4vJV+XCZeuP7Xys1ML\/OPsOeyK7iRmIk\/Kc1iqXJ90Hhbo207XcpLgmf1fej\/PtruOZNm7a9rqWmxTX5g\/5unBHYqZNzJaXmzaHXcu1xfX1Ef\/g6Qf14s6L8i8\/21nz5d3L+sjpj4zX6tmDs2MdhbPk4\/4MinoM9GlZ+6W2tkO+Llwoxsa6rjJOyuu65mK\/YwZWIQPDofT442Oc++pX9eU33tCpz39+zFX+4XAVeru0PvAQGPea\/bNeHuZigqZjHzr4UJm7L6MegV\/miDMaG5eBOClv3CqNA4oZWGAGXnhBGg6l4VAaDit88K23KqnhUBoOpeFQIm6B3VrFphIl4uls\/6xZo+W+nfuqB71GZiWYsJnQKyN+bEUG4qQ8p9U87QMIXctNimvyd+X9ON+eU9rG1U7bXtdyk+La\/CFfFy4UMx7wuitJIg0G2r\/vPl3b3ZUGA2kwmPmops1h13JtcX19oXi4C9cvyH8IjETx9PWTJ55EPYJLuqS8\/NMMFtq\/fv264oNeM0jmnKr4vjnVG6uNGYgZ2KYMJIl09ar2n31WT999dyWxNRisfxZmPIJJD3ZZc6lS2fK4ylsEZkS50RmIk\/KcVu+0DyB0LTcprsnflffjfHtOaRtXO217XctNimvzh3xduFDMeMAbpFT3lE+fls6ckYqiHllRSMOhlOe1PeXntDnsWq4trq8vFA\/366d\/XVyWVrnYJe1SPfIPnp9VIXEQnynTcRfajw96HTeL8y0fJ+X55jfWHjOwVRk4dXgoUA26KKSiqNRqMuaBsPPnpTvuqCdsJm24LJOGQynPtQ1LpkzuZKty4SGwVKlY8O1pT\/xuuVCdPzjKaXZLrGlFMxAn5RVdMbFbMQMbnYGikIpCGg6lS5dUPaHtT9jYWSZlmZTnUlFsTEo4C7bBnNd5Mfn+3MHPVb9r5n7zQAPZJetEidx4xWWjMxAn5Tmt3mU+UMGQmtrvyvtxvk0b88S07XUtNymuzR\/ydeFCMfPM4TLq3t\/Z0UOnT2v\/qaekvT0pSfp3oyikopDyXNWEzaTNBM2ZdXmWfVjWL+wsk7JMynOpKFrb6Zr7tri+vlC8cbxQZE9lfka9tgn41OGpimGirpTy44IuaFD+aQaLtb+xD3rNIEfLriJOysteA7H9mIFNyECei4ny3E\/9lL7E75S\/9S0d+RlUmkpvvy1dv149EFZN2Ht70mAgJUmvDOzs70t5rtCEXd3HZsJeg8vimTINyj+VC2fKvPO6VAWPrXJJlQpbcdmaDMRJeU6repkPVDCkpva78n6cb9PGPDFte13LTYpr84d8XbhQzDxzuPC681zKc507OJBKWaEojnYjSaTBQMoyKcvqCZqJumnCHgyOlp9kFYVE28OhqkmbybmcpHfvvFNn3vte7bz73aoulWeZNBxKxDp1tq2jvr5QvM9d1VUl5Z\/K5aWdl8RPov7Bzj8oLQl+1pet\/fbVcelari2ury8Uf4QbDrX7iU\/ozN\/5O9r52Me0KS+uCUzKHddSDIsZiBmIGXj4YSlNpTSV0lT7jzyiX\/3+75dK\/QiIm5StJJEGAynLpCyrJ+yrVyV3wr5yRdrbk6h\/MJhU41F\/UUhFIZU7c3fCVnlJvDrDPnNG4wk7yyRvwj5a2ewsd+L91e\/\/1ariRIlcXnG5PQO8kIZ16WATXlwTJ+XbV3VkYgZiBrpmIE0lJsoR9r\/whfJs74SQLq807VpjOC5JpMFASlMpy+o23Qmbs236wIQNBoNwPU1sUUhFIbGD5x42OD96Utwmbewsk7JMYsIuiqbauvG0VZ7JDx5\/RnvD5EiZC1kh+OogYTg84ouGl4EkkQbzf3GN1+rczK2ZlOeWwYaKl\/lABV1qar8r78f5Nm3ME9O217XcpLg2f8jXhQvFzDOHm1h3MIdJUl7rLZGmUpZJWVafZXOGDcoJ++DXfk1vXryow89+VhoM6vg+CSoKKc9VnWHbhM2ZtTNhv\/nJT0pZJuW5lOcK9fUI55zpZY8XGuR1h9JyDs4ulfpwKA2HOnJvvqSP8+9I+z0q6lquLa6vLxQf5N75Th0+\/3z1wprHTp6s5Dq\/uGajJ+V7771Xhh7fvxgaMxAzsEkZSBJpMNB3yklZWVZP2OVErdGEXe3A98pL4iUOH3hAB2fP9ht9UUh5rrufflrVpM0ZdQnuY59+6CHtvO99qs54s0w7v\/zL2r127Uj9h6dOSYOBnv\/8A3rgtVO68sxAGgyOxKy8URRSUYiH8AD6GGVuyM\/OSy\/VY89zKc+l4VAaDqXhsMrL9\/9qeem+1LFFrp56SneVOd0pJTbgHvKJJ5+s7iFXdZSJ2X3llSrHPGT43M2bJbPe\/zZ2UmYy\/va3vy2wjFV05IGEHh3oWq6OO6OdnZ1g7U3+rrwf59vBRmdITtte13KT4tr8IV8XLhQzw5RtRVXT5rCxXJJIg4FUTgJg58UX64nDnbDtsniaSoNBrzxXE1SeS0w25Rk2DySdfOwx7fzBPyjBlbVVMQ8\/LJX4R\/\/Pl3T4C79QX54fDEpv+a+c7JTnUp5LeS5RrgmjcciVjz9eHxSUkknNHoyqDhTKgweemr8NXAVwsHv\/\/aoelKPfXB3wcaa8H1+CB+mqAxEeqCttgVEbHJxUY+cgBa7sT9WHUpKX8WRb2hzc7Hzxi7cOdMrcwTH28eRNXsr0VP\/yXMpHDxlWxPp+rPSkzMQaSi08CPmMs8mYONPNF2XMQMzA9BlgmwLT17AmJZNEGgykNJWY5Jicm+5jp6k0GEw\/sHLSYRKqJi2b0MpJpqoQySRmYNJqQllPNXm5smkCh6fuEJjwQqg6FD\/mmYGVnJTZ4EFo4PBMsgCdGKQLOIMbZ1yUx8\/A66+\/rmeeeUbI49d2tIZorW4G2M7YpgD66vZ0zj1Lkvq+dJpKkybs8rK4BoM6fs7dWuvqk6TOUZLcLgcDaTCQBgNpMJAGgzqGASeJlKYaP\/kPt8ZYyUmZDR74eWUn4PLoxqEbKAePnIT9\/X3NA0V5lJmVG+u18v5Rn\/q7lpsU1+Tvyvtxvs1k\/PWvf30uuSNffntwXdC13KS4Nn\/I14VzYyZ9L1fRzzbFNmZ9Q4cz25Vd1tU0MW4O+5TvWq4trrOvvKXEm82KwUCZpGvvf3\/18NF+eWm8KO9lf\/GjH9W1y5el0q9yIRaU6lT\/KNuI++4T\/52mC\/5rTSawIMq+7QdQlAcW\/9cf\/+O69vGPV29r2y\/v847x7LP1+EpZlFcTfrYcL+NjvC6Kq1dVjf3ll7V\/48YRkJfswgVd+8Y3ZGWKq1c15sq6+R\/IirL+rLzMf+1nfkYa5c\/uXfPO9T9yeKhTH\/6wqsvx5T54qoQuudBKTsqzyIntMJp2GqdOndLZs2f1WHl\/56GHHtKs8d73vldMWn3r71puUlyTvyvvx\/k242I9IWedO+rz24Prgq7lJsW1+UO+ydxDcmPIG98\/vofkcVPAeBgX4+uyvvrGuDnsU7Zruba4vr5QPNzXfv3X9dhXv6pf\/bf\/tlrtTKi8mvTe8v4r0sDkSQASDn8I+Brx+7+vhzzwhPJDv\/3bCqLs20MBvLc8AP+\/\/8f\/0GP\/8l\/qoa997Sh++qf10Ajv\/Tt\/R1\/9T\/+pGt9DH\/mIXLz3x39c1dgD+1zy4u8vJ3G\/+k\/+CempkedSPrqnXEr0arKuvWv1ubGTMmuBidmA7YKdx5e\/\/GU999xzETEHS\/kO8P1zv5OboMftqvv+5Oy5c9UqP3dwIJ4aBryi1MCZHwGnHnhAX\/pH\/2gp39FV3j+e\/cxnpDSV0lRKUylNpTSV0lRKUylNpTQVD9CRx3XBRk\/Kk1YCO5Bz5YYRcU4xB4vPAd+\/Sd\/RVfeH+se44vdp8veJPFn+mJh9jCfl8qpezOft+TzFb86vXKmfVG+TaWppXgu51ZPyWqyh2MmYgZiBzcxAeW+0OpNLUylNpTSV0lRKUylNpTSV0nTtzvQ2c2UtblRrNSlzKdq9R4wOt7h0xZZiBrY3A2xrbHOWAXQ4szdHLmgkaTr5LI8zwDRdUIdiM6uQgbWalEkYOwF2BgAdLiJmIGZgMRlgm2PbA+iLaTW2EjOwPRlY6Um5aaOHB9uzmuJIYwZWJwNse2B1ehR70iUDMWY9MrDSk\/IyUsgZgGEZ7W9Km+RwU8ay6HGQO8Oi255XezYe5Lza2PR6Y+6mX8PkzjB9LYspuRWTMisjlE54YD50zgCAcVHWGSA3tXb0Ex4cZaPlZ6ApR\/DA4tH5\/hmMX0VJX0P9ggfmQ1+H8Vh\/FyXJS6gteBDyRe5WBuoc3bJNgweubd8\/pPGrKjd6UmbFgFDy4VlBAJ0YdCS26djbDHIBQjmAJ08A3WJc3bhtleQChA9puTAAAAlESURBVMYPT+4AOjGmYwO4VQP9AqF+wTMGgE4MOhLbdOxtBXkAofHDkyOAbjGubty2SnIBQuOHJ3cAnRjTsQHcKmOjJ2VWBvBXACvG5dHhLM63jd9GSS6AP3by5fLocADdj99Wm1wAf\/x+noiBIw4doK8i6Bvw+0b\/XR4dzuJ82\/htk+QB+OMmVy6PDgfQ\/fhttckF8Mfv54kYOOLQAfqqg0l51fu4sP7ZClxYgxvakOXR5IYOcy7D2sScbeKY5rLyWyq1HJpsCY0uLwPrlrM4KTsrkCMpViBw6Kj2yAA5BBQxiR7RLQPkjO8fQO9WarWjGAfjAavd09XsHfkD9M4kekS3DJAzvnsAvVup5UVt7qQ8ZU5ZaYYpq4jFygyQw1LEf1NkgNyBKYqubBHGY1jZTq54x8jfindxZbtH7sDKdtDpWJyUnWRENWYgZiBmIGYgZmCZGdjKSZkjJi5lWOLR4cxeYbkyXSNf5M06hA5ndpTtGSBX5Myi0OHMXkdJ\/xmH9R0dzuwo2zNArsiZRaHDmR1lewbIFTmzKHQ4s9dFbuWkzMphZbHSADpcRL8MkDfyB9D7lY7R5IzcAfRNyAjjYDwAfRPGtMgxkDNyB9AX2fYmtEXOyB1AX8cxbcWk3LRy4ME6rrhF97kpT\/DgWP3ZgsJNOYIH65iCpn7Dg3Uc0yL73JQjeLDIvqxjW005ggfrOCb6vBWTMgONiBmIGYgZiBmIGVj1DMRJedXXUOzfumYg9jtmIGYgZqB3BuKk3DtlsUDMQMxAzEDMQMzAfDIQJ+X55DXWGjOwmRmIo4oZiBmYawbipDzX9C63cp5AXG4PYuuLzEBc3\/PPdszx\/HO8Si0sY33HSXmVvgEz7AtfpnV+AnGGqdiaqljfrPetGfDsBtqpJnJLjjsFx6CNyADrm\/W+yMGs3aQ8KUGu39UXmdRQW\/TFhxuHz7VnrbfV3+br249Z1jVN27QP\/LIhzo2Z5Ce2Swxx02Le9Tf1a1K7rt\/Vm+pbFE9ffLht43PtWett9bf5+vZjlnVN0zbtA79siHNjJvmJ7RJD3LSYd\/3T9qut3NpNym2DWVUfXwyOuHzAz6PP1Etbft3wPrcpNmNjzIbQuIjpw4diV51j\/E3jXPW+9+kfY2SsPuD71NM1lnppy4+H97mVsY\/ZEcbGmA2h6ojpw4diV51j\/E3jnEffN25SJoHzSNS861xEv2ljkV+ueecs1r+YDPC9WUxLs21lEf2mjbhNzXa9bXttGzcphzYQOIOtcGx0JHB1bABnMBtpMJ8r8bl2m84GbX4rh\/Thxrg+42clqduvyzhXmk4sOkB3AWcw3myTxiPhTJqO7QLehfng0JEAPQTy7fux4S0e24XxvnRj0Lv6m2K7lieOOlzAzRO05dcPZzAfNjoSuDo2gDOYjTSYz5X4XLtN99clsZT3AQ+aeHyzAPX79RjnStOJRQfoLuAMxptt0ngknEnTsV3AuzAfHDoSoIdAvn0\/NrzFY7sw3pduDHqDX\/iA+V19EkeswWKRxpmEWxY2blL2E0mS+YIYsC0G3eVNN4nfYpHY5kNiwxuw4c02CYfPYHxIEuvCYijr8uhw5u8qpy1H\/bRHeYAO0AE6MQAdzoANzDYJR7wB23zGmXR9FgOHH9uk6dh9QX2UdwHn1wPnxqDDWRw6nAs480+SxLpl0eEoh8R2AYdvUaC9pvZdn6tbPJzbT2zzIbFDfpdDt1jiAVwTiHVhcZRzeXQ483eV05ajftqjPEAH6ACdGIAOZ8AGZpuEI96AbT7jTLo+i4HDj23SdOy+oD7Ku4Dz64FzY9DhLA4dzgWc+SdJYt2y6HCUQ2K7gMO3DGz8pEyimxLr+ly9S3xTTBNP\/QZWOGiKNZ4Yypg9K0md1N23Psq5ZXzbfE28+ZvktOWa6lslfp5jm2fdoRy2tef6XD1UD1yXGOJCoKyB7zMIxbkcMZRxuVno1EndfeuinFvGt83XxJu\/SU5brqm+VeLnObZ51h3MoUNu\/KTsjHWiykblYmKBKQNY4YC2mqrAR4zrh3Ph+tZBd\/uOvow+k1NrG4k9TT8o62KaOtrKuHWjWyz9xXZhvlWUbj\/R59VH8gLa2sBHjNsHOBeubx10t+\/oy+gzObW2kdjT9IOyLqapo62MWze6xdJfbBfmW4aMk\/Io66wQVo6LkatREEs5ApDY6D7w+VybTXyoLjgfbfW0+aiHdtpiZu2jTR+zbmNR9fnjwG5ru2+uqc+H1e\/zfeu2euYt6Zff10ltEk854pDY6D7w+VybTXyoLjgfbfW0+aiHdtpiZu2jTR+zbmNR9fnjwG5ru2+uqc+H1e\/zfeu2emYhN2hSDqdj2uROWy7ci+4s7fIF6VKC2FAc5Zt8bvykuC51uPWZ3rVc1zird5XltGPpWs7iTDblAj\/rtck\/C542pqln2nLTtOWWod2uOSHWLWs65Zt8FoOcFNelDurx0bVc1zi\/\/lW0px1L13IWZ7IpB\/hZr03+WfNrOSmTpBBCySGZbiz2ceL8stRH\/UjfZzY+YnzAW4wr\/Th8xPo8HL7jwK0D3W0De5q6KefX43P44brWTyxlXMB1Le\/GUY56kCEen8GPIR7O\/Cbh8AF0403ChXzG4zPAWTmTcPiRxpmEw3ccWF2+DNVJe24c9nHi\/LLUR\/1I32c2PmJ8wFuMK\/04fMT6PBy+48CtA91tA3uauinn1+Nz+OG61k8sZVzAdS3vxlGOepAhHp\/BjyEezvwm4fABdONNwoV8xuMzwFk5k3D4kcaZhMO3DKzdpEyymkAC8SGB6UgDPMBGuoAzwKMjgatj9wXlfbh14MNGukCHB+gu4JpAHF8w82Ob7kvXh24gDt2V6MB4dGC2K9EBfoDuAs4Ab3qTJMaFGwfv2r7u+5tseBduPfBmo7sw3qTrQzceiW0w25WmWwwSzoDtwnjWN7zZXSVlmkAd+JDAdKQBHmAjXcAZ4NGRwNWx+4LyPtw68GEjfcCDJh6fD2LJsfHYpvvS9aEbiEN3JTowHh2Y7Up0gB+gu4AzwJveJIlx4cbBu7av+\/4mG96FWw+82egujDfp+tCNR2IbzHal6RaDhDNguzCe9Q1v9iLk2k3Ki0hK1zZYYWDRK61r\/1a1X137H+P6ZWAT1jfbE1jVsaxqv\/p9U2J01wwsY33\/\/wAAAP\/\/XWNhjQAAAAZJREFUAwD+cXK\/wz+99AAAAABJRU5ErkJggg==","height":234,"width":388}}
%---
%[output:40ff28ed]
%   data: {"dataType":"text","outputData":{"text":"[23:23:14][INFO]  Top-5 hits (Tanimoto >= 0.30):\n","truncated":false}}
%---
%[output:91bdcc14]
%   data: {"dataType":"text","outputData":{"text":"    Rank            Name                           SMILES                   Tanimoto\n    ____    _____________________    ___________________________________    ________\n\n     1      \"IBUPROFEN\"              \"CC(C)Cc1ccc(C(C)C(=O)O)cc1\"                 1 \n     2      \"FLURBIPROFEN\"           \"CC(C(=O)O)c1ccc(-c2ccccc2)c(F)c1\"         0.4 \n     3      \"KETOPROFEN\"             \"CC(C(=O)O)c1cccc(C(=O)c2ccccc2)c1\"     0.3947 \n     4      \"AMPHETAMINE SULFATE\"    \"CC(N)Cc1ccccc1.O=S(=O)(O)O\"            0.3514 \n     5      \"AMPHETAMINE\"            \"CC(N)Cc1ccccc1\"                        0.3333 \n\n","truncated":false}}
%---
%[output:14cbbcd5]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAeUAAAElCAYAAADa2PrWAAAQAElEQVR4AeydTahd2ZWYd9eTS082HRAVIQtcQWiUDAMBaaSJRiFo2AFrYAuUjLpBhGBwZhkIEhAEPAwtxdgEDULSNhG4oeENO3QZQtOGDGKIOm6pSoIyxnKVfpqopNR3n1b11r7n3N\/zf77C6+29fvbea3373rP07nuW3nvjfxKQgAQkIAEJDILAe8n\/JCABCUhAAhIYBIFpNuVBoDUJCUhAAhKQwHYEbMrb8TJaAhKQgAQk0BoBm3JraBvf2A0lIAEJSGDiBGzKE79gy5OABCQggfEQsCmP566mmalVSUACEpDAVwRsyl+hcCIBCUhAAhLol4BNuV\/+nj5NAlYlAQlIYCcCNuWdsL276MKFC2mdvLtiN40zdlvZ7SryLKUqA2Kq7Pva8n3z+b775uvb2jc\/I59zHpLb8jm+kNwe8\/Axhm3VSFydxLoqf\/jqxlhT5y\/tEZ+PZcwU9ap6sQ25VvILGXKeQ8\/NptzADT148CDlwpa5zhzbnISac+HNOqf6m6wVdsGSebk3tvAzoucx6NhD0HN\/3Tzivxrfvs7z+NK36d75HuvmXZyxLoeu\/dTMmTEyH7Jw7+QaMuRch56bTXnoN5Tlxws+Uxuf8sZqfNO3G5J7uT+2t+61Q7l21YJt9l21T\/iqzm76jDirHDk7P4s5tohjji10RnTszBnRmYegYw+9yXHV3pyJH2G+67n7rt\/1XNdJoAsCNuUuKHvGgoAP0wWGyX\/p4p4bPGPy92GB4yJgU+7wvvjuIJfyaHzYGEPQQ7CVc2xIbs\/1sDNiD0EPyW35vPTjC1sTY74f81xi\/7Ch5\/PQGZHcxxxbLthCSnuuMycuxqp5+BhDiAsJW4zYmTOGoIdgi3k50nxKm\/oygZwhc4SoGJmHlDb0UoitsmGvEmLX2YnJpSp+E1vdHtjL9bkt5owheTy20JkjpZ7b8KEj5TzXw48tJGyMIeGb+2hT7ugVwAuPB2wu2Mrjsa2LiTV5LHMk1jKPOEb08DGiY0fQkXKOThy+EHTsTQp7xv4xYuOM0Ms5OkIcksdhzyX3E4ee++vmxCL4GRHmpbAfvhD0bWJYV8bX6ey9TXzdPpvYOauUTdZVxbBPnjdzbFWxu9rYD2FvZJN98vhYw1jaw1a3Z5U\/9mBNzIkLwYZvG2FNrGdEr1xfYySedSHoeSg6En58uY4dHTuCjpRzYrCHoBOTC7bwM6Ln\/rnObcod3Twvuk2O2jSOvcrYUicG4cVe+tCx468T\/MTlfnTsua2NOedssi9xyKrYdf5Va9f5YFHuj449X4st13eZs2cT+2x6NmeVsunaPuIi17bOZv8m995lv13W5DmvW48fiTVVrzn82COmHPERk9vRsZe2XHd+TMCmfMyhk6+8KHPp5NARHBJvWNnUXxZs4FQfMVwPuZMdYy5hYxyawDrPlfk+OTa93z65dLEWXqU0cO4strApd3TNvEB5Y+bS0dGDOSYYVCWUc2FObFXcHG2wgMlYaq\/Kl\/yrZNeaqs7YZS9yYq8Q9NiHeS7EhK9qJDZiGNHzOPRciMn9U5rndebzKdXYVi025bbIrtm3yzckb4ryPHTsq9LET1weg449t20yX7UO3yZ77Bqzzf7bxJIPLMo16Njx7yvr9uIcYvJz0LFjY0RnHoKOPfQmxzb3jjz3OYO1sQ8jOixCsCHYGZuSpvbbZp+q2CrbqhrhUq5Bx163Dh8xpb\/KVsbMUi+KHl1T5mKrJK8r9+f2PufxQo3c0LvMh\/PibEb08nxs+JDw5Tbs6OFbNRKby6p1+PJY5tjy\/dGxI7l9k3m+lvXosY45thD08MWILfxhy8fcTxx67l83Z82qGPyl5PGcl\/vRt\/HnsdvO83OZl2dX2fIziCcmt5Vz\/Lmwpoyp04ldtzb3M2evch12bPhWCTFVsWHHF4Ktbi9i8MXIHGENthB07AjzsDOiY88FG74Q9NxfNScm4hnRq+JyGzHE5oItj3FeTWB0TZkyuNxSsCO8CHIfOvYuhfOrzsMegp85Y0ipY89tdfMyrkoPG3sg6FWCD8l96CG5vW4esflYFYs\/7MxzCXs+hj9s6DHPx9wec8aQPJZ52BlDZ8wFHxK2fI4NPQQ9F+y5zjy35XN8ueCrkjyGeR6DXso6f1V8aSv1fM+YV8WUtlJnbWkLHV8p4ctHYnI9n+MLwc6ckWcD81Kw4y\/t2DYR1lXFYc+lKiZsVXHY8DOGoOcSdkbsMTIPwRYSNkZsjFWCL6TOX9ojPsbcjy3XmVfZsM9NRtmUly\/p2MKbqbxYdOzHEX6VgAQkcEwgng08H3LBfhzhVwl0T2CUTTl\/A3WPzBMlIIGpEKABlzKV2so6qLO0qQ+PwCibMi+uEBr0rlgfPXqUhizm5v34GvA14Gtg+9fArj1hCOtG15Rpxjk49F0aMy\/0733ve+ny5cuKDHwN+BrwNTCh18C3v\/3txTdcea8Yy3x0TbkpsDTljz76KN2+fTvdu3evVm7evLk4cl3cqj2a8P2Lf3s3\/fN\/88eV8q\/+\/Q9r82\/i7HKP\/ZjUsy7PGZsul+q7lYtcNn0vN\/FaYQ+e7YsH9wi\/jK4p7\/Jd8ap7uXjxYrp06VKt4Gc946q4tn3f+kf\/JH34D+ul7fPz\/WExBCZ5TkOYy6X6fSQXuWz6\/mzitRJ78Iwao4yuKZeQadJ8hI2dEZ15CDr20Lcdv\/nNb6bvfOc7iXHbtVONh4VM3r1dNLlAYVnksswEi1yg8K7IJKXRNWUaLI02BD2\/VvTwMaLn\/m3nvEi++93vbrts0vEyqb5eucilmkC11dfLMheZjLApc4002hD0UsLHWPrUJSCBTQkYJwEJdE1gdN8pdw3I8yQgAQlIQAJdEbApd0XacyQggUEQMAkJDJmATXnIt2NuEpCABCQwKwI25Vldt8VKQALTJGBVUyFgU57KTVqHBCQgAQmMnoBNefRXOPMCbtxI6erVYzk6mjkMy5fAtAjMsRqb8hxv3ZolIAEJSGCQBGzKg7wWk5KABCQggWkSWF2VTXk1H70SkIAEJCCBzgjYlDtD7UESkIAEJCCB1QTG2pRXV6VXAhKQgAQkMEICNuURXpopS0ACEpDANAnYlId0r+YiAQlIQAKzJmBTnvX1W7wEJCABCQyJgE15SLcxzVysSgISkIAENiRgU94QlGESkIAEJCCBtgnYlNsm7P7TJGBVEpCABFogYFNuAapbSkACEpCABHYhYFPehZprJDBNAlYlAQn0TMCm3PMFePyeBO7eTen+\/WO5cmXPzVwuAQlIoF8CNuV++Xu6BCTQNgH3l8CICNiUR3RZpioBCUhAAtMmYFOe9v1anQQkME0CVjVRAjbliV6sZUlAAhKQwPgI2JTHd2dmLAEJSGCaBKwq2ZR9EUhAAhKQgAQGQsCmPJCLMA0JSEACEpgkga2KsilvhctgCUhAAhKQQHsEbMrtsXVnCUhAAhKQwFYERtOUt6rKYAlIQAISkMAICdiUR3hppiwBCUhAAtMkYFPu9V49fG8CN26kdPXqsRwd7b2dG0hAAhLok4BNuU\/6ni0BCUhAAhLICNiUMxhOmyHgLhKQgAQksBsBm\/Ju3FwlAQlIQAISaJyATblxpG44TQJWJQEJSKB9AqNuyhcuXEhIiQlbSOlTl4AEJCABCQyVwKibchVUmvGDBw9SCHpVnDYJSCAlGUhAAsMiMNqmTLOl8eY4q2zEYM\/jnEtAAhKQgASGSGC0TbkpmI8ePUrIq1ev0hAl6nzz+nV6XSdvXi\/Chph\/2zm9efMmhbR9lvsP8z2yfC\/mOUcmPMeRxcNwxF9G2ZT5zpfvgJvgfu3atXT58uV069at9PDhwyXhkp88ebJkr4ptw\/b48eNFmc+fP0\/Pnj2rlJcvXi5iiG0jh3LPvpnk+Xz++efpxYsXC\/n1p5\/2dk\/kNCQu5DMUkcvyc4W7kcsyl32Y3LlzZ\/Es55m+eCCO9MvomnKTDZk7u337drp37166fv16Onfu3JKcPXs2nTlzZsleFduG7fTp06SZDg8P06lTpyrl5MmTi5iu8uybSc7569\/4RqJ+5IMPPujtnshpSFzIZygil+XnCnezLRfWTF32YcIznGf5zZs3F8\/DsX4ZXVMGNI05JHTGXeTixYvp0qVL6fz584vGR\/PLhUbIAz+3dT2nroMTJ9KJGjk4OCCkMv82ch0Ck6jr4L33EvUjX3v\/\/c4YxPn5OCQueV59z+VyWPm6lMsyl32Y8AznWc4zffFAHOmX98aWNx9b50L+6IyKBCQgAQmMiYC5lgRG15TLAnKd5sx30LkNHXtucy4BCUhAAhIYIoFJNWUA04BpxCHo2Icsn708\/u3pIec42Nzu3k3p\/v1juXJlsGmamAQkMB8C+1Q6+qZc1XSxhewDp6u1v3\/4Xjr65Yv00188W5Kf\/+r4N6u7ysVzJCABCUigPwKjb8r9oWv25I+fvkpV8snTL5o9yN0kIAEJSGCwBIbblAeLzMQkIAEJSEAC7RCwKbfD1V0lIAEJSEACWxOwKW+NbK8FLpaABCQgAQnUErAp16LRIQEJSEACEuiWgE25W97TPM2qJCABCUigEQI25UYwuokEJCABCUhgfwI25f0ZusM0CViVBCQggc4J2JQ7R+6BEpCABCQggWoCNuVqLlolME0CViUBCQyagE150NdjcmsJ3LiR0tWrx3J0tDbcAAlIQAJDJmBTHvLtmJsEJLAJAWMkMBkCNuXJXKWFSEACEpDA2AnYlMd+g+YvAQlMk4BVzZKATXmW127REpCABCQwRAI25SHeijlJQAISmCYBq1pDwKa8BpBuCUhAAhKQQFcEbMpdkfYcCUhAAhKYJoEGq7IpNwjTrSQgAQlIQAL7ELAp70PPtRKQgAQkIIEGCQyoKTdYlVtJQAISkIAERkjApjzCSzNlCUhAAhKYJgGbcsv36vYSkIAEJCCBTQnYlDclZZwEJCABCUigZQI25ZYBT3P7AVV1925K9+8fy5UrA0rMVCQgAQlsT8CmvD0zV0hAAhKQgARaIWBTbgWrm46RgDlLQAIS6JuATbnvG\/B8CUhAAhKQwFsCNuW3IBwkME0CViUBCYyJgE15TLdlrhKQgAQkMGkCNuVJX6\/FSWCaBKxKAlMlYFOe6s1alwQkIAEJjI6ATXl0V2bCEpDANAlYlQRSsin7KpCABCQgAQkMhIBNeSAXYRoSkIAEpkjAmrYjMMqmfOHChZRLWfIqXxmrLgEJSEACEhgKgdE1ZRrugwcPUi7YAijzOl\/EOE6IwI0bKV29eixHRxMqzFIkIIHhEmgvs9E15VUooiHnMTRo7LnNuQQkIAEJSGCIBEbXlGmyQwRpThKQgAQkIIF9CfTZlPfKne9+Q\/Zp1I8ePUrIq1evUh8SEN68fp1eV8mb14uQWj9r3sb0kX\/fZ7558yaF9J2L5\/fzHpK73HkN8BxHFg\/MEX8ZbVOmEYfQnHe9g2vXrqXLly+nW7dupYcPHy4JvgAf\/QAAEABJREFUl\/zkyZMle1XsLrbHjx8vUn\/+\/Hl69uzZkrx88XKlnzURw1675LDtmraZbJPP559\/nl68eLGQX3\/6aWv3tElOQ+KySb5dxchl+bkCe7ksc9mHyZ07dxbPcp7pi4fmSL+Mtik3xfv27dvp3r176fr16+ncuXNLcvbs2XTmzJkle1Xswlaxxyr76dOnF6UcHh6mU6dOLcnJkydX+lkTMW3mmdfQNpP8rHXzr3\/jG4n6kQ8++KC1e1qXB\/4hcSGfoYhclp8r3I1clrnsw4RnOM\/ymzdvLp6ZY\/0yuqa8z3fFVZd08eLFdOnSpXT+\/PlEYywlml5pb1Inr4MTJ9KJCjk4OMCd6vysOXgb02ROq\/bqgsmq83PfwXvvJepHvvb++5V3mMe3OR8Slzbr3HZvuRxWvi7lssxlHyY8w3mW80xfPDRH+uW9keZdmTYfZ5dNGx175QKNmxIwTgISkIAEOiDQa1OmYVbVWGcnlgaLPxds+BDmdT78igQkIAEJSGCoBHptyrtCofHmUu6zylfGqs+YgKVLQAISGBiBXppyfCcLi5jnI00VnyIBCUhAAhKYE4FemjJNFwE0YynYFQlIYCcCLpKABEZMoJemPGJepi4BCUhAAhJojUCvTZnvkFurzI0lIIHpELASCcyEQK9NGcb5z5LzOT5FAmsJ3L2b0v37x3LlytpwAyQgAQkMmUCvTZkmzHfLVTJkaOYmAQlIoAECbiGBJQK9NuWlbDRIQAISkIAEZkzApjzjy7d0CUhAAo0TcMO9CPTalPnYmo+w96rAxRKQgAQkIIGJEOi1KUdDZixlInwHVcbRL1+kn\/7iWaX8\/FfH\/0TkoBI2GQlIQALDINBZFr02Zb5TrpPOCMzooI+fvkp18snTL2ZEwlIlIAEJDJNAr015mEjMSgISkIAEJNAPgU6bclli+ZF1rpex6hKQgAQkIIGpE+i1KVd9dA1w7IzKdgTqfmbsz4u342i0BCQggb4I9NqUq4qmIfMdc5VvmLbhZOXPi4dzF2YiAQlIYBcCg2vKuxThGglIQAISkMAUCNiUp3CLLdQwmi1v3Ejp6tVjOToaTdomKgEJSKCKQK9NmY+pq4SPsKuS1SYBCUhAAhKYMoFemzLNt0qmDNza+iTg2RKQgASGTaDXpjxsNGYnAQlIQAIS6JbAIJpy\/hF2t+V7mgTGT8AKJCCB6RDovSnTkPOPsNGng9dKJCABCUhAApsT6LUp04BpyHm66Nhzm3MJSGBuBKxXAvMk0GtTnidyq+6cwO\/9XkpVkvxPAhKQwLAI2JSHdR9m0zQBmvGbNylVCb6mz3M\/CawgoEsC6wj02pSrPqrmo2vs6xLXL4G1BGi6NOO6QHzE1Pm1S0ACEuiYQK9NmVppwDTiEHTsigQkIAEJjIGAOTZJoPemTDE04hB0RQISkIAEJDBHAoNoynMEb80dEFj38TQfXRPTQSoeIQEJjItAX9n20pTjo+q86Cpb7ncugZ0I0HRpvlWCb6dNXSQBCUigHQKdN2Wab9VH1WHD306p7jpJAnfvpnT\/\/rFcuVJdIs23SqqjtUpAAhLojUC7Tbkoi4ZL8y3M76j4iXvHqCIBCUhAAhKYAYFOm\/IMeFri0AjkH1uTW64zx6ZIQAISGAgBm\/L2F+GKsRCg6eYfW5c6Pmxjqcc8JSCByROwKU\/+ii1QAhKQgATGQqDTprzJz4v5eTJxYwE4mTwtRAISkIAEeifQaVOmWhoujZd5Kdjxl\/YqndiQ0h92xtKnPmMCfFw94\/ItXQISGD6BzpsySGi8NMxSsONfJ6wjNgQ91jAPOyN6+BxnRuDNm\/TVvw6Viv\/4WTJioy7AqEpAAn0S6KUpUzANsxTs64Qmy7qquCofsdir4rXNgABNFylLxYaUdnUJSEACPRLorSn3WPM7Rz969Cghr169Sn1IJPPm9ev0ukrevF6E1PpZsy5mnf\/LPRaHrPnSB5+2z+Q76bbPaHz\/nl6r1tHPM0Lum3HnOY6seYwN3j26psx3vjlVvgsubbl\/3fzatWvp8uXL6datW+nhw4dLwiU\/efJkyV4Vu4vt8ePHixSfP3+enj17tiQvX7xc6WfNuph1fvbg\/MVBK77wcCDftpnswnHXNX\/94EFjdzslLrvyrFonl+XnCpzkssxlHyZ37txZPMt5pq94jA3eNbqmnBPdtyGz1+3bt9O9e\/fS9evX07lz55bk7Nmz6cyZM0v2qthdbKdPnyaNdHh4mE6dOrUkJ0+erPR\/\/etfTwhr6mLwIev8xLx\/8v3FOUe\/fJF++otnS\/LzX71c5AiLtpnswnEIa+Sy\/P7hXjbk0tp7jByGKHJZfr3sw4RnOM\/ymzdvLp5lY\/0y2qbcREPm0i5evJguXbqUzp8\/v2g6NMdcaFg0tdzW9Jw8Dk6cSCcq5ODgAHcq\/Qdf2hHWMBJUxuBD1vmJOXFwgi3Sx09fVconT79Y+Km9Cyacs5H84R+mwz\/4g2P58z+vvMPDL\/+wUylf\/kFoozM2jBsUlw1zbrL+ur3kclj5upTLMpd9mPAM51nOM33xsBrpl1E25aYa8kjvzLS3IRC\/Yc0vdZWCb5u9jJXAtgSMl8CWBEbXlFc1ZH62jD9ngI49tzmfCQGaLo24rlx8xNT5tUtAAhLomMDomjJ8aLSlYEdowLkPHbsiAQlIQAKNEHCTFgmMrinTZKskZ5T7c7tzCUhAAhKQwJAJjK4pDxmmuQ2MwLqPp\/nompiBpW06EpBADwQGcqRNeSAXYRotEaDp0nyrBF9Lx7qtBCQggV0I2JR3oeaacRGg+VbJuKowWwlIYAYEGm7KMyBmiRKQgAQkIIGWCNiUWwLrthKQgAQkIIFtCdiUNyBmiAQkIAEJSKALAjblLih7hgQkIAEJSGADAjblDSBNM8SqJCABCUhgaARsykO7EfNpngD\/d6h811LPfc4lIAEJ9EjAptwjfI9ugMDduyndv38sV66kBnZ0CwlIQAK9EbAp94begyUgAQlIQALvErApv8tDTQIDJGBKEpDAXAjYlOdy09YpAQlIQAKDJ2BTHvwVmaAEpknAqiQggWUCNuVlJlokIAEJSEACvRCwKfeC3UM7JcA\/RpEfWOq5z7kE9iLgYgnsR8CmvB8\/V0tAAhKQgAQaI2BTbgylG0lAAhKYJgGr6o6ATbk71p4kAQlIQAISWEnAprwSj04JSEACEpgmgWFWZVMe5r2YlQQkIAEJzJCATXmGlz6pkm\/cSOnq1WM5Onq3NP7hiU3k3VVqEpCABHojsG9T7i1xD5bAWgL8X59yYUGVjl2RgAQkMAACNuUBXIIpdECA75hpyPlR6Nhzm3MJSEACPRKwKVfB1yYBCUhAAhLogYBNuQfoHikBCUhAAhKoImBTrqIyTdu8q4qPqvm4Ohfs8yZj9RKQwIAI2JQ7uIyjX75IP\/3Fs0r5+a9edpCBRywI0IBLWTj8IgEJSGAYBGzKHdzDx09fpTr55OkXHWQw4SNWlZZ\/R7xqvmoPfRKQgAQ6JGBT7hC2R3VMoPyuuE7vOC2Pk4AEJFBHwKZcR0a7BPoj4MkSkMBMCdiUZ3rxli0BCUhAAsMjYFMe3p2YUVsEqn6u3NZZ7rtMQIsEJLCWgE15LSIDJkGAhlz1M2XskyjQIiQggSkQsClP4RatQQIS6IuA50qgUQKjbcoXLlyoBYEvpDZIxzQI3L2b0v37x3LlyjRqsgoJSGC2BEbZlGm4dTeG78GDBykEvS5WuwQkIAEJVBDQ1BuB0TXlVU0WH804p4mOPbc5nyEBfp7Mz49LwT5DHJYsAQkMk8DomjJNFhkmTrMaNAEacCmDTtjkJCCBBgmMYqvRNeWmqT569Cghr169Sk1L5Prm9ev0uk7evF6E1cUsnF9+qfMv9l2zx+t1fnJbE\/NlCov\/Nc2oy\/1S+V3yW73LHDyr+feZTGXKa4DnOLJ4UI34y+yb8rVr19Lly5fTrVu30sOHD5eES37y5MmSvSq2tD1+\/Hjx0nj+\/Hl69uxZpbx8cfwPUtTFYGcTxl33WHcG+66L4XzyoKZ9mJSMutJPfO1r6a8fPKiUpnIYI5emal+1j1yWnyvwkssyl32Y3LlzZ\/Es55nOs2qssnVTHmuhdXnfvn073bt3L12\/fj2dO3duSc6ePZvOnDmzZK+KLW2nT59eHHt4eCqdOlUtJ0+efBtzWBnz\/sn3V\/rZd90ef+evzmGTPSIPWOzDpGTUlQ7Ets8aI5e2mbC\/XJafK3JpngnPcJ7lN2\/e5O0+Wpl9U7548WK6dOlSOn\/+fDo8PFySaFhVvk1svDIOThykEydOVMrBwQEh6aDGf+LgxGr\/l+sO1uzxd\/7986DmfZmwR9cCxLbPHCOXtpmwv1yWnytyaZ4Jz3Ce5TzTeb+PVSbVlPkFsPI3rdGxr74gvZMnwC948TPkyRdqgRKQwJgJTKopcxE0YBpxCDp2ZYYEaMK5gCDXY45dkYAEJDAAAqNtyquaLb6QATDuLYXZH8x3x5vI7EEJQAISGAqB0TbloQA0DwlIQAISkEBTBGzKTZF0n44IFMfcuJHS1avHcnRUOAuVj6vDxBwJ3VECEpDAAAjYlAdwCabQAQEaMB9lc1TM0ZljUyQgAQkMgIBNeQCXYAoSkIAEJCABCNiUoaBIQAISkIAEBkDApjyASzAFCUyTgFVJQALbErApb0vM+HESiJ8f8zNk5lSRz9EVCUhAAj0TsCn3fAFjPD7+ytDR5U4zRiLxfB42RwmsIaBbAm0SsCm3SXeie9OUP\/zww8Xf5T3REi1LAhKQQC8EbMq9YB\/voUe\/fJH+21\/+Lt37i0\/Tn\/zVZ+mnv3j2lfz8V8f\/DOVgq+Pj6hCSZM6oSGD2BAQwFAI25aHcxEjy+Pjpq\/Tot\/8v\/c1v\/jZ9\/NtXCT3kk6dfDLcKGjAfVyORJXPsoTtKQAIS6JmATbnnC\/B4CUhAAhJoh8AYd7Upj\/HWzFkCEpCABCZJwKY8yWu1KAlIQAISGCOB9U15jFWZswRKAvHz4\/gZMiOCvYxVl4AEJNATAZtyT+A9tiECd++mdP\/+sVy5Ur9pNGCacC71K\/RIQAIS6JzAXJty56A9UAISkIAEJLCOgE15HSH90yDAd8d8tzyNaqxCAhKYKAGb8pQudiC18BeM5H+pSMx7\/ctFoiEzljIQbqYhAQlIwKbsa6BxAvGXiZRjr3+5CN8p10njBNxQAhKQwG4EbMq7cXNVdwQ8SQISkMBsCNiUZ3PVFprKj61DF40EJCCBgRCwKe95EZ+9fL3nDi7vhAANeEgfX3dStIdIQAJjI2BT3vPGfv\/wvTTIX2zas642l9fx4hfCev1lsDaLdm8JSEACGxCwKW8AaV1I+QtNoff6i03rku7RH3yqRpn1eDH7H+0OEpDAngRsynsCdPlICPDRNR9hjyRd05SABOZJwKY8z3ufR9U04VyoOtdjjl2RQB0B7RLokIBNuUPYHtUCgRs3Urp69XP2FMEAAArMSURBVFiOjt49gO+ON5F3V6lJQAIS6I2ATbk39B4sAQlIoDcCHjxQAjblgV6MaUlAAhKQwPwI2JTnd+fzqZifGXdY7YkTJxLS4ZEeJQEJ5AQmMLcpT+ASLWEYBGjIH374oY15GNdhFhIYJQGb8iivzaSHRCD+MpT\/+pdP072\/+DT9yV99lvyLUIZ0Q+YigfEQqGjK40neTCWwlgAfYa+TtZusDvjqL0H57av0N7\/52\/TxlyM2\/yKU1dz0SkACywRsystMtEyJgP+XqCndprVIYPIEZtOUJ3+TFigBCUhg5AT4vYyRl7B3+pNsyhcuXEghexNyAwlIQAISaIzAqn9Zj6b89\/7+ucTY2IEj22hyTZlm\/ODBgxSCvs+dPHnyJP3oRz9KjPvs087afnb9\/DdP0v\/82Q8TYz8ZbHbq0f9+vviFK37pqpT416jil7Tq\/JuddBwFjzFwOc62u6+8d3wPLfOeK5dV\/7LeT\/7H\/0n\/+T\/9x\/To0aNlYDOxTKop04BpxvndoWPPbdvMeeP8+Mc\/3mbJ5GM\/e9uUh14ov2xVJ\/FLWOv829Q4Fi7b1NRErO+haopz5lL3vvvof\/3f9IMf\/KAa2Eysk2rKu9wZfyJbJbHnqpjfffpJ+t2vlyXW1PlZsy5mnb+PPYLJ06Lmdbmu80ctEbfJ+PLlyxRSxjfBvdyzSo9zSi5VsXO0BZc51r6q5jlzifcM7\/lcmmASe4x1nG1T\/ta3vpUuXryYrl27li5fvlwr+Llcxrq4\/\/BH\/zT98b\/+Z0uCnTWMVX5s+FbFrPP3scd\/+Xf\/EiSJkfND1uW6zs8+EQOTTeSHP\/tZ+u8ffbSQ73\/\/++\/cI3uxZ5XgY3\/GVX5i1knsAQ\/AMLIn9nVr5+DnvQMXxjnUu2mN8JgrF94bvEdK4b2zLxO48mznGc9eY5NZN+Xbt2+ne\/fuKSNm8I\/\/9E\/TP\/izP1vIH\/3kJwO+S19nvtd8DXT1GuDZPrZmHPnOtikDgD9JXbp0KSky8DXga8DXwHReAzzbecaPUWbdlMd4YeYsgSERMBcJSKBZApNqylW\/ac1vXmNvFpu7SUACEpCABJonMKmmDB4aMI04BB27IgEJSGAzAkZJoD8Ck2vKoKQRh6ArEpCABCQggTEQmGRTHgN4c5SABCTQJQHPGgcBm\/Kae4qPwRnXhE7WTe0hq4qMGMZVcVPxUWfIupo2jVu3z9D9USfjqlzxh6yKm4ovamWsqwlfLnVxU7fDYOo1rqrPpryCDi+O+BicEX1F+CRd1EztIehVhWKPGEb0qrip2KiPOkPQp1LbrnXAIHgwolfthR1\/CHpV3FRs1Be1MqKXtWHDlwu2Mm7q+vY1T4+ITbnmTnlx8AbJ3ejYc9uU59RKzXmN6NhzGzr23DbleVW91I+9qm7s+Kt8U7FV1UjN2PMa0bHntlLPfWOf19WLfey1NZ2\/TI6J2pSPOfhVAhLokAAP4JAOj\/WoARPgD2fIgFPsJDWacicHech0CZRvJB62pW261a+uTBbLfIIJrxEEfTlqXpbgAIsQbPOiYLUQsClDQWmMAA8UHybHOGVxzKH86uujJJJSvFZgE4JtOVLL1AlMtylP\/eYGWB8PER4oA0ytt5RgEkISzBkVCUhAAlUEbMpVVLRtTYBmY0N+Fxs8csGLzqhIQAISqCJgU66i8qWNhyeN5svpV\/9Dx\/6VoftJpydSKzXnh6JjX2fL\/VObUz8c8rrQsee2Oc2pHQZ5zejYcxs69tw25XlVvdSPfcp1W9vuBGzKK9jxxuENFIK+InySLmqO+hnRo1D0fI6eS\/imOMIhrxU96sQe8zmNMKD2EPSoH1vMsaOHoIdviiP1Ra2M6FEnOnNszHPBhk+ZFwGb8pr75o0RsiZ0su6onzEvMnTGKsljV85H6sxrzkvAnusxr7OHfwojNYbk9WArdWxIbp\/qnDpD8hqxhc48l7DPbYTB3GrO67Up5zScS0ACEpCABHokYFPuEb5HT5qAxUlAAhLYmoBNeWtkLpCABCQgAQm0Q8Cm3A5Xd5VAawT4ZaDWNl+3sf5REfC1MqrrWiRrU15g8IsExkGAh+zcfxFmHDc1jCx5rfCaGUY2ZrEJAZvyJpQmEsObc500USpnNLFPV3uQby5dndvUOXnu5bw8Y52f+DIm1\/EjuS2f41slxK7yN+3jvJDYGz3mb8dBDeQXMqjETKYTAjblTjAP4xD+1JwLWeU6c2xzEh5+1J0LtiEyIC\/yrMoNe5XksbE+j8OWx8Q8j8nn4WfM7TGv24\/4VT78bQh5sW+MzIcsMCLXkCZyZS\/2bWIv92ifgE25fcazO4GHQJtFN\/WAYZ+qXLHha7OGrvemHuoqz8WGr7Tvqtftxxn4dt3XdTsQcMkoCdiUR3ltJt02gSk1kL4bYt\/nt\/1acX8JNEnAptwkzYnsxUM0l7IsfNgYQ9BDsJVzbEhuz\/WwM2IPQQ\/Jbfm89OML26qRxksssioOHzG5YMtlnS9iI67UsYet6xEOXZ9Zd14dh9zOPJe6vdbZ6\/bAXq7NbTFnDMnjsYXOHCn13IYPHSnnuR5+bCFhYwwJ30THyZdlU578FW9XIG9sHtK5YCt3wbYuJtbkscyRWMs84hjRw8eIjh1BR8o5OnH4QtCxr5M8njVIuQZbxMWILeKYh50RPXwxYkPwI9hzHRs69l2EtaXssk+sKfdCD9+qkThqiZhSD3vdyFrW5H507Nhijh6CDd82wppYz4je5Hr2Q9gbYe9cx4aOHUFHyjkx2EPQickFW\/gZ0XO\/83ERsCmP675az5Y39SaHbBrHXmVsqROD8DApfejY8dcJfuJyPzr23LZqTnxIvo459nJt2Kr8+LDna7AhYcOf69jRsTPfVlhbyrZ75PHlXui5nzm5lpLH4ct11rQhu5yxy5o893Xr8SOxpooFfuwRU474iMnt6NhLW647HyGBLGWbcgbD6TEB3vS5HFun95Uaq6qqevBVxe1r4\/xS9t1zl\/XksMs6OJVS7sPeueBHZ9xFOI\/1ueyyz1jW5HXGfCy5m+duBGzKu3Gb7Cre+Dz4cplssT0XljPO53VpEcP91Pnr7Luuq9tvUzvnlsJabIx1gj\/qZETPY9FzISb3T2me15nPt6kRPqzdZo2x\/RGYUFPuD+KUT+YN3VV9PDjK89Cxr8oBP3F5DDr23FbO8RNX2rHhw86IzjyXsFX58WHP48s5fuJKe5WtjNlFX3Uevl327GNNU3y22acqtsq2igeMyzXo2OvW4SOm9FfZyhj18RKwKY\/37lrJPB4EvPER9FYOqtmU8zg3BL0MxRb+8OU2fOjhWzUSR3wu2PI16LmfObaIYY4tBD18q0biYk2M2DZZU8bE+nLM49i79GPLY4YyJy9yZcxzQseeC7Y8Jp8Thx4jc4Q12ELQsSPMw86Ijj0XbPhC0HN\/1ZyYiGdEr4rLbcQQmwu2PGbVnHXbxK\/aS183BGzK3XDe+ZQ2F9a9WbGHcD5zxpBSx57b6uZlXJUeNvZA0KsEH5L70ENy+7p5rImxKj58MZYxYWes8pW20InPJeyrRuJ52EYMep1ETIxlXNjzkZhcr5pvErPvurozsOdSdU7YquKw4WcMQc8l7IzYY2Qegi0kbIzYGKsEX0idv7RHfIy5H1uuMw8br5GYY1fGQcCmPI57MksJvEPAh+07OFQqCPgaqYAyAtP\/BwAA\/\/+DkuB9AAAABklEQVQDAMDKvIIX0ZVGAAAAAElFTkSuQmCC","height":234,"width":388}}
%---
%[output:0f35e493]
%   data: {"dataType":"text","outputData":{"text":"[23:23:15][INFO]  GPU free memory: 4.5 GB\n","truncated":false}}
%---
%[output:5035f545]
%   data: {"dataType":"text","outputData":{"text":"[23:23:15][INFO]  Safe batch size: 441907 molecules (3620.1 MB\/batch)\n","truncated":false}}
%---
%[output:5f3193cf]
%   data: {"dataType":"text","outputData":{"text":"[23:23:15][INFO]  Number of batches needed for screening 1 million: 3\n","truncated":false}}
%---
%[output:56d9677d]
%   data: {"dataType":"text","outputData":{"text":"[23:23:15][INFO]  Number of batches needed for screening 10 million: 23\n","truncated":false}}
%---
%[output:76ce6c50]
%   data: {"dataType":"text","outputData":{"text":"[23:23:15][INFO]  R01: Complete -- GPU-accelerated similarity screening finished\n","truncated":false}}
%---
