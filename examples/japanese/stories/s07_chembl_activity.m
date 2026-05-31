%[text] # S07: ChEMBL から生物活性データを取得する
%[text] EasyMolKit 応用ストーリー -- レイヤー 2
%[text] 
%[text] 薬の「効くか効かないか」を示すデータが数十万件も無料で公開されているとしたら？
%[text] ChEMBL は European Bioinformatics Institute が管理する生物活性データベースで、査読済み学術誌から抽出した IC50 測定値を REST API で提供しています。EGFR 阻害剤だけでも数千件のバイオアクティビティデータが収録されており、このスクリプトでは ChEMBL からデータを取得し、高活性化合物を絞り込み、構造類似度で既存薬エルロチニブとの関係を探るワークフローを体験します。
%[text] ### ストーリー
%[text] あなたは計算創薬ラボの大学院生です。
%[text] 指導教員から、公開バイオアクティビティデータベース ChEMBL を使ってEGFR（上皮成長因子受容体）の既知阻害剤の全体像を調査するよう依頼されました。EGFR は非小細胞肺がん（NSCLC）などの固形腫瘍で過剰発現する受容体チロシンキナーゼです。2004 年以降、複数の小分子 EGFR 阻害剤がFDA から承認されています:
%[text] - ゲフィチニブ   (Iressa,   2003, 第 1 世代)
%[text] - エルロチニブ   (Tarceva,  2004, 第 1 世代)
%[text] - アファチニブ   (Gilotrif, 2013, 第 2 世代, 共有結合)
%[text] - オシメルチニブ (Tagrisso, 2015, 第 3 世代, T790M 選択的) \
%[text] ### あなたの課題:
%[text] 1. ChEMBL で EGFR を特定し、ターゲット識別子を確認します。
%[text] 2. EGFR に対する IC50 バイオアクティビティデータを 50 件取得します。
%[text] 3. IC50 分布を調べ、効果量の単位を理解します。
%[text] 4. 高活性阻害剤（$&dollar&;\\text{IC}\_{50} \\le 100\\text{ nM}&dollar&;$）にフィルタします。
%[text] 5. 活性化合物セットを SDF ファイルに保存します。
%[text] 6. エルロチニブに対する類似度スクリーニングを行い活性化合物を順位付けます。
%[text] 7. 上位ヒットを可視化し、別のクエリ化合物と比較します。 \
%[text] ### 学習目標
%[text] - emk.db.searchChemblTarget で ChEMBL のタンパク質ターゲットを検索できます
%[text] - emk.db.getChemblActivity で IC50 バイオアクティビティデータをダウンロードできます
%[text] - 活性単位（nM）と創薬における効果量の閾値を理解できます
%[text] - バイオアクティビティデータと構造類似度解析を連携させられます
%[text] - emk.io.writeSdf でキュレーション済み化合物セットをエクスポートできます \
%[text] ### 前提条件
%[text] - F03（フィンガープリント）と F04（類似度）の完了
%[text] - 推奨: S04（バーチャルスクリーニング）と S06（PubChem 検索）
%[text] - RDKitインストール済み（`emk.setup.install()` を一度だけ実行しておく）
%[text] - 追加Toolbox不要（MATLAB だけで動きます）
%[text] - インターネット接続が必要（ChEMBL REST API） \
%[text] **所要時間**: 35〜50 分 | 実行方法: Ctrl+Enter でセクションを一つずつ実行
%[text] 
%[text] ## データ
%[text] - ローカルデータファイル不要。化合物はリアルタイムで以下から取得:
%[text] - ChEMBL REST API  https://www.ebi.ac.uk/chembl/api/data/ \
%[text] ## 参照文献
%[text] - Mendez D et al. (2019) ChEMBL: towards direct deposition of bioassay data. *Nucleic Acids Res* 47:D930-D940. doi:10.1093/nar/gky1075 〔Open Access〕
%[text] - Paez JG et al. (2004) EGFR mutations in lung cancer: correlation with clinical response to gefitinib therapy. *Science* 304:1497-1500. doi:10.1126/science.1099314 〔要機関アクセス〕
%[text] - Yun CH et al. (2008) The T790M mutation in EGFR kinase causes drug resistance by increasing the affinity for ATP. *Proc Natl Acad Sci USA* 105:2070-2075. doi:10.1073/pnas.0709662105 〔Open Access〕
%[text] - Willett P (2006) Similarity-based virtual screening using 2D fingerprints. *Drug Discov Today* 11:1046-1053. doi:10.1016/j.drudis.2006.10.005 〔要機関アクセス〕\
%[text] ## 注意
%[text] このストーリーはライブ ChEMBL REST API にクエリを送信します。インターネット接続が必要です。ChEMBL のデータベース更新に伴い結果がわずかに異なる場合があります。
%%
%[text] ## セクション 0: セットアップ
%[text] パスと Python 環境を初期化します。
%[text] **常にこのセクションを最初に実行してください。**
% Resolve project root (works for Desktop, MCP, and MATLAB Online)
sDir = fileparts(mfilename('fullpath'));
if strlength(sDir) > 0
    addpath(genpath(fullfile(sDir, '..', '..', '..', 'src')));
elseif ~isempty(which("logInfo"))
    addpath(genpath(fileparts(fileparts(which("logInfo")))));
end
projectRoot = resolveProjectRoot();
addpath(genpath(fullfile(projectRoot, 'src')));
emk.setup.initPython();
%[text] Python/RDKit プロセスをウォームアップします（最初の呼び出しは少し時間がかかります）。
mol_warmup = emk.mol.fromSmiles("C");   % メタン -- 軽量
clear mol_warmup;
logInfo("S07: セットアップ完了");
%%
%[text] ## セクション 1: ターゲット -- 腫瘍学の創薬における EGFR
%[text] ### タンパク質キナーゼターゲットと EGFR の成功事例
%[text] EGFR（上皮成長因子受容体 / 遺伝子名: EGFR）は、細胞の増殖や生存のスイッチを握る「膜貫通型受容体チロシンキナーゼ」です。がん治療、特に非小細胞肺がん（NSCLC）の創薬において、最も劇的な成功を収めた標的タンパク質のひとつとして知られています。
%[text] 【EGFR阻害剤の進化と世代交代の歴史】
%[text] がん細胞の遺伝子に変異（エクソン19欠失やL858R変異）が起きると、このスイッチが「常時ON」に固定され、細胞が暴走して腫瘍が拡大します。小分子阻害剤は、この活性部位にある「ATP結合ポケット」に嵌まることでシグナルを遮断します。
%[text] 第 1 世代（可逆的阻害剤） ： エルロチニブ、ゲフィチニブ
%[text] - 変異型EGFRに劇的に効くものの、長期間使用するとポケットの入り口の原子が変わる「T790M（ゲートキーパー変異）」が起き、薬が阻害できなくなる（薬剤耐性）という課題に直面しました。 \
%[text] 第 2 世代（不可逆的共有結合阻害剤） ： アファチニブ、ダコミチニブ
%[text] - ポケット内のアミノ酸残基（Cys797）と強力な「共有結合」を形成することで、T790M変異による耐性を部分的に克服しました。 \
%[text] 第 3 世代（変異選択的阻害剤） ： オシメルチニブ
%[text] - 耐性の原因である「T790M変異型」をピンポイントで強力に阻害するよう専用設計され、正常な野生型EGFRへの副作用を抑えることに成功しました。現在は第一選択の標準治療薬となっています。 \
%[text] 【ChEMBLデータベースの意義】
%[text] これら歴代の阻害剤開発の裏には、世界中の研究者が蓄積してきた膨大な「化合物と活性のデータ」があります。それを手動で丁寧に集約（キュレーション）したのが、欧州バイオインフォマティクス研究所（EMBL-EBI）の「ChEMBL（ケムブル）」です。
%[text] 2024年時点で、240万件を超える化合物、2000万件以上のバイオアクティビティ、15,000件のターゲットタンパク質が登録されており、今回のターゲットであるヒトEGFR（CHEMBL203）にも、過去の実験論文から数千件もの貴重な活性データが紐づいています。
%[text] 【生物活性の指標：IC50 とは】
%[text] アッセイ（試験）で標的タンパク質の働きを「50%阻害（半分にブロック）するために必要な化合物の濃度」を指します。単位は「nM（ナノモラー：$&dollar&;10^{-9}&dollar&;$ M）」で表され、数値が小さい（より薄い濃度で効く）ほど、阻害剤として強力であることを意味します。一般的な創薬研究における活性の目安は以下の通りです：
%[text] - $&dollar&;\\text{IC}\_{50} \\le 100\\text{ nM}&dollar&;$  ： 高活性（次のステップである「リード化合物」へ進める優秀な品質）
%[text] - $&dollar&;100\\text{ nM} \< \\text{IC}\_{50} \\le 1000\\text{ nM}&dollar&;$ ： 中程度の活性（構造の最適化・改造が必要なレベル）
%[text] - $&dollar&;\\text{IC}\_{50} \> 1000\\text{ nM}&dollar&;$  ： 弱活性（1 $&dollar&;\\mu\\text{M}&dollar&;$ の壁を超えられず、不採用となることが多いライン）
%[text] - $&dollar&;\\text{IC}\_{50} \> 10000\\text{ nM}&dollar&;$ ： 不活性（10 $&dollar&;\\mu\\text{M}&dollar&;$ 以上必要であり、薬としての効果はないと判定） \
EGFR_NAME   = "EGFR";
ERLOTINIB_SMILES = "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1";
ERLOTINIB_NAME   = "Erlotinib";

logInfo("ストーリーターゲット: %s", EGFR_NAME);
logInfo("参照化合物: %s", ERLOTINIB_NAME);
logInfo("  SMILES: %s", ERLOTINIB_SMILES);
%[text] **✏️ やってみよう 1**
%[text] Q: IC50 とは何の略で、単位は？
%[text]    IC50 が低いほど薬物の効果量として優れているのはなぜか？
%[text]    期待値: IC50 = 半数阻害濃度; 単位は nM。
%[text]    IC50 が低いと、より低い濃度で 50% 阻害を達成できる（= より強力）。
%[text] Q: エルロチニブは第 1 世代 EGFR 阻害剤です。SMILES を見て以下を特定してみましょう:
%[text]      (a) アルキン基（-C\#C-）-- EGFR バックポケット結合のファーマコフォア。
%[text]      (b) キナゾリンコア（窒素 2 個の二環系）-- ヒンジ領域に結合するヒンジバインディングモチーフ。
%[text]      (c) メトキシエトキシチェーン 2 本（-OCCOC-）-- 水溶性と選択性の決定基。
%[text]    ヒント: mol = `emk.mol.fromSmiles(ERLOTINIB_SMILES)`; `emk.viz.draw2d(mol)`
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 2: ChEMBL で EGFR を検索する
%[text] ### ChEMBL ターゲット識別子
%[text] ChEMBL に登録されている膨大なタンパク質には、それぞれ `CHEMBL{数字}` という形式の「ChEMBL ターゲット ID（一意識別子）」が割り当てられています（例：ヒトEGFRなら `CHEMBL203`、ヒトCOX-2なら `CHEMBL230`）。このIDは、そのタンパク質に関するあらゆる論文データや化合物データを引き出すための不変のマスターキーとなります。
%[text] 通常ならWebサイトで検索して探すこのIDを、EasyMolKitでは `emk.db.searchChemblTarget` 関数を使い、タンパク質の正式名（Preferred Name）をキーワードにしてプログラムから直接問い合わせ・取得することができます。まずは、私たちが解析したいヒト（Homo sapiens）のEGFRに紐づく正確なターゲットIDを特定してみましょう。
logInfo("ChEMBL で EGFR ターゲットを検索中...");
targetTbl = emk.db.searchChemblTarget( ...
    "Epidermal growth factor receptor", ...
    TargetType="SINGLE PROTEIN", MaxRows=5);

logInfo("ChEMBL ターゲット検索結果:");
disp(targetTbl(:, ["TargetChEMBLID","PreferredName","Organism"]));
%[text] CHEMBL203 がヒト EGFR であることを確認します。Organism == "Homo sapiens" の行が目的のターゲットです。
humanIdx = find(targetTbl.Organism == "Homo sapiens", 1, "first");
if isempty(humanIdx)
    % フォールバック: 既知の ChEMBL ID を直接使用する。
    logWarn("結果にヒト EGFR が見つからなかった; CHEMBL203 を直接使用");
    EGFR_CHEMBL_ID = "CHEMBL203";
else
    EGFR_CHEMBL_ID = targetTbl.TargetChEMBLID(humanIdx);
    logInfo("EGFR ターゲット確認: %s (%s, %s)", ...
        EGFR_CHEMBL_ID, ...
        targetTbl.PreferredName(humanIdx), ...
        targetTbl.Organism(humanIdx));
end
%[text] **✏️ やってみよう 2**
%[text] Q: EGFR の ChEMBL ターゲット ID は？
%[text]    期待値: "CHEMBL203"
%[text] Q: 検索結果に複数の EGFR エントリが返された。なぜか？
%[text]    期待値: ChEMBL は複数の生物種（ヒト・マウス・ラットなど）のターゲットを
%[text]    収録している。ヒト（Homo sapiens）のエントリが最多のバイオアクティビティ
%[text]    データを持ち、治療的に重要なエントリ。
%[text] Q: 別のターゲットを検索してみよう。
%[text]    ヒト COX-2（シクロオキシゲナーゼ-2）の ChEMBL ID は？
%[text]    ヒント: `emk.db.searchChemblTarget("cyclooxygenase-2", MaxRows=10)`
%[text]    期待値: "CHEMBL230"
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 3: EGFR の IC50 バイオアクティビティデータをダウンロードする
%[text] ### ChEMBL アクティビティエンドポイント
%[text] ターゲットのID（マスターキー）が特定できたら、次はそのタンパク質に対して過去に行われた試験の「測定データ（バイオアクティビティレコード）」をダウンロードします。1つのレコードには、試験された化合物情報（ChEMBL化合物IDやSMILES構造式）に加え、実験条件、測定値、単位、そして数値の関係性を示す「関係演算子（Relation）」が含まれています。
%[text] 【関係演算子（Relation）の種類と意味】
%[text] - `'='` ： 活性が定量的に測定された正確なデータ（例：$&dollar&;\\text{IC}\_{50} = 5.2\\text{ nM}&dollar&;$）
%[text] - `'\<'` ： スクリーニングの上限を超えて強く効いたデータ（例：$&dollar&;\\text{IC}\_{50} \< 10\\text{ nM}&dollar&;$）
%[text] - `'\>'` ： 活性が弱すぎて、実験の最大濃度でも効果が測れなかったデータ（例：$&dollar&;\\text{IC}\_{50} \> 10000\\text{ nM}&dollar&;$） \
%[text] これらが混ざったままだと、このあとの統計解析やグラフプロットが正しく計算できません。そのため、EasyMolKitの `emk.db.getChemblActivity` 関数は、自動的に最も信頼性の高い `'='`（正確な等号測定値）のみを厳選し、かつ単位を `nM` に統一、RDKitで扱える綺麗な構造（正規化SMILES）を持つ行だけをフィルタリングした状態でMATLABテーブルとして返してくれます。
%[text] 現在、ChEMBLにはEGFRのデータが5,000件以上登録されていますが、本演習ではAPI通信をスムーズに行うため、件数を「最大50件（`MaxRows=50`）」に制限してリアルタイムにデータを取得してみましょう。
logInfo("ChEMBL から EGFR IC50 データを取得中（最大 50 件）...");
actTbl = emk.db.getChemblActivity(EGFR_CHEMBL_ID, ...
    ActivityType="IC50", MaxRows=50);

% Fallback: use ChEMBL ID as display name when preferred name is absent
emptyName = strlength(actTbl.Name) == 0 | actTbl.Name == "";
actTbl.Name(emptyName) = actTbl.MoleculeChEMBLID(emptyName);

logInfo("活性データセット: %d 化合物の IC50 データ", height(actTbl));
logInfo("カラム: %s", strjoin(actTbl.Properties.VariableNames, ", "));
logInfo("IC50 統計（nM）:");
logInfo("  最小値  : %.1f nM", min(actTbl.Value_nM));
logInfo("  中央値  : %.1f nM", median(actTbl.Value_nM));
logInfo("  最大値  : %.1f nM", max(actTbl.Value_nM));
%[text] 最初の数行を表示する
disp(actTbl(1:min(5, height(actTbl)), ...
    ["MoleculeChEMBLID","Name","Value_nM"]));
%[text] **データ注意**: 同一化合物（MoleculeChEMBLID が同じ）が複数行返されることがある。
%[text] ChEMBL の各行は独立したアッセイ測定を表しており、同じ化合物が異なるアッセイ条件
%[text] （細胞系・生化学的・異なる実験室など）で複数回測定されている場合に重複が生じる。
%[text] これは ChEMBL のデータ特性であり、実際の解析では `unique` や `groupsummary` で化合物ごとに集約することが多い。
%[text] **✏️ やってみよう 3**
%[text] Q: IC50 \< 100 nM（高活性）の化合物は何件か？
%[text]    ヒント: sum(actTbl.Value\_nM \<= 100)
%[text] Q: データセット内で最も活性の高い化合物は？
%[text]    ヒント: actTbl(actTbl.Value\_nM == min(actTbl.Value\_nM), :)
%[text] Q: IC50 の代わりに Ki 値を取得してみよう。Ki 測定値は存在するか？
%[text]    ヒント: kiTbl = `emk.db.getChemblActivity(EGFR_CHEMBL_ID, ActivityType="Ki", MaxRows=25)`
%[text]    注: Ki（阻害定数）は結合親和性の熱力学的指標; IC50 は機能的指標。
%[text]    両者とも単位は nM。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 4: IC50 分布を調べる
%[text] ### 効果量分布と対数スケール思考
%[text] 先ほど取得した $&dollar&;\\text{IC}\_{50}&dollar&;$ の統計値（最小値と最大値）を見ると、データが「$&dollar&;0.5\\text{ nM}&dollar&;$（超強力）」から「$&dollar&;3,000,000\\text{ nM}&dollar&;$（ほぼただの水）」まで、数百万倍という膨大な桁数にまたがって分布していることがわかります。このようにスケールの幅が広すぎるデータを生の数値（$&dollar&;\\text{nM}&dollar&;$）のままグラフ化しようとすると、最も詳しく見たいはずの「1〜100 $&dollar&;\\text{nM}&dollar&;$ 付近の高活性な化合物たちの差異」がすべてゼロ付近の壁に押しつぶされて判別できなくなってしまいます。
%[text] そこで医薬化学者は、$&dollar&;\\text{IC}\_{50}&dollar&;$ の値をモル濃度（$&dollar&;\\text{M}&dollar&;$）に変換した上で、負の常用対数をとった「$&dollar&;\\text{pIC}\_{50}&dollar&;$」という指標を標準的に使用します。 $&dollar&;\\text{pIC}\_{50} = -\\log\_{10}(\\text{IC}\_{50}\\text{ モル濃度}) = 9 - \\log\_{10}(\\text{Value\\\_nM})&dollar&;$
%[text] 【$&dollar&;\\text{pIC}\_{50}&dollar&;$ スケールのメリットと直感的な対応】
%[text]  「数値が大きいほど強力」になる：生の $&dollar&;\\text{IC}\_{50}&dollar&;$ は数字が小さいほど強いという直感に反する性質がありますが、$&dollar&;\\text{pIC}\_{50}&dollar&;$ は大きくなるほど強力です。
%[text]  桁数を等間隔（線形）に扱える：10倍の活性差が「1」の差としてきれいにプロットされ、統計解析やヒストグラムでの視覚化が格段に容易になります。
%[text] - $&dollar&;\\text{pIC}\_{50} \> 7&dollar&;$ （$&dollar&;\\text{IC}\_{50} \< 100\\text{ nM}&dollar&;$） ： 優れた活性を持つリード化合物
%[text] - $&dollar&;\\text{pIC}\_{50} \> 8&dollar&;$ （$&dollar&;\\text{IC}\_{50} \< 10\\text{ nM}&dollar&;$）  ： 既存の医薬品に匹敵する高活性
%[text] - $&dollar&;\\text{pIC}\_{50} \> 9&dollar&;$ （$&dollar&;\\text{IC}\_{50} \< 1\\text{ nM}&dollar&;$）  ： 共有結合阻害剤などにみられる超強力な分子 \
%[text] MATLAB のベクトル演算を活かして、ダウンロードしたテーブルに `pIC50` の列を瞬時に追加し、その効果量の分布を美しいヒストグラムで可視化してみましょう。
pIC50 = 9 - log10(actTbl.Value_nM);
actTbl.pIC50 = pIC50;

logInfo("pIC50 統計:");
logInfo("  最小値  : %.2f (IC50 = %.1f nM)", min(pIC50), max(actTbl.Value_nM));
logInfo("  中央値  : %.2f (IC50 = %.1f nM)", median(pIC50), median(actTbl.Value_nM));
logInfo("  最大値  : %.2f (IC50 = %.1f nM)", max(pIC50), min(actTbl.Value_nM));
%[text] 効果量クラス別の化合物数
nHigh    = sum(actTbl.Value_nM <= 100);
nMod     = sum(actTbl.Value_nM > 100  & actTbl.Value_nM <= 1000);
nWeak    = sum(actTbl.Value_nM > 1000);

logInfo("効果量分布:");
logInfo("  高活性         (IC50 <= 100 nM)    : %d 件", nHigh);
logInfo("  中程度         (100 < IC50 <= 1000): %d 件", nMod);
logInfo("  弱活性/不活性  (IC50 > 1000 nM)   : %d 件", nWeak);
%[text] pIC50 値のヒストグラム（Base MATLAB のみ使用）
figure("Name", "EGFR IC50 分布", "Position", [100 100 580 420]);
histogram(pIC50, 10, "FaceColor", [0.2 0.5 0.8]);
xlabel("pIC50  (-log_{10}[IC50 in M])");
ylabel("化合物数");
title(sprintf("EGFR IC50 分布  (n=%d, ChEMBL より)", height(actTbl)));
xline(7, "r--", "LineWidth", 1.5, "Label", "pIC50=7 (100 nM)");
xline(8, "g--", "LineWidth", 1.5, "Label", "pIC50=8 (10 nM)");
grid on;
%[text] **✏️ やってみよう 4**
%[text] Q: 野生型 EGFR に対するエルロチニブの IC50 ~2 nM に対応する pIC50 は？
%[text]    ヒント: pIC50 = 9 - log10(IC50\_nM) の式を使ってみましょう。
%[text] Q: 分布は偏っているか？偏っているとすれば、どちらの方向に、なぜか？
%[text]    期待値: 通常は右偏り（超強力な化合物は少なく、中程度の活性化合物が多い）。
%[text]    ChEMBL の IC50 データセットは同じシリーズの活性・低活性化合物を
%[text]    両方報告するため、対数正規分布を示すことが多い。
%[text] Q: ヒストグラムは横軸に pIC50 を表示している。
%[text]    生の IC50 値（nM 単位）を線形スケールの棒グラフにすると
%[text]    情報量は増えるか減るか？
%[text]    なぜ医薬化学者は対数スケールを好むのか？
%[text]    期待値: 線形スケールでは全値が低い方に圧縮される;
%[text]    pIC50（対数スケール）は分布を広げ、強力な化合物間の差異を視覚化できる。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 5: 高活性化合物をフィルタする（$&dollar&;\\text{IC}\_{50} \\le 100\\text{ nM}&dollar&;$）
%[text] ### 活性カットオフの設定
%[text] 集まった分布データの中には、薬としての効果が薄い「不活性な化合物」も多く含まれています。プロジェクトを効率的に進めるためには、リソースを注ぎ込むべき有望な分子だけを絞り込むための境界線である「活性カットオフ（閾値）」を設定する必要があります。
%[text] 【創薬ステージにおける一般的なカットオフ基準】
%[text] - 一次スクリーニング（数万件のふるい落とし） ： $&dollar&;\\text{IC}\_{50} \< 10,000\\text{ nM}&dollar&;$ ($&dollar&;10\\text{ }\\mu\\text{M}&dollar&;$)
%[text] - ヒット・トゥ・リード（本格的な検討の出発点） ： $&dollar&;\\text{IC}\_{50} \< 1,000\\text{ nM}&dollar&;$ ($&dollar&;1\\text{ }\\mu\\text{M}&dollar&;$)
%[text] - リード最適化（薬としての洗練を目指す目標）  ： $&dollar&;\\text{IC}\_{50} \\le 100\\text{ nM}&dollar&;$
%[text] - 臨床開発候補（実際の薬として勝負できるレベル）： $&dollar&;\\text{IC}\_{50} \< 10\\text{ nM}&dollar&;$ \
%[text] 本セクションでは、リード最適化の合格基準である「$&dollar&;\\text{IC}\_{50} \\le 100\\text{ nM}&dollar&;$」をカットオフ基準に設定します。MATLAB の論理インデックス機能を用いて、50件の雑多なデータからこの厳しいエリート基準をクリアした高品質な化合物だけをスパッと抽出（フィルタリング）し、さらに活性が高い順（最も強力なものを先頭）に並び替えてみましょう。
IC50_CUTOFF_NM = 100;

activeMask = actTbl.Value_nM <= IC50_CUTOFF_NM;
activeTbl  = actTbl(activeMask, :);

logInfo("活性カットオフ: IC50 <= %.0f nM", IC50_CUTOFF_NM);
logInfo("活性化合物: %d / %d 件", height(activeTbl), height(actTbl));

%[text] IC50 昇順（最も強力なものを先頭）でソートする
activeTbl = sortrows(activeTbl, "Value_nM", "ascend");
logInfo("EGFR 阻害剤 上位 5 件（IC50 順）:");
disp(activeTbl(1:min(5, height(activeTbl)), ...
    ["MoleculeChEMBLID","Name","Value_nM","pIC50"]));

%[text] **✏️ やってみよう 5**
%[text] Q: 100 nM カットオフ適用後に残る活性化合物数は？
%[text]    非常に少ない場合は 1000 nM に緩和してみよう。
%[text] Q: 最も活性の高い化合物の IC50 は？
%[text]    PubChem や ChEMBL でその構造と名前を調べられるか？
%[text]     ヒント: disp(activeTbl(1, ["MoleculeChEMBLID","Name","Value_nM","SMILES"]))
%[text] Q: IC50\_CUTOFF\_NM を 1000 に変更してこのセクションを再実行する。
%[text]    何件の化合物が該当するか？
%[text]    フラグメントベースまたは HTS 由来の EGFR データセットとして
%[text]    その数は妥当か？
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 6: 活性化合物を SDF に保存する
%[text] ### SDF（構造データファイル）形式
%[text] せっかく綺麗にキュレーション（選別）した高活性な化合物のリストも、MATLABの画面上だけで終わってしまっては他の研究者と共有できません。そこで、化学の業界標準フォーマットである「SDF（Structure-Data File：構造データファイル）」形式としてローカルに書き出します。
%[text] SDF ファイルは、単なるテキストファイルですが、内部は以下の2つの重要な情報が1つのレコードとしてセットになっています：
%[text] 1. MOL ブロック：原子の種類、それらがどう結合しているか、2次元や3次元の座標情報（結合テーブル）
%[text] 2. データフィールド：その化合物に紐づくテキストメタデータ（例：ChEMBL ID、$&dollar&;\\text{IC}\_{50}&dollar&;$ の値、アッセイ条件など） \
%[text] このファイル形式は非常に汎用性が高く、ChemDraw（構造式描画）、Schrodinger MaestroやMOE（分子モデリング）、KNIME（データフロー構築）など、世界中のあらゆる化学ソフトウェアでそのまま読み込むことができます。「データ解析担当（計算化学者）が魅力的なヒット化合物を厳選し、SDFファイルにエクスポートして実験担当（医薬化学者）へ合成や発注を依頼する」という、実際の製薬現場で毎日行われている標準的な業務引き継ぎのワークフローを、`emk.io.writeSdf` 関数を使って体験してみましょう。
if height(activeTbl) == 0
    logWarn("保存する活性化合物がない -- SDF エクスポートをスキップ");
else
    % 活性化合物の SMILES を RDKit Mol オブジェクトに変換する。
    logInfo("%d 件の活性化合物を変換中...", height(activeTbl));
    validMols = {};
    validIdx  = [];

    for i = 1:height(activeTbl)
        smi = activeTbl.SMILES(i);
        if ~emk.mol.isValid(smi)
            logWarn("  [%d] 無効な SMILES のためスキップ: %s", i, smi);
            continue;
        end
        mol = emk.mol.fromSmiles(smi);
        validMols{end+1} = mol; %#ok<AGROW>
        validIdx(end+1)  = i;   %#ok<AGROW>
    end

    logInfo("変換成功: %d / %d 件", numel(validMols), height(activeTbl));

    if numel(validMols) > 0
        runDir  = makeRunDir("Prefix", "s07_egfr");
        sdfPath = fullfile(runDir, "egfr_actives.sdf");
        emk.io.writeSdf(validMols, sdfPath);
        logInfo("活性化合物 %d 件を保存: %s", numel(validMols), sdfPath);

        % メタデータを CSV として SDF の隣に保存する
        csvPath = fullfile(runDir, "egfr_actives_metadata.csv");
        writetable(activeTbl(validIdx, :), csvPath);
        logInfo("メタデータ CSV を保存: %s", csvPath);
    end
end

%[text] **✏️ やってみよう 6**
%[text] Q: SDF ファイルをテキストエディタで開いて、以下を確認してみましょう:
%[text]    (a) MOL ブロック（"M  END" で終わる）
%[text]    (b) レコード区切り（"$$$$"）
%[text]    (c) ファイルに座標が含まれているか？
%[text]        RDKit は SDF ファイル書き込み時に 2D 座標を自動生成する。
%[text] Q: SDF を MATLAB に再読み込みするにはどうするか？
%[text]    ヒント: reloaded = `emk.io.readSdf(sdfPath)`
%[text]    （お使いのバージョンの EasyMolKit に `emk.io.readSdf` が存在する場合）
%[text] Q: SDF ファイル内のレコード数を数えてみましょう:
%[text]    ヒント: txt = fileread(sdfPath); sum(contains(strsplit(txt, newline), "$$$$"))
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 7: エルロチニブに対する類似度スクリーニング
%[text] ### 類似度ガイドによるヒット選択
%[text] 厳選された高活性阻害剤のセットを手に入れた私たちが次に立てるべき問いは、「これらの活性化合物の中に、すでに臨床で成功している既存薬（エルロチニブ）と似た形をしたものは含まれているか？」という点です。
%[text] すでに承認されている医薬品と構造が似ている（骨格を共有している）化合物を見つけることには、大きなメリットがあります：
%[text] - 安全性の予測 ： 重篤な毒性や予期せぬ副反応（オフターゲット毒性）を起こしにくい骨格である可能性が高い
%[text] - 知見の流用 ： 過去の膨大な研究データ（構造活性相関：SAR）をヒントにできるため、開発のスピードが上がる
%[text] - 合成の容易さ ： 既存の製造・合成ルートを応用（アナログ展開）しやすい \
%[text] コンピューターにこの「構造の似ている度合い」を計算させるため、各原子の周囲2結合分の局所的な化学環境を2048ビットのデジタル指紋としてスキャンする「ECFP4（モルガン半径=2）フィンガープリント」を計算します。そして、それらの指紋の重なり具合を `0.0`（全く別物）〜 `1.0`（指紋の一致）の範囲で点数化する「タニモト類似度（Tanimoto Metric）」を算出します。
%[text] `emk.similarity.rankBy` 関数を用いて、基準となる「エルロチニブ」に対する類似度スコアが高い順にChEMBLの活性化合物をソートし、ランキングテーブルを作成してみましょう。上位にランクインした化合物がどのような2次元構造を持っているか、`emk.viz.draw2d` で視覚的に確かめることで、リガンドベース創薬の強力なアプローチを体感します。 
if ~exist("validMols", "var") || numel(validMols) == 0
    logWarn("有効な活性 mol が存在しない -- セクション 6 を先に実行してください");
else
    logInfo("%s に対する類似度スクリーニングを実行中...", ERLOTINIB_NAME);

    % クエリフィンガープリント（エルロチニブ）を計算する
    mol_erl = emk.mol.fromSmiles(ERLOTINIB_SMILES);
    fp_erl  = emk.fingerprint.morgan(mol_erl, Radius=2, NBits=2048);

    % 全活性化合物のフィンガープリントを計算する
    logInfo("%d 件の活性化合物の ECFP4 フィンガープリントを計算中...", numel(validMols));
    dbFps = cell(1, numel(validMols));
    for i = 1:numel(validMols)
        dbFps{i} = emk.fingerprint.morgan(validMols{i}, Radius=2, NBits=2048);
        logProgress(i, numel(validMols), "フィンガープリント計算");
    end

    % エルロチニブに対するタニモト類似度で順位付けする
    rankResult = emk.similarity.rankBy(fp_erl, dbFps, Inf, Metric="tanimoto");

    % 結果テーブルを構築する
    rankedIdx    = rankResult.Indices;
    rankedScores = rankResult.Scores;
    nRanked      = numel(rankedIdx);

    metaTbl = activeTbl(validIdx, :);
    rankTbl = table( ...
        (1:nRanked)', ...
        metaTbl.MoleculeChEMBLID(rankedIdx(:)), ...
        metaTbl.Name(rankedIdx(:)), ...
        metaTbl.Value_nM(rankedIdx(:)), ...
        rankedScores(:), ...
        VariableNames=["Rank","ChEMBLID","Name","IC50_nM","Tanimoto_Erlotinib"]);

    logInfo("エルロチニブに最も類似する活性化合物 上位 10 件（タニモト, ECFP4）:");
    disp(rankTbl(1:min(10, nRanked), :));
%[text] **注**: 同一 ChEMBLID が複数ランクに登場することがある（Section 3 の重複行に対応）。
%[text] 各行は独立した IC50 測定値であり、最もよく研究された化合物で特に見られる。

    logInfo("タニモトスコア統計:");
    logInfo("  最大値  : %.3f", max(rankedScores));
    logInfo("  中央値  : %.3f", median(rankedScores));
    logInfo("  最小値  : %.3f", min(rankedScores));

    % タニモト > 0.4（スキャフォールド類似閾値）の化合物数
    nSimilar = sum(rankedScores > 0.4);
    logInfo("エルロチニブとタニモト > 0.4: %d / %d 件", nSimilar, nRanked);

    % 上位ヒットの構造を可視化する
    logInfo("上位ヒットを可視化: %s (IC50 = %.1f nM, タニモト = %.3f)", ...
        metaTbl.MoleculeChEMBLID(rankedIdx(1)), ...
        metaTbl.Value_nM(rankedIdx(1)), ...
        rankedScores(1));
    emk.viz.draw2d(validMols{rankedIdx(1)}, ...
        Title=sprintf("%s IC50=%.1f nM Tan=%.3f", ...
            metaTbl.MoleculeChEMBLID(rankedIdx(1)), ...
            metaTbl.Value_nM(rankedIdx(1)), ...
            rankedScores(1)));
end
%[text] **✏️ やってみよう 7**
%[text] Q: エルロチニブに最も類似する化合物のタニモトスコアは？
%[text]    構造的に非常に近い（タニモト \> 0.7）か、中程度の類似か？
%[text] Q: 上位化合物を ChEMBL または PubChem で調べてみよう。
%[text]    既知の EGFR 阻害剤か？名前にエルロチニブとの関係が示唆されるか？
%[text]    ヒント:
%[text]      topChEMBLID = rankTbl.ChEMBLID(1);
%[text]      disp(rankTbl(1, ["ChEMBLID","Name","Tanimoto"]))
%[text] Q: クエリ化合物をエルロチニブとゲフィチニブで比較する。
%[text]    ゲフィチニブ SMILES: "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"
%[text]    同じ活性セットをゲフィチニブに対して順位付けし、上位 10 件を比較してみましょう:
%[text]    ヒント:
%[text]      mol\_gef = emk.mol.fromSmiles("COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1");
%[text]      fp\_gef  = `emk.fingerprint.morgan(mol_gef, Radius=2, NBits=2048)`;
%[text]      res\_gef = `emk.similarity.rankBy(fp_gef, dbFps, 10, Metric="tanimoto")`;
%[text]    エルロチニブとゲフィチニブは同じ上位化合物を特定するか？
%[text]    EGFR 阻害剤の化学多様性について何が言えるか？
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 8: まとめと重要ポイント
%[text] ### このストーリーで行ったこと
%[text] 1. ターゲットの特定: `emk.db.searchChemblTarget` を使って、ChEMBLターゲットデータベースからEGFR（CHEMBL203）を特定しました。
%[text] 2. データのダウンロード: `emk.db.getChemblActivity` を使い、ChEMBLの査読済み文献からEGFRの$ &dollar&;\\textrm{IC}\_{50}&dollar&;$ レコードを50件取得しました。
%[text] 3. 効果量の解析: $&dollar&;\\textrm{IC}\_{50}&dollar&;$ 分布を調べ、対数スケールで直感的に理解するために $&dollar&;\\textrm{pIC}\_{50}&dollar&;$ に変換し、創薬における効果量の閾値を学びました。
%[text] 4. 活性フィルタリング: $&dollar&;{\\textrm{IC}}\_{50} \\le 100\\;\\textrm{nM}&dollar&;$ の活性カットオフを適用し、高品質なリード化合物のみを選別しました。
%[text] 5. データの保存: `emk.io.writeSdf` を使用し、化学構造交換の業界標準フォーマットであるSDFファイルに活性化合物セットを保存しました。
%[text] 6. 類似度スクリーニング: ECFP4フィンガープリントを計算し、既存薬エルロチニブに対するタニモト類似度で活性化合物を順位付けするリガンドベースのヒット選択戦略を体験しました。 \
%[text] ### 導入した重要概念
%[text] - ChEMBL: 薬物様小分子のキュレーション済みバイオアクティビティデータベース
%[text] - $&dollar&;\\textrm{IC}\_{50}&dollar&; / &dollar&;\\textrm{pIC}\_{50}&dollar&;$: 効果量の指標と対数スケール思考
%[text] - 活性カットオフ: 創薬ステージに応じた適切な閾値の選択
%[text] - SDFフォーマット: 化学構造とメタデータを一括管理する業界標準フォーマット
%[text] - 構造類似度: フィンガープリントを用いたリガンドベースのヒット選択 \
%[text] ### 次のステップ
%[text] - S04 バーチャルスクリーニング: FDA承認薬ライブラリに対する完全なリガンドベース・バーチャルスクリーニング（LBVS）の実施
%[text] - A01 化学空間マッピング: 主成分分析（PCA）を用いたEGFR阻害剤の多様性の可視化
%[text] - A03 QSAR 回帰: 化合物記述子から $&dollar&;\\textrm{pIC}\_{50}&dollar&;$ を予測する機械学習モデルの構築 \
logInfo("S07 完了 -- EGFR バイオアクティビティ解析終了。");
if exist("runDir", "var")
    logInfo("ファイル保存先: %s", runDir);
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
