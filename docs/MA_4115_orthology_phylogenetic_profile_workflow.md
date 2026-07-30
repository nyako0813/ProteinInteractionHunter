# MA_4115 orthology / phylogenetic profile real-data workflow

作成日: 2026-07-30

開始 commit: 6dcb95857f4c65990bdc50c45104be27d5ebdf11
対象: Methanosarcina acetivorans C2A, GCF_000007345.1, WP_011024006.1 (MA_RS21465 / MA_4115, 200 aa)

## 目的と境界

本 pilot は real-data panel、query-only RBH audit、small/full OrthoFinder、formal converter、三状態 profile、coverage-only 非伝播監査を再現可能にする。synthetic fixture は実データへ流用していない。RBH は operational proxy であり formal orthology の真値ではない。結果確認後の panel 変更や threshold 調整は行っていない。coverage-only run では orthology、profile、scoring、rank、tier を無効のまま維持した。

query FASTA は /home/nyako/projects/ProteinHunter_v5/data/databases/target/methanosarcina_acetivorans/ncbi_dataset/data/GCF_000007345.1/protein.faa。regular/readable、SHA256 2eb0af39dece3e0c30bc673d4e10fab5888de221f8b4db6c71b658053062e6a2、4,627 records、query は一意であることを実行前に確認した。

## 既存 engine 契約

Orthology formal columns は protein_id、reference_id、ortholog_id、reference_organism、identity、query_coverage、subject_coverage、evalue、orthogroup、relationship、paralog_ambiguity、source、source_record_id。source=local_table が正式接続で computed は未実装。one_to_one、one_to_many、many_to_one、many_to_many を保持でき、duplicate identity、unknown ID、conflicting duplicate、malformed row は拒否する。OrthoFinder にない identity/coverage/evalue は推測せず空欄。rule mvp1g-orthology-v1、component orthology、category evolutionary。schema 変更は不要で converter の追加だけで接続できた。

Profile required columns は protein_id、species_id、presence、optional は taxonomic_group、source、source_record_id。presence は true、false、空欄の三状態。duplicate protein/species は拒否、順序 deterministic。rule mvp1h-phylogenetic-profile-v1 は jointly informative = shared presence + shared absence + discordant、similarity = (shared presence + shared absence) / jointly informative。既定閾値は shared 2、informative 3、similarity 0.8。missing は denominator に入らない。

両 component は evolutionary category と cap 2.0 を共有するが、同じ OrthoFinder source 由来で非独立。今回は scoring 自体を実行していない。disabled local_table も validation と manifest fingerprint の対象にする後方互換な最小変更を加え、engine には伝播させない。

## Tool readiness

取得日 2026-07-30。official docs、official repository、installed help を source of truth とした。

- NCBI Datasets 18.33.1 と同 release の dataformat。isolated path /home/nyako/tools/ncbi-datasets/2026-07-30。datasets SHA256 8459ef1e87433f7b1198f5703c8cc10b55f1904cd448cf7e996c0892b141cd1f、dataformat SHA256 8450cf7cbdb0ed7fece567405732cd1ff838b5352faa58abbccdeff56e1ff0e8。公式: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/ と https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/
- OrthoFinder 3.1.5、GPL-3.0、isolated path /home/nyako/tools/orthofinder/3.1.5/.venv/bin/orthofinder。archive SHA256 d74b5dbf9348e7ffdb735f1fda4ff9c49f504261811355519cf9dafcf80e7739。旧 repository でなく現行 OrthoFinder/OrthoFinder、https://orthofinder.github.io/OrthoFinder/download_and_install/、https://orthofinder.github.io/OrthoFinder/tutorials/guide-to-results/ を確認。
- bundled DIAMOND 2.0.13、GPL-3.0。query audit と default search backend に使用。公式: https://github.com/bbuchfink/diamond。system BLASTP 2.12.0+ も確認したが不採用。
- bundled binary checksums: DIAMOND d4049613dd6b7ca114f3f55eb62027aa4625fb184f5fbf510b8a9b477624240e、FAMSA 2.2.3 20b52ee56519facaae159b61a40f424835ff6f0e6b7492bf632628b8c06f7ed3、FastTree 2.1.11 6fe112c32a2e0b5a2cc6baeeb5f564300f5f650b0dfe96b54a06f0c2ef7a97a5、MCL 14-137 3a8790a79e3709393255c1c2337318b16db07efdfaf50f5b24cb4a2f9d2a6b03。default の DIAMOND search、FAMSA MSA、FastTree tree、MCL clustering を固定。system Python は変更していない。

OrthoFinder は species tree、orthogroup、pairwise relation、HOG、core/assign と output reuse を提供する。MMseqs2 は今回だけのためには導入しなかった。RBH は query audit には有用だが gene-tree based orthology の代替にしない。

## Panel policy と members

policy は orthology_panel_policy.v1.yaml、search policy は orthology_search_policy.v1.yaml として結果を見る前に固定。versioned RefSeq、complete/chromosome、annotation、protein count 500..12,000、one representative per species、contig/MAG/atypical/obsolete 除外、query 固定。CheckM missing だけでは reference isolate を除外しない。bacterial outgroup は含めない。

30 targets 中27 assemblies。layers は Methanosarcina 7、other Methanosarcinales 5、other methanogens 7、non-methanogenic Archaea 8。fixed manifest SHA256 6edef15a3175d96ec3ac27b2293c03f75ddc79c359c28321a9a7cd078846c074、panel SHA256 5df96d00ccf8bc14ed11b4d2f6fda7d50b649ebe933d51740f9ab2d469ebeeab、panel audit 9cfdc978082c94a078ce704c4b9047948621854a3d6e451bea0b62352c60db4a、panel metadata 0efea0129ef863f37eaf98306ca4c66bb54d0d08d41162a2d01b13cc3587acb0。

| Layer | Versioned accessions |
|---|---|
| Methanosarcina | GCF_000007345.1, GCF_000970025.1, GCF_000970285.1, GCF_000970205.1, GCF_000970085.1, GCF_000969885.1, GCF_000969905.1 |
| Methanosarcinales | GCF_000013725.1, GCF_000025865.1, GCF_000328665.1, GCF_000217995.1, GCF_000204415.1 |
| Other methanogens | GCF_001458655.1, GCF_000016525.1, GCF_000011005.1, GCF_000091665.1, GCF_002945325.1, GCF_000015825.1, GCF_000013445.1 |
| Non-methanogenic Archaea | GCF_000008665.1, GCF_004799605.1, GCF_000025685.1, GCF_000018465.1, GCF_000012285.1, GCF_008245085.1, GCF_000009965.1, GCF_000195915.1 |

policy eligible assembly が返らず除外した targets は Methanolobus psychrophilus、Methermicoccus shengliensis、Methanomassiliicoccus luminyensis。selected の Methanococcus maripaludis と Methanoculleus marisnigri JR1 は taxonomy conflict で manual review。metadata missing は0。手動推測で修正していない。

## Download、QA、normalization

    datasets download genome accession --inputfile orthology_panel_accessions.txt --include protein --filename orthology_panel.zip --no-progressbar

ZIP は 14,079,514 bytes、SHA256 5c693b8953985b3c13c2302f70b85bfc89bfdc489dc67de9543257febcacafe4、3.92 s、peak RSS 17,072 KiB、CRC OK。catalog/md5、27 directories、27 protein.faa、expected list 一致を確認。downloaded query FASTA は期待 SHA と一致。

prepare_orthology_proteomes.py は species ID GCF_000007345_1、protein ID GCF_000007345_1__WP_011024006.1 を生成し、raw ID、header、assembly、checksums を mapping に保持。deterministic sort、UTF-8/LF。27 species、73,631 proteins、duplicate IDs、empty、zero、invalid、internal stop は0。non-standard、extreme length、exact duplicate sequence は audit に保持し silent truncation なし。query round-trip 1/1。mapping SHA256 a3784d2ecb5f6b1011c8acc56784dc06c582ff2ffbb6c58faf43a03ae6fa804d、audit e5824ee00375896c0a6cb949a24c941c0382e74ecb42b4ad6c2f6739771054b2。

## Query-only RBH

固定 policy は more-sensitive、E-value 1e-5、top25、max HSP1、bitscore40、identity20 percent、alignment60 aa、query/subject coverage0.5、length ratio0.5..2.0、tie1 percent、secondary0.9、unique reciprocal best hit。species ごとに search を分離し DIAMOND default filtering を維持。

26 species 中 unique_RBH proxy 22、no_detectable_homolog 4、ambiguous/multi-copy/fragment 0。no-hit は GCF_000008665.1、GCF_000018465.1、GCF_000012285.1、GCF_000195915.1。16.20 s、peak RSS123,128 KiB。audit SHA256 1ad9c16cf4316be56f12cfdf7558bffbe187755c23e593bf019de0ad0cd6a8d3、summary 52e6ebcd27345bc9c40bc8a5cc09cc281beccd7df8895c5b11384be043cba4d3。formal source にはしていない。

## Small/full OrthoFinder

small 14 species command:

    orthofinder -f small_input -t 4 -a 2 -n MA4115_small_v1 -o small_results

exit0、6:25.56、peak RSS728,336 KiB、disk263 MiB。41,143 genes、34,539 assigned、6,604 unassigned、4,472 groups、216 single-copy。MA_4115 は OG0001023、query copy1、comparison orthologues11、group copies は全て1、paralog0。ID round-trip と schema checks は合格。

full command:

    orthofinder -f normalized -t 4 -a 2 -n MA4115_full_v1 -o full_results

start 2026-07-30T12:39:38Z、end12:55:51Z、exit0、16:13.13、peak RSS621,832 KiB、disk540 MiB。27 species、73,631 genes、65,235 assigned、8,396 unassigned、5,782 groups、153 all-species single-copy、2,706 groups は一 species 以上で multi-copy。species tree SHA256 4fed39791d9b562bf99c9140e0f7d8f27d925291e47846bfcb1b1307c1e06768、query pairs file d712629141adeb4895ee3e9c360c4336ce530c83d75fa710ca6e6844c5489d0f、Orthogroups.tsv 7b3e8300ec6658799890361648ac4bb1e0a474246c2b6369d38fb17e2464bbd8。warning は multiple potential species-tree roots、selected outgroup は GCF_000018465_1 と GCF_000012285_1。

MA_4115 は full OG0000982、group total23、23 species present、全 copy1。comparison ortholog species22、one_to_one22、paralog0、ambiguity0、fragment-only0。RBH と species 22/22、exact IDs 22/22 一致、conflict0。残る4 valid species の not-detected は biological absence の確定ではない。

## Converter と orthology coverage

convert_orthofinder_orthology.py は official consolidated pair file、Orthogroups.tsv、species mapping、raw ID round-trip を検証し、全 relation と ambiguity を保持。formal SHA256 c47f0cea2ac43d643e42700b9f1dda240ca6b0f7214bc84a74c7d5b63ac6d4bc、audit 668852a1d3f5773d387bdd9abd9c69d1c4fa24d303d785bdbe23c0d89596cb89、coverage 826625c5b119f328a1401a7b62fab70be641cad8d82b181e6d73b6f6c41645d9、metadata baf486230ac75bb4fdb2b87668d376a974c96e131a1743329205b9b9d0c5c116。54,467 records、duplicates0、malformed0、unknown IDs0。

| Metric | Count |
|---|---:|
| M. acetivorans proteins | 4,627 |
| assigned / unassigned | 4,231 / 396 |
| external ortholog / none | 4,056 / 571 |
| one-to-one | 3,456 |
| one-to-many | 1,085 |
| many-to-one | 2,261 |
| many-to-many | 1,144 |
| paralog ambiguity | 2,951 |
| fragment-only | 0 |

relation metrics は相互排他的ではなく、assigned と external ortholog を分離した。

## Profile construction と coverage

build_phylogenetic_profiles.py は4,627 x 26 = 120,302 cells。mapping は present_unique=true、present_multi_copy=true、present_ambiguous=missing、fragment_only=missing、not_detected=false only for valid/evaluable proteome、species_missing/proteome_invalid/not_evaluated=missing。missing を false にしない。

formal columns は protein_id、species_id、presence、taxonomic_group、source、source_record_id。formal SHA256 7447e4668f1ad4e1efa85a43ef1947edcdff23518df4255644290d021d6a6295、audit 03fb56621a65a7c93e262e94b4a4fe50224c3fc324ca1cdeff735ffe75a8d467、pair audit 2d7ccddd6328a7dd578e2634290f9a813b01780b91fc65f2be059f789ef4ce4f、coverage dfd02667737706465e290fb4346ae33dd614e411d74ef376cb64a6a02bb39ef6、metadata cf803b69ed69416290a1b8c67b34e19581da940adda5b89a5975dd61ed50b162。

states: unique31,346、multi-copy1,866、ambiguous14,414、not-detected72,676。missing cells0、uncertain14,414。informative minimum を満たす profiles4,623、不足4。query informative26、present22、not-detected4、uncertain0、comparable candidates4,623、not comparable4。

similarity values4,625、min0、Q1 0.2、median0.347826、Q3 0.576923、max1.0。全 threshold pass586、self 除外585、0.79..0.81 は118。similarity 1.0 の例は WP_011020101.1、WP_011020164.1、WP_011020274.1、WP_011020594.1、WP_011020654.1。interaction truth とは解釈しない。pass pair で shared absence が numerator の過半を占めるものは0だが metric の inflation risk は残る。

## Coverage-only A/B

pilot_orthology_coverage.yaml は local tables を validation/manifest のみで読み、orthology/profile/domains/localization/scoring/tiers disabled。validation は orthology54,467 records / 4,056 proteins / unknown0、profile120,302 observations / 4,627 proteins / 26 species / unknown0。candidate pairs4,627、non-self4,626、included3,795、flagged831、excluded1。

localization coverage baseline と run_id 除外で candidate TSV 4,627全行、JSONL 4,627全行が一致し、warning summary は bytes 一致。candidate order、disposition、fragment、coordinates、context、operon、annotation、domain/localization は不変。Excel は同じ22 sheets/dimensions、orthology/profile sheets は header-only。manifest 差は run/config と追加 optional input fingerprints。全 bundle で orthology/profile arrays 空、statuses not_run、score/rank/tier 空。

## Readiness、限界、次段階

Orthology coverage-only readiness は ready。source/version、panel、ID round-trip、relations、ambiguity、loader、provenance を監査可能。ただし tree/root uncertainty、assembly/annotation error、remote homology false negative、paralog reconciliation error、pair output に identity/coverage/evalue がない制限が残る。

Profile coverage-only readiness は ready。三状態は保持。ただし species は等重みで taxonomic non-independence を補正せず、orthology と同じ source。shared-absence metric、valid proteome の non-detection を false とする operational assumption、14,414 ambiguous cells が主要リスク。

Formal scoring readiness は not ready / not authorized。root sensitivity、manual taxonomy flags、ambiguous groups、high-similarity candidates の independent evidence、correlation policy を別途審査する。threshold、weight、cap、tier 条件は変更していない。

tracked 対象は scripts、tests、policies、fixed accessions、config、本書のみ。ZIP、proteomes、normalized FASTA、DIAMOND DB/search、OrthoFinder raw results/trees、formal real-data TSV、Excel、cache/logs は非追跡。quality gates は pytest 393 passed、Ruff all checks passed、mypy 127 source files success、git diff --check success。large-file audit と push verification も commit 前後に実施する。
