# nf-core/sarek（中文使用说明）

> 面向第一次接触 **Linux / Nextflow / nf-core/sarek** 的用户。
> 这份 README 基于当前仓库内容整理，目标是“能按步骤跑起来”。

## 1. 项目简介

`nf-core/sarek` 是一个用于 DNA 测序变异检测的 Nextflow 流程，支持：

- **Germline（胚系）变异**分析
- **Somatic（体细胞）变异**分析
- **肿瘤-正常配对**与肿瘤单样本场景
- WGS / WES / 靶向数据（取决于输入和参数）

它与 nf-core 的关系：

- 本项目就是 nf-core 官方维护的 `sarek` 流程仓库。
- 使用 Nextflow DSL2 编写，按模块调用工具，便于复现和迁移。

当前仓库中可见的常用step/toosl（按类别举例）：

- 比对mapping/markduplicates/prepare_recalibration/recalibrate：`bwa-mem`、`bwa-mem2`、`dragmap`、`GATK` 相关步骤
- 小变异调用variant_calling：`Mutect2`、`HaplotypeCaller`、`Strelka`、`freebayes`、`lofreq` 等
- 结构变异variant_calling / 拷贝数：`manta`、`tiddit`、`ascat`、`cnvkit`、`controlfreec`
- 注释annotate：`snpEff`、`VEP`、`SnpSift`、`bcftools annotate`
- 质控汇总：`MultiQC`

通俗理解：

- **Manta / TIDDIT**：偏向检测结构变异（大片段变化）。
- **ASCAT / CNVkit / Control-FREEC**：偏向检测拷贝数变化（CNV）。
- **Mutect2 / Strelka 等**：偏向检测 SNP/Indel（小变异）。

---

## 此流程涉及的软件较多，并且涉及大量的数据库，流程支持任何一步开始。为目前测试的最复杂的流程。进行变异检测软件的选择

<table>
  <thead>
    <tr>
      <th>Tool</th><th>WGS</th><th>WES</th><th>Panel</th><th>Germline</th><th>Tumor-Only</th><th>Somatic (Tumor-Normal)</th>
    </tr>
  </thead>
  <tbody>
    <tr><td>DeepVariant</td><td>✓</td><td>✓</td><td>✓</td><td>✓</td><td>-</td><td>-</td></tr>
    <tr><td>FreeBayes</td><td>✓</td><td>✓</td><td>✓</td><td>✓</td><td>✓</td><td>✓</td></tr>
    <tr><td>GATK HaplotypeCaller</td><td>✓</td><td>✓</td><td>✓</td><td>✓</td><td>-</td><td>-</td></tr>
    <tr><td>GATK Mutect2</td><td>✓</td><td>✓</td><td>✓</td><td>-</td><td>✓</td><td>✓</td></tr>
    <tr><td>lofreq</td><td>✓</td><td>✓</td><td>✓</td><td>-</td><td>✓</td><td>-</td></tr>
    <tr><td>mpileup</td><td>✓</td><td>✓</td><td>✓</td><td>✓</td><td>✓</td><td>-</td></tr>
    <tr><td>Strelka</td><td>✓</td><td>✓</td><td>✓</td><td>✓</td><td>-</td><td>✓</td></tr>
    <tr><td>Manta</td><td>✓</td><td>✓</td><td>✓</td><td>✓</td><td>✓</td><td>✓</td></tr>
    <tr><td>indexcov</td><td>✓</td><td>-</td><td>-</td><td>✓</td><td>-</td><td>✓</td></tr>
    <tr><td>TIDDIT</td><td>✓</td><td>✓</td><td>✓</td><td>✓</td><td>✓</td><td>✓</td></tr>
    <tr><td>ASCAT</td><td>✓</td><td>✓</td><td>-</td><td>-</td><td>-</td><td>✓</td></tr>
    <tr><td>CNVKit</td><td>✓</td><td>✓</td><td>-</td><td>✓</td><td>✓</td><td>✓</td></tr>
    <tr><td>Control-FREEC</td><td>✓</td><td>✓</td><td>✓</td><td>-</td><td>✓</td><td>✓</td></tr>
    <tr><td>MSIsensorPro</td><td>✓</td><td>✓</td><td>✓</td><td>-</td><td>✓</td><td>✓</td></tr>
  </tbody>
</table>


## 2. 输入数据说明

### 2.1 你至少需要准备什么

最基础要有：

1. 一个样本表（CSV，传给 `--input`）
2. 原始数据或中间数据（FastQ / BAM / CRAM 等，取决于 `--step`）
3. 参考基因组与已知位点资源（可用 `--genome` 或本地 `--fasta` + 其他文件）

### 2.2 BAM / CRAM / 索引是什么

- **BAM**：二进制比对结果文件。
- **CRAM**：更省空间的比对结果格式。
- **BAI / CRAI**：对应 BAM / CRAM 的索引文件（让工具可快速随机访问）。

如果从 `markduplicates` 或更后步骤启动，通常需要 BAM/CRAM 及其索引。

### 2.3 肿瘤-正常配对是什么意思

在样本表中：

- 同一个 `patient`
- 不同的 `sample`
- 用 `status` 区分：`0` 正常、`1` 肿瘤

这样流程能把同一患者的正常样本与肿瘤样本配对做体细胞分析。

可以做什么：
小变异：mutect2 或 strelka（也有人两者都跑）
结构变异：manta,tiddit
CNV：ascat / cnvkit / controlfreec
MSI：msisensor2 或 msisensorpro
注释：vep 或 snpeff

### 2.4 样本表（CSV）做什么

`--input` 指向的 CSV 告诉流程：

- 有哪些样本
- 每个样本对应哪些文件
- 是哪种输入类型（FASTQ / BAM / CRAM / VCF 等）
- 是否配对、是否多 lane

仓库自带示例：`assets/samplesheet.csv`。

### 2.5 `params.yaml` / `custom.config` / `nextflow.config` 的区别

- **`params.yaml`**：放“流程参数”（如 `input`、`outdir`、`tools`、`fasta`）。推荐新手使用。
- **`custom.config`**：放“Nextflow 运行层配置”（如资源、executor、容器参数）。
- **`nextflow.config`**：仓库默认总配置，定义了很多默认参数与 profile。

> 注意：nf-core 建议不要用 `-c` 去传业务参数，业务参数优先放命令行或 `-params-file`。

---

## 3. 运行前准备（新手清单）

请按下面清单逐项确认：

- [ ] 已安装可用版本的 Nextflow（本仓库 README 徽章显示当前要求为 `>=25.10.2`）
- [ ] 已准备容器运行方式（常见为 `singularity` / `apptainer` 或 `docker`）
- [ ] 已准备本地参考数据（FASTA、dbSNP、known_indels、intervals 等）
- [ ] 已规划容器缓存目录（避免反复拉镜像）
- [ ] 已规划工作目录 `work/` 和结果目录 `result/`
- [ ] 磁盘空间充足（WGS 通常非常占空间）
- [ ] 所有路径尽量使用**绝对路径**

---

## 4. 推荐目录结构（新手版）

```text
/project/sarek_run/
├── input/
│   ├── samplesheet.csv
│   └── bam_or_fastq/
├── reference/
│   ├── genome.fa
│   ├── genome.fa.fai
│   ├── genome.dict
│   ├── dbsnp.vcf.gz
│   ├── dbsnp.vcf.gz.tbi
│   ├── known_indels.vcf.gz
│   └── intervals.bed
├── config/
│   ├── params.yaml
│   └── custom.config
├── work/
├── result/
├── singularity_cache/
└── logs/
```

- `input/`：输入样本表和原始/中间数据
- `reference/`：本地参考与注释资源
- `config/`：参数文件与运行配置
- `work/`：Nextflow 任务中间目录（排错最关键）
- `result/`：最终发布结果目录（`--outdir`）
- `singularity_cache/`：容器缓存（强烈建议长期保留）
- `logs/`：你手工保存的运行日志

---

## 5. 最基础运行方法（一步一步）

下面给一个“能直接照抄改路径”的最小流程。

### 5.1 准备输入样本表

最简 FASTQ 示例（与仓库示例格式一致）：

```csv
patient,sample,lane,fastq_1,fastq_2
PATIENT1,SAMPLE1,L001,/abs/path/S1_R1.fastq.gz,/abs/path/S1_R2.fastq.gz
```

若是肿瘤-正常配对，建议补 `status`（`0` 正常，`1` 肿瘤）和 `patient` 对应关系。

### 5.2 准备 `params.yaml`

```yaml
input: '/project/sarek_run/input/samplesheet.csv'
outdir: '/project/sarek_run/result'
step: 'mapping'
(根据当前仓库的参数定义，step 可选：mapping（默认）, markduplicates,prepare_recalibration,recalibrate,variant_calling,annotate)

# 二选一：
# 1) 用 iGenomes 名称（默认是 GATK.GRCh38）
genome: 'GATK.GRCh38'
# 2) 或显式提供本地参考（更适合离线/内网）
# fasta: '/project/sarek_run/reference/genome.fa'
# dbsnp: '/project/sarek_run/reference/dbsnp.vcf.gz'
# known_indels: '/project/sarek_run/reference/known_indels.vcf.gz'
# intervals: '/project/sarek_run/reference/intervals.bed'

# 按需选择工具，可多个逗号分隔
tools: 'manta,tiddit,ascat'
(tools可选ascat,bbsplit,bcfann,cnvkit,controlfreec,deepvariant,freebayes,haplotypecaller,indexcov,lofreq,manta,merge,mpileup,msisensor2,msisensorpro,muse,mutect2,ngscheckmate,sentieon_dedup,sentieon_dnascope,sentieon_haplotyper,sentieon_tnscope,snpeff,snpsift,strelka,tiddit,vep,varlociraptor)

# 常见离线场景建议
igenomes_ignore: true
save_output_as_bam: false
```

### 5.3 准备 `custom.config`（可选）

仓库里没有现成 `custom.config`，新手可自己建一个最小版，例如：

```groovy
process {
  cpus = 8
  memory = '32 GB'
  time = '24h'
}

singularity {
  enabled = true
  autoMounts = true
}
```

### 5.4 设置环境变量（以 Singularity/Apptainer 为例）

```bash
export NXF_HOME=/project/sarek_run/.nextflow
export NXF_WORK=/project/sarek_run/work

export NXF_DISABLE_CHECK_LATEST=true
export NXF_SINGULARITY_CACHEDIR=/data/person/wup/liusy/sarek/singularity_cache
export SINGULARITY_CACHEDIR=/data/person/wup/liusy/sarek/singularity_cache
export APPTAINER_CACHEDIR=/data/person/wup/liusy/sarek/singularity_cache
```

### 5.5 执行运行

在仓库目录（`/data/person/wup/liusy/sarek`）下：

```bash
conda activate nextflow
cd /data/person/wup/liusy/sarek
```
```bash
nextflow run . \
  -profile singularity \
  -params-file /data/person/wup/liusy/sarek/test/params.yaml \
  -c /data/person/wup/liusy/sarek/test/custom.config
```
**WGS的manta、TIDDIT、ASCAT**
input: /data/person/wup/liusy/sarek/wgs_manta_bqsr.csv
```bash
nohup nextflow run . \
  -ansi-log false \
  -profile singularity \
  -params-file /data/person/wup/liusy/sarek/test/params_wgs_manta.yaml \
  -c /data/person/wup/liusy/sarek/test/custom_wgs.config \
  > sarek.log 2>&1 &
```
**WES数据的mutect2**
```bash
nohup nextflow run . \
  -ansi-log false \
  -profile singularity \
  -params-file /data/person/wup/liusy/sarek/test/params_wes_mutect2.yaml \
  -c /data/person/wup/liusy/sarek/test/custom_wes_mutect2.conf \
  > /data/person/wup/liusy/sarek/sarek.log 2>&1 &

```
**WES数据variant_calling: ascat, manta, tiddit**
```bash

nohup nextflow run . \
  -ansi-log false \
  -profile singularity \
  -params-file /data/person/wup/liusy/sarek/test/params_wes_manta.yaml \
  -c /data/person/wup/liusy/sarek/test/custom_wes_mutect2.conf \
  > /data/person/wup/liusy/sarek/sarek.log 2>&1 &

```

如果你跑的是官方远程仓库版本，可用(不推荐，网不好)：

```bash
nextflow run nf-core/sarek -r <版本号> -profile singularity -params-file params.yaml
```

### 5.6 中断后继续跑

```bash
nextflow run . \
  -profile singularity \
  -params-file /data/person/wup/liusy/sarek/test/params.yaml \
  -c /data/person/wup/liusy/sarek/test/custom.config \
  -resume

#nohup后台进行
nohup nextflow run . \
  -ansi-log false \
  -profile singularity \
  -params-file /data/person/wup/liusy/sarek/test/params.yaml \
  -c /data/person/wup/liusy/sarek/test/custom.config \
  -resume \
  > sarek.log 2>&1 &
```

`-resume` 会复用已成功任务，避免重算。

---

## mutect2得到的wgs/wes的vcf.gz文件进行annovar注释
提供sarek_mutect2_annovar_manifest.tsv
| pair_id | tumor_sample | normal_sample | vcf |
| :--- | :--- | :--- | :--- |
| HP_tumor_vs_HP_normal | HP_tumor | HP_normal | /data/person/wup/liusy/sarek/test/result/variant_calling/mutect2/HP_tumor_vs_HP_normal/HP_tumor_vs_HP_normal.mutect2.filtered.vcf.gz |
| XHS_tumor_vs_XHS_normal | XHS_tumor | XHS_normal | /data/person/wup/liusy/sarek/test/result/variant_calling/mutect2/XHS_tumor_vs_XHS_normal/XHS_tumor_vs_XHS_normal.mutect2.filtered.vcf.gz |

```bash
cd /data/person/wup/liusy/wgs/scripts
```

```bash
sbatch --array=1-${N}%8   run_annovar_sarek_mutect2_array.slurm   sarek_mutect2_annovar_manifest.tsv   /data/person/wup/liusy/sarek/annotation/annovar/   /data/person/wup/liusy/sarek/tmp_annovar_mutect2
```



## 6. 重要参数解释（新手版）

以下是最常用且最容易混淆的参数：

- `input`
  - 做什么：指定样本表 CSV 路径。
  - 何时改：每次换批次样本都要改。

- `outdir`
  - 做什么：指定最终结果输出目录。
  - 何时改：每次新项目建议新目录。

- `step`
  - 做什么：指定从哪一步开始，常见有 `mapping`、`markduplicates`、`prepare_recalibration`、`recalibrate`、`variant_calling`、`annotate`。
  - 何时改：你有中间产物，想从流程中段接着跑时。

- `tools`
  - 做什么：指定使用哪些工具（逗号分隔），涵盖变异检测/注释等。
  - 何时改：按分析目标选择（如结构变异、CNV、小变异等）。

- `fasta`
  - 做什么：本地参考基因组 FASTA 路径。
  - 何时改：不用 `--genome` 或离线运行时基本都要改。

- `dbsnp`
  - 做什么：dbSNP 已知位点 VCF。
  - 何时改：做 BQSR、注释或相关调用时通常需要。

- `known_indels`
  - 做什么：已知 indel 位点文件（常用于重校准相关步骤）。
  - 何时改：走 GATK 重校准链路时应提供。

- `intervals`
  - 做什么：目标区域（BED/interval）文件。
  - 何时改：WES/Panel 场景通常要改。

- `genome`
  - 做什么：使用 iGenomes 里的参考名（默认 `GATK.GRCh38`）。
  - 何时改：你想直接用 iGenomes 预置参考时。

- `igenomes_ignore`
  - 做什么：是否忽略 iGenomes 配置。
  - 何时改：你完全使用本地参考时常设为 `true`。

- `save_output_as_bam`
  - 做什么：让预处理输出保存为 BAM（否则常见是 CRAM）。
  - 何时改：下游工具只接受 BAM 或团队规范要求 BAM 时。

- `-profile`
  - 做什么：选择运行环境配置（如 `singularity`、`docker`）。
  - 何时改：根据你实际计算平台和容器方案。

- `-resume`
  - 做什么：断点续跑，复用缓存。
  - 何时改：中断重跑时几乎都应加。

---

## 7. 本地/离线运行建议（强烈推荐）

1. **尽量使用本地参考数据**
   - 参考文件大且会反复读取，放本地更稳、更快。

2. **尽量使用本地容器缓存**
   - 避免每次拉取镜像，减少网络相关失败。

3. **不要过度依赖远程 S3 路径**
   - 网络抖动、权限、TLS 问题都会导致任务失败或重试。

4. **路径尽量写绝对路径**
   - 在集群/容器环境中，相对路径最容易引发“文件找不到”。

5. **网络差时先检查**
   - 镜像拉取是否成功
   - 参考文件是否可访问
   - DNS / 代理 / 防火墙设置

6. **减少在线下载失败的方法**
   - 预下载参考与注释文件
   - 预热容器缓存（先拉好镜像）
   - 尽量固定运行参数与 profile，减少缓存失效

---

## 8. 常见报错与排查（新手重点）

### 8.1 `Singularity image pull failed`

通常表示：镜像拉取失败（网络/权限/仓库连接问题）。
先查：

1. `NXF_SINGULARITY_CACHEDIR` 是否可写
2. 节点是否能访问镜像源
3. 是否有代理/TLS 限制

### 8.2 `TLS handshake timeout`

通常表示：网络连接慢或被中间网络设备阻断。
先查：

- 当前节点外网连通性
- DNS 解析与代理配置
- 是否可以改为本地缓存避免在线拉取

### 8.3 `connection reset by peer`

通常表示：远端连接被重置（网络或服务端问题）。
先查：

- 是否偶发（可重试）
- 是否只在下载阶段报错
- 是否有并发过高导致连接不稳定

### 8.4 `process terminated with exit status (1)`

通常表示：任务命令执行失败，但“原因不在这行字里”。
先查：

1. 失败任务目录下的 `.command.err`
2. 同目录 `.command.out`、`.command.sh`
3. 顶层 `.nextflow.log`

### 8.5 `Execution is retried (1)`

通常表示：Nextflow 按策略自动重试一次。
这**不等于最终失败**，要看最后状态是否成功。

### 8.6 `Submitted / Cached / Re-submitted / Error` 区别

- `Submitted process`：任务已提交执行。
- `Cached process`：命中缓存，直接复用历史成功结果。
- `Re-submitted process`：任务重试后再次提交。
- `Error`：任务最终失败。

### 8.7 流程中断后怎么继续

- 原目录下直接加 `-resume` 重跑。
- 不要随意删 `work/`，否则缓存会丢失。

### 8.8 如何去 `work` 目录查失败任务

一般步骤：

1. 从日志里找到失败任务的 `work/xx/xxxx...` 路径
2. 进入该目录
3. 重点看：`.command.sh`、`.command.err`、`.command.out`

### 8.9 如何看 `.nextflow.log`

```bash
tail -n 200 .nextflow.log
```

或者按关键词搜索：

```bash
rg -n "ERROR|WARN|Caused by|work/" .nextflow.log
```

---

## 9. 如何判断流程是否真的完成

可以同时看 4 件事：

1. **终端/日志出现成功结束信息**（如完成总结，无未处理错误）
2. `result/` 下有预期目录（如 `pipeline_info/`、`multiqc/`、对应工具输出目录）
3. `.nextflow.log` 最后没有持续失败堆栈
4. 关键样本关键步骤产物存在且非空

`Cached process` 的意思是“之前跑过且成功，这次直接复用”，不是失败。

如果要核对某个样本某一步是否完成：

- 看 `result/` 对应样本目录是否有预期文件
- 回查日志中该任务是否最终 `COMPLETED` 而非 `ERROR`

---

## 10. 给新手的实用建议

1. 先用 1~2 个小样本试跑整链路。
2. 运行中尽量不要频繁改 `-profile`。
3. `singularity_cache` 尽量长期保留。
4. 中断后优先 `-resume`，不要直接重头跑。
5. 出现一次 retry 不一定失败，先看最终状态。
6. 报错先看失败任务目录的 `.command.err`，再看 `.nextflow.log`。

---

## 补充链接

- 官方使用文档入口：`docs/usage.md`（仓库内提示以 nf-core 网站最新版为准）
- 输出说明：`docs/output.md`
- 引用文献与工具引用：`CITATIONS.md`
- 变更记录：`CHANGELOG.md`

如果你是首次部署，建议先使用测试 profile 验证 Nextflow + 容器环境可用，再上真实数据。
