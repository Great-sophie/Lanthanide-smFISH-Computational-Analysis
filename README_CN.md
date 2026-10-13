# Chapter 3｜HER2 smFISH 计算成像扩展

对应论文：*Lanthanide Complex for Single-Molecule Fluorescent In Situ Hybridization and Background-Free Imaging*，*Analytical Chemistry* (2024)。
DOI：https://doi.org/10.1021/acs.analchem.3c04530

该仓库整合两类**独立**分析：

1. **MDA-MB-231、BT474、SKBR3 的 Eu 红通道共享尺度局部对比度图。** 从 RGB 8-bit 导出图提取红通道，采用两尺度高斯差分计算局部对比度，并在三张图之间使用统一的显示尺度。所得结果是显示信号的空间特征，不等于 HER2 mRNA 的绝对表达量，也不验证 231 必然最弱。
2. **SKBR3 细胞在 25%、50%、75%、100% dUTP-N3 条件下的 TrackMate 单细胞斑点计数。** 实验者确认每份 `<梯度>-<编号>.csv` 对应一个细胞，统计每细胞候选斑点数、中位数、均值、样本方差及 QUALITY 阈值敏感性。原始梯度标签不交换。文件名与文件夹梯度冲突的 CSV 会被排除并审计。

原始 TIFF/CSV 不包含在本仓库；运行需要本地真实实验数据。单细胞不是独立生物学重复。只有获得真正同一细胞的 ImageJ 计数，才进行两种算法的配对比较。

## Ubuntu 运行

```bash
conda activate pathfm
python -m pip install -r requirements.txt
python -m pytest -q tests/

GRADIENT_DATA="/media/sophie/EXTERNAL_USB/时间分辨FISH/Time-resolved FISH PPT/Figure3 SKBR3 dutp-N3梯度/fig3 0710"
MAPPING_DATA="$GRADIENT_DATA/Her2 mRNA mapping"

bash scripts/run_all.sh "$GRADIENT_DATA" "$MAPPING_DATA"
```

检查 `figures/` 中的三张新图和 `outputs/spot_counts/file_audit.csv`，再决定是否上传图像。`outputs/` 被 Git 忽略。

主文档见 `README.md`。
