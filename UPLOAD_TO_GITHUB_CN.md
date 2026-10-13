# GitHub 上传指南

1. 在 GitHub 新建一个空仓库：`Lanthanide-smFISH-Computational-Analysis`（不初始化 README）。
2. 下载并解压本项目 ZIP 到 `~/Chapter3_smFISH_Analysis/`。
3. 在项目根目录执行 `python -m pip install -r requirements.txt` 和 `python -m pytest -q tests/`。
4. 按 `README_CN.md` 设置真实数据路径，运行 `bash scripts/run_all.sh "$GRADIENT_DATA" "$MAPPING_DATA"`，确认真实的三张结果图和审计表。
5. 公开上传前，核查图像分享权限和 `git status`，避免包含原始实验数据。
6. `git init -b main && git add . && git commit -m "Add HER2 smFISH computational imaging analysis"`。
7. `git remote add origin git@github.com:Great-sophie/Lanthanide-smFISH-Computational-Analysis.git`。
8. `git push -u origin main`。

如果仓库已带 README，则**不要在本地再次创建独立历史**。先使用 `git clone` 克隆，再将本项目文件复制进去提交，避免分支冲突。
