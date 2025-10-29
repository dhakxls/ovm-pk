---
name: OVM-PK Architecture
alwaysApply: true
---
# Project Architecture
- Package: `src/ovmpk/` with docking/, fetchers/, prep/, utils/
- Scripts: env/bootstrap + MD runners in `scripts/`
- Configs: `configs/` (default/fast_test/prod_test)
- Data: inputs in `data/input`, working in `data/work`, results in `data/output`
- Tests: `tests/` (fetchers → prep → docking → pose_selection → md_* → others)
