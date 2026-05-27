# TyWRF-Core 实施计划

## Summary
- 项目名：`TyWRF-Core`
- 项目目录：`/home/zzy/Projects/TyWRF-Core`
- 目标：实现 CPU-first、CUDA-ready 的 WRF-compatible 台风积分器，复用 WPS/real.exe 前处理和现有后处理。
- v1 只支持当前 PGWRF/KROSA 子集；先保证 6 小时循环核心场小 RMSE，再做性能优化和 CUDA 迁移。

## 主 Agent 与子 Agent 管理策略
- 主 agent 由我承担：负责架构冻结、任务拆分、子 agent 调度、代码审查、冲突处理、最终集成和对你汇报。
- 子 agent 只接收我分发的明确任务，不直接和你沟通；你只看我的阶段汇报。
- 每个子 agent 负责独立写入范围，避免并行改同一文件。
- 子 agent 最终必须提交：改动文件列表、设计说明、测试结果、遗留风险。
- 我集成时统一做：
  - 编译检查；
  - 单元测试；
  - 6 小时 WRF reference 对比；
  - 代码风格和 CUDA-ready 约束检查；
  - 汇总报告。

## 技术栈
- 语言：C++20 主框架，Fortran wrapper 调用 WRF 物理。
- 并行：OpenMP CPU reference；后续 CUDA kernel 替换核心算子。
- 构建：CMake。
- I/O：NetCDF-C / NetCDF-Fortran，读取 `wrfinput_d01`、`wrfbdy_d01`、`wrffdda_d01`，写 WRF-compatible core `wrfout`。
- 测试：C++ 单元测试 + Python/pytest 回归对比脚本。
- 平台：优先 x99-wg Linux；设计必须兼容后续 9950X3D WSL2 和 CUDA GPU。
- 外部依赖优先复用现有 PGWRF/WRF 编译环境，不在 v1 引入复杂新框架。

## 初始目录结构
```text
/home/zzy/Projects/TyWRF-Core/
  CMakeLists.txt
  README.md
  docs/
    architecture.md
    wrf_compatibility.md
    cuda_ready_layout.md
    validation_plan.md
  include/tywrf/
    field_view.hpp
    grid.hpp
    state.hpp
    namelist.hpp
    kernel.hpp
  src/
    core/
    dynamics/
    physics_bridge/
    nest/
    io/
    diagnostics/
  bindings/
    wrf_physics/
  tests/
    unit/
    regression/
  tools/
    compare_wrfout.py
    extract_reference.py
    run_6h_cycle_test.py
  third_party/
```

## 功能边界
- v1 支持当前 PGWRF 子集：
  - d01 `10 km`；
  - d02 `2 km` moving nest；
  - 60 eta levels；
  - d01 lateral boundary；
  - d01 spectral nudging；
  - two-way nesting feedback；
  - 6 小时输出。
- v1 输入：
  - `namelist.input`
  - `wrfinput_d01`
  - `wrfbdy_d01`
  - `wrffdda_d01`
- v1 输出 WRF-compatible core `wrfout`，至少包含：
  - `Times, XLAT, XLONG, HGT`
  - `U, V, W, PH, PHB, T, MU, MUB, P, PB`
  - `QVAPOR, QCLOUD, QRAIN, QICE, QSNOW, QGRAUP, QNICE, QNRAIN`
  - `PSFC, U10, V10, T2, Q2`
  - `RAINC, RAINNC`
- v1 物理：
  - 先通过 wrapper 调用 WRF 现有 Thompson、RRTMG、YSU、MM5 surface layer、Noah LSM、KF。
  - 不重写物理参数化。
- v1 不支持：
  - 任意 namelist；
  - 完整 220 个 `wrfout` 变量；
  - bitwise 复现 WRF；
  - best-track nudging；
  - 直接 GPU 版；
  - 非当前 PGWRF 配置族。

## CUDA-Ready 设计约束
- 核心状态采用 SoA 连续数组，不使用 AoS。
- 所有三维场统一扁平布局：
  ```text
  idx = ((j * nz) + k) * nx + i
  ```
  其中 `i` 为连续维。
- 所有场显式 halo，CPU 和未来 CUDA 共用同一 indexer。
- 核心 kernel 只接收 POD view：
  ```text
  pointer + shape + stride + halo
  ```
- 核心积分循环禁止：
  - 嵌套 `std::vector`
  - 虚函数调度
  - NetCDF I/O
  - 日志输出
  - host-only 容器依赖
- CPU OpenMP loop 顺序保持 `j-k-i` 或 tiled `j/k`，`i` 最内层。
- 每个动力算子拆成可替换 kernel：
  - advection
  - pressure gradient
  - acoustic small step
  - mass update
  - diffusion
  - boundary update
  - spectral nudging
  - nest interpolation
  - nest feedback

## 子 Agent 分工
- Agent A：WRF I/O 与 schema
  - 负责 `src/io/`、`tools/extract_reference.py`
  - 输出最小变量 schema、NetCDF 读写接口、WRF-compatible metadata 方案。
- Agent B：数据布局与 kernel API
  - 负责 `include/tywrf/field_view.hpp`、`grid.hpp`、`state.hpp`
  - 输出 CUDA-ready field layout、halo、indexer、OpenMP loop 模板。
- Agent C：动力核心 skeleton
  - 负责 `src/dynamics/`
  - 实现 kernel 调度框架、时间步循环 skeleton、空 tendency 回归测试。
- Agent D：WRF physics bridge
  - 负责 `bindings/wrf_physics/`、`src/physics_bridge/`
  - 确认 WRF physics 调用入口、状态变量 staging buffer、Fortran/C++ ABI。
- Agent E：moving nest 与 nudging
  - 负责 `src/nest/`
  - 实现 d02 moving nest 位置更新接口、parent-child interpolation/feedback skeleton、`wrffdda` 读取和插值接口。
- Agent F：验证体系
  - 负责 `tests/regression/`、`tools/compare_wrfout.py`
  - 实现核心变量 RMSE、台风中心、MSLP、Vmax、降水对比报告。

## 实施阶段
1. Bootstrap
   - 创建 `/home/zzy/Projects/TyWRF-Core`
   - 初始化 Git、CMake、README、基础目录。
   - 建立 build 命令和空测试。
2. WRF Compatibility Baseline
   - 读取当前 KROSA `wrfinput/wrfbdy/wrffdda`。
   - 写出最小 core `wrfout`。
   - 确认现有 Python 后处理能读取。
3. Data Layout Freeze
   - 冻结 SoA field layout、halo、indexer、kernel API。
   - 写 CPU scalar 与 OpenMP kernel smoke tests。
4. Dynamics Skeleton
   - 建立 d01/d02 时间步循环、parent-child step ratio、boundary/nudging 调用点。
   - 暂用 zero tendency 验证 I/O、循环和输出时间轴。
5. Physics Bridge
   - 接入 WRF Fortran 物理 wrapper。
   - 完成 staging buffer 和状态回写。
6. First 6h Cycle
   - 从 WRF reference 初始场跑 6 小时。
   - 输出 d01/d02 core `wrfout`。
   - 与 WRF reference 做 RMSE 和台风诊断对比。
7. Optimization Preparation
   - 记录 kernel timing。
   - 标记 CUDA 迁移优先级。
   - 建立 CPU baseline 性能报告。

## 验收标准
- 编译：
  - `cmake --build build -j` 成功。
  - OpenMP CPU 版可运行。
- I/O：
  - 能读取当前 PGWRF KROSA 输入文件。
  - 能写出 WRF-compatible core `wrfout_d01/d02`。
- 6 小时循环核心场：
  - `U/V/T/PH/MU/P/QVAPOR` normalized RMSE ≤ `5%`。
  - `W` 和水凝物 normalized RMSE ≤ `10%`。
  - `PSFC/U10/V10/T2/Q2` normalized RMSE ≤ `10%`。
- 台风诊断：
  - d02 中心位置误差 ≤ `20 km`。
  - MSLP 误差 ≤ `5 hPa`。
  - Vmax 误差 ≤ `5 m s-1`。
- 工程约束：
  - 核心 kernel 无 NetCDF/I/O/logging。
  - 核心数组连续内存。
  - `i` 维连续访问。
  - CPU/OpenMP kernel API 可直接映射 CUDA kernel 参数。

## 首次执行命令
```bash
ssh x99-wg 'mkdir -p /home/zzy/Projects/TyWRF-Core'
ssh x99-wg 'cd /home/zzy/Projects/TyWRF-Core && git init'
```

## 默认假设
- 项目名固定为 `TyWRF-Core`。
- 目录名使用大小写一致的 `/home/zzy/Projects/TyWRF-Core`。
- v1 以当前 PGWRF/KROSA 配置为唯一目标。
- v1 先 CPU，CUDA-ready，不直接写 CUDA kernel。
- 我作为主 agent 统一协调所有子 agent，所有对子 agent 的调度、验收和集成都由我完成。
