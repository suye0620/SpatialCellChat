# SparseChatArray 测试总结

## 1. 测试范围

`test-SparseChatArray.R` 主要验证 `SparseChatArray` 的基础数据结构和方法，包括：

- 构造器及对象类型
- 维度和维度名称
- `[` 与 `[[` 索引
- `print()`
- `marginSums()`
- `as.data.frame()`
- 不支持的算术运算是否正确报错
- `t()` 逐 layer 转置

## 2. 性能测试设置

性能测试针对 `marginSums(arr, margin = 3L)`，即将所有 layer 沿第三维求和。

| 项目 | 设置 |
|---|---|
| 单层矩阵维度 | `10000 x 10000` |
| layer 数量 | `200` |
| 稀疏度 | `99.99%` |
| 非零密度 | `0.01%` |
| 每层非零元数量 | 约 `10000` |
| 随机种子 | `42` |
| benchmark 重复次数 | `3` |
| 计时工具 | `microbenchmark` |

测试数据在 benchmark 函数外部生成，随后以已经构造好的 `SparseChatArray` 对象作为参数传入：

```r
benchmark_result <- run_sparse_chat_benchmark(
  arr = benchmark_arr,
  iterations = benchmark_config$iterations
)
```

## 3. 比较对象

测试始终比较以下两种实现：

```r
marginSums(arr, margin = 3L)
```

以及 R 原生逐层归约：

```r
Reduce(`+`, unclass(arr))
```

其中：

- `marginSums()` 使用 `cpp_sum_layers` C++ 实现；
- `Reduce(+)` 作为 R 实现对照；
- 两种实现都会使用 `microbenchmark` 重复计时；
- 测试会检查两种实现的结果是否一致。

## 4. 记录的结果指标

benchmark 返回并打印以下信息：

- 输入矩阵维度
- layer 数量
- 实际稀疏度和密度
- 各 layer 非零元数量
- 输入总非零元数量
- 输出非零元数量
- 输入对象大小
- 输出矩阵大小
- C++ `marginSums` 的计时结果
- R `Reduce(+)` 的计时结果
- 两种结果是否一致
- 数据生成阶段的 `generation_time`

## 5. 实际性能结果

本次 benchmark 使用 `3` 次重复测量，结果如下。

### 输入与输出规模

| 指标 | 结果 |
|---|---:|
| 实际密度 | `9.999e-05` |
| 实际稀疏度 | `0.99990001`，约 `99.99%` |
| 每层非零元数量 | `9999`，共 `200` 层 |
| 输入总非零元数量 | `1,999,800` |
| 输出非零元数量 | `1,979,965` |
| 输入对象大小 | `30.81792 MB` |
| 输出矩阵大小 | `22.69849 MB` |
| 数据生成耗时 | `0.42 s` elapsed |

### `marginSums()` C++ 实现

```text
Unit: milliseconds
min      102.6307 ms
lq       105.1252 ms
mean     248.7493 ms
median   107.6197 ms
uq       321.8085 ms
max      535.9974 ms
```

### R `Reduce(+)` 实现

```text
Unit: seconds
min       10.80463 s
lq        11.44451 s
mean      12.22149 s
median    12.08438 s
uq        12.92992 s
max       13.77547 s
```

### 性能比较

以中位数作为主要比较指标：

- `marginSums()`：`107.6197 ms`
- `Reduce(+)`：`12.08438 s`
- C++ 实现约为 R `Reduce(+)` 的 **112.29 倍**速度

以均值计算，C++ 实现约为 R 实现的 **49.13 倍**速度。C++ 均值受到第 3 次测量较高耗时（`535.9974 ms`）影响，因此中位数更适合作为本次 3 次重复测试的代表值。

```r
results_equal
# [1] TRUE
```

两种实现的结果一致，说明在本次输入下，C++ 跨 layer 求和实现与 R `Reduce(+)` 的数值结果相同。

## 6. 当前确认结论

- 测试数据生成和性能测试过程已经分离；benchmark 函数只接收 `arr`，不在函数体内创建测试数据。
- 比较逻辑已固定，不再设置 `compare_reduce` 开关。
- `marginSums()` 和 `Reduce(+)` 的结果会进行一致性检查。
- `t.SparseChatArray()` 已使用 `pbapply::pblapply()` 对各 layer 执行转置，并保留 layer 顺序和名称。
- 项目 `renv` 中的 `microbenchmark` 和 `pbapply` 均可加载。

## 7. 运行状态说明

本节性能结果来自用户提供的实际 `microbenchmark` 输出。

在当前环境中直接重新运行完整测试脚本时，基础检查曾在 `marginSums(arr, 3)` 处因 `cpp_sum_layers` 未注册而中止；随后单独尝试通过 `Rcpp::sourceCpp()` 编译 `src/SpatialChat_Rcpp.cpp` 时，C++ 编译本身已启动并生成目标文件，但被 WorkBuddy 环境对临时文件 `tmp.def` 的清理失败阻断。

因此：

- 本文档中的性能数字来自已提供的实际 benchmark 结果；
- 当前环境的复跑问题属于 Rcpp 动态编译和临时文件清理流程；
- 该环境问题不能解释为 `SparseChatArray` 算法性能或结果错误。

## 8. 总结

`SparseChatArray` 的基础接口和 `marginSums` 性能测试框架已经确定。在本次 `10000 x 10000 x 200` 的高稀疏输入上，C++ `marginSums()` 与 R `Reduce(+)` 结果一致；按中位数计，C++ 实现约快 `112.29` 倍，按均值计约快 `49.13` 倍。

本结果说明，当前 `cpp_sum_layers` 实现适合用于大量稀疏 layer 的跨第三维求和。后续如需形成更稳定的性能结论，建议增加重复次数，并测试不同矩阵维度、layer 数量和稀疏度。