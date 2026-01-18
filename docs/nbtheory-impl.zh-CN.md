# OpenFHE 数论工具 nbtheory-impl.h 中文说明

本文件汇总 `src/core/include/math/nbtheory-impl.h` 中各函数/模板的作用、参数与典型使用场景，便于在密码学与同态计算相关代码中正确调用。

说明约定：
- IntType 表示大整数或本地整数类型（如 `BigInteger` 或 `NativeInteger`），需支持文中所用的基本算术与模运算接口。
- IntVector 表示整数向量容器，其元素类型为 `IntVector::Integer`。
- usint 为无符号整型（库中常用别名）。

目录：
- 随机、素性测试与分解
- 生成元与原根函数
- 互素与欧拉函数相关
- 取整/分圆多项式与多项式工具
- 便捷素数搜索器

---

- 函数：`template <typename IntType> static IntType RNG(const IntType& modulus)`
  - 作用：在区间 [0, modulus) 内生成均匀随机数（分块拼接，避免偏差）。
  - 参数：modulus 上界（不含）。
  - 返回：随机整数。
  - 使用场景：Miller-Rabin 见证、Pollard-Rho、随机生成元候选等。
  - 注意：要求 `modulus > 0`；内部使用 PRNG 与 32 位分块方式。

- 函数：`template <typename IntType> static bool WitnessFunction(const IntType& a, const IntType& d, usint s, const IntType& p)`
  - 作用：Miller–Rabin 素性测试中的单次见证计算。
  - 参数：a 见证；d 与 s 满足 p-1 = 2^s · d；p 被测数。
  - 返回：true 表示“合数”见证；false 表示“可能为素数”。
  - 使用场景：`MillerRabinPrimalityTest` 内部调用。

- 函数：`template <typename IntType> static IntType FindGenerator(const IntType& q)`
  - 作用：在素数模 q 的乘法群 Z_q^* 中找一个原根（生成元）。
  - 参数：q 为素数。
  - 返回：生成元 g。
  - 使用场景：构造 NTT 原根、原始 m 次单位根等。
  - 注意：通过分解 q-1 并检验候选；依赖 `PrimeFactorize` 与随机数。

- 函数：`template <typename IntType> IntType FindGeneratorCyclic(const IntType& q)`
  - 作用：在任意循环群（不一定是素数模）上寻找生成元。
  - 参数：q 群的模数/规模（可组合数）。
  - 返回：生成元 g。
  - 使用场景：任意 cyclotomic 设置下的生成元查找。
  - 注意：使用 φ(q) 与其素因子约束判断，且保证候选与 q 互素。

- 函数：`template <typename IntType> bool IsGenerator(const IntType& g, const IntType& q)`
  - 作用：判断 g 是否为模 q 的生成元（循环群）。
  - 参数：g 候选，q 模。
  - 返回：true 表示 g 为生成元。
  - 使用场景：验证给定原根或生成元正确性。

- 函数：`template <typename IntType> IntType RootOfUnity(usint m, const IntType& modulo)`
  - 作用：返回模 `modulo` 的最小 m 次原始单位根（primitive root of unity）。
  - 参数：m 循环阶（通常为 2 的幂），modulo 素数模（要求 (modulo-1) 可被 m 整除）。
  - 返回：单位根 r，满足 r^m ≡ 1 (mod modulo) 且为原始根。
  - 使用场景：NTT/FFT、RNS-NTT 同态运算准备。
  - 注意：
    - 先找生成元，再指数降幂得到 m 次原根；
    - 对所有与 m 互素的幂次遍历，选取最小原根，保证不同上下文一致性；
    - 若不满足 (modulo-1) 可被 m 整除，抛出异常。

- 函数：`template <typename IntType> std::vector<IntType> RootsOfUnity(usint m, const std::vector<IntType>& moduli)`
  - 作用：对一组 `moduli` 分别计算它们的 m 次原根。
  - 使用场景：DCRT 多塔 NTT 根预计算。

- 函数：`template <typename IntType> IntType GreatestCommonDivisor(const IntType& a, const IntType& b)`
  - 作用：扩展欧几里得算法求 gcd(a,b)。
  - 使用场景：判互素、Pollard-Rho 辅助。

- 函数：`template <typename IntType> bool MillerRabinPrimalityTest(const IntType& p, const usint niter)`
  - 作用：Miller–Rabin 概率素性测试。
  - 参数：p 被测整数；niter 迭代次数（见证数）。
  - 返回：true 认为素数；false 非素。
  - 使用场景：素数搜索、NTT 质数校验。

- 函数：`template <typename IntType> const IntType PollardRhoFactorization(const IntType& n)`
  - 作用：返回 n 的一个非平凡因子（Pollard Rho）。
  - 参数：n 正整数。
  - 返回：非平凡因子（若 n 为偶数先返回 2）。
  - 使用场景：素因子分解、FindGenerator 等辅助。

- 函数：`template <typename IntType> void PrimeFactorize(IntType n, std::set<IntType>& primeFactors)`
  - 作用：递归分解 n，输出其所有素因子（无重）。
  - 使用场景：生成元/原根判定、群阶分解等。

- 函数：`template <typename IntType> IntType FirstPrime(uint32_t nBits, uint64_t m)`
  - 作用：找第一个满足 bit 长度为 `nBits` 且 `q ≡ 1 (mod m)` 的素数 q（从 2^nBits 往上偏移对齐 m）。
  - 参数：nBits 位数；m 模约束。
  - 返回：满足条件的最小素数。
  - 使用场景：构造 NTT 友好素数（通常令 m=2n）。
  - 注意：若超出 `NativeInteger` 的最大位长上限会抛出异常。

- 函数：`template <typename IntType> IntType LastPrime(uint32_t nBits, uint64_t m)`
  - 作用：找最后一个满足 bit 长度为 `nBits` 且 `q ≡ 1 (mod m)` 的素数 q（向下搜索）。
  - 注意：返回值 MSB 必须等于 nBits，否则抛出异常。

- 函数：`template <typename IntType> IntType NextPrime(const IntType& q, uint64_t m)`
  - 作用：找到下一位满足 `≡ 1 (mod m)` 的素数（从 q+m 开始向上找）。
  - 使用场景：构造模数链、DCRT 多塔质数序列。

- 函数：`template <typename IntType> IntType PreviousPrime(const IntType& q, uint64_t m)`
  - 作用：找到上一位满足 `≡ 1 (mod m)` 的素数（从 q-m 开始向下找）。
  - 使用场景：回退/补齐某一层的 NTT 质数。

- 函数：`template <typename IntType> IntType NextPowerOfTwo(IntType n)`
  - 作用：返回不小于 n 的二的幂的幂指数（ceil(log2(n)))。
  - 使用场景：计算 NTT 尺寸、指定位宽。

- 函数：`template <typename IntType> std::vector<IntType> GetTotientList(const IntType& n)`
  - 作用：返回 [1, n) 中与 n 互素的所有整数列表。
  - 使用场景：遍历 m 次单位根所有与 m 互素的幂以筛原始单位根；一般性 totient 列表。

- 函数：`template <typename IntVector> IntVector PolyMod(const IntVector& dividend, const IntVector& divisor, const typename IntVector::Integer& modulus)`
  - 作用：多项式模约简：`dividend mod divisor`（系数模 modulus）。
  - 使用场景：多项式环上操作的基础工具（如模 Φ_m(x)）。

- 函数：`template <typename IntVector> IntVector PolynomialMultiplication(const IntVector& a, const IntVector& b)`
  - 作用：朴素多项式乘法（系数模 a.getModulus()）。
  - 使用场景：小规模测试或构造型函数；大规模通常用 NTT。

- 函数：`template <typename IntVector> IntVector GetCyclotomicPolynomial(usint m, const typename IntVector::Integer& modulus)`
  - 作用：返回 m 次分圆多项式 Φ_m(x) 在系数域 modulus 下的表示（系数被映射为模数域元素）。
  - 使用场景：构造环 `Z_q[x]/(Φ_m(x))` 或 `Z_q[x]/(x^n+1)`（当 m=2n）。

- 函数：`template <typename IntVector> typename IntVector::Integer SyntheticRemainder(const IntVector& dividend, const typename IntVector::Integer& a, const typename IntVector::Integer& modulus)`
  - 作用：综合除法（Horner 法）求 `dividend(x)` 在 `x=a` 处的值（系数模 modulus）。
  - 使用场景：快速求值、内插/检验。

- 函数：`template <typename IntVector> IntVector SyntheticPolyRemainder(const IntVector& dividend, const IntVector& aList, const typename IntVector::Integer& modulus)`
  - 作用：对一组点 `aList` 批量求值（返回对应余数列表）。
  - 使用场景：批量点值计算。

- 函数：`template <typename IntVector> IntVector PolynomialPower(const IntVector& input, usint power)`
  - 作用：将多项式做“幂次稀疏化”映射：把每个系数移至 `i*power` 的位置，等价于 `x -> x^power` 的替换效果（非乘方乘法）。
  - 使用场景：构造型操作、指数映射（如升维/采样中的结构化多项式）。

- 函数：`template <typename IntVector> IntVector SyntheticPolynomialDivision(const IntVector& dividend, const typename IntVector::Integer& a, const typename IntVector::Integer& modulus)`
  - 作用：综合除法：对 `(x - a)` 进行多项式除法，返回商（系数模 modulus）。
  - 使用场景：快速多项式除法（一次线性因子）。

---

使用建议与常见场景
- NTT 场景：
  - 使用 `FirstPrime/NextPrime` 构造满足 `q ≡ 1 (mod 2n)` 的质数链；
  - 通过 `RootOfUnity(2n, q)` 获取 2n 次原根（或 n 次原根视实现），再建立 NTT 表；
  - 在 DCRT 多塔下，对每个 qi 单独构建原根，组合为多塔 NTT。

- 素性与分解：
  - `MillerRabinPrimalityTest` 较快判断素性；
  - 需要分解群阶或检验生成元时，用 `PrimeFactorize` 配合 `FindGenerator/IsGenerator`。

- 多项式工具：
  - `GetCyclotomicPolynomial` 用于构造环模；
  - `PolyMod`、`SyntheticRemainder`、`SyntheticPolynomialDivision` 支持基础代数运算与求值。

注意事项
- `RootOfUnity` 要求 `(q-1) % m == 0`，否则抛异常；
- 在 `NativeInteger` 情况下，`FirstPrime/LastPrime` 对位长有硬性上限；
- 多数函数依赖 IntType 的接口（如 `ModExp / Mod / RShiftEq / ComputeMu` 等），需确保所选整数类型实现了这些接口；
- 随机函数 `RNG` 基于库内 PRNG；在参数生成或安全性敏感场景，请使用库内一致的安全随机源。

示例片段（思路）
- 构造 NTT 质数与原根
  - 给定 n=2^k，令 m=2n，调用 `FirstPrime(bits, m)` 得到 q；
  - 调 `RootOfUnity(m, q)` 获取 2n 次原根；
  - 据此构造 NTT/INTT 所需的 twiddle 因子。

- 生成元验证
  - 对给定 q，先分解 q-1：`PrimeFactorize(q-1, S)`；
  - 随机挑 g，用 `IsGenerator(g, q)` 验证；不通过则重试。

- 多项式求值
  - 使用 `SyntheticRemainder(P, a, mod)` 快速求 P(a)；对多点用 `SyntheticPolyRemainder`。

本说明仅针对 nbtheory-impl.h 中公开的工具函数，具体类型实现（如 IntType 的大整数接口）请参考相应头文件与实现。