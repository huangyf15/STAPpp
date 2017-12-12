# MKL PARDISO 求解器

update:2017/12/12

添加稀疏求解器有两项主要工作：

+ 用 CSR 矩阵替代 Skyline 矩阵

+ 用 MKL PARDISO 求解题替代 LDLT 求解器

## CSRMatrix

从 SkylineMatrix 替换成 CSRMatrix 矩阵需要进行的工作有：

+ 完成一个`CSRMatrix`类。

+ 在 `CDomain::AllocateMatrices()`里面，更换或重写`CDomain::CalculateColumnHeights()`和`CDomain::CalculateDiagonalAddress()`
`SkylineMatrix::Allocate()`三项。

+ 需要值得注意的是，计算CSR储存方法中的`columns`是比较困难的。
为了知道`columns`，实际需要知道每一行有多少个非零元，而为了去重，
在`AllocateMatrices()`阶段就需要全部储存这些非零元，但此阶段无法按顺序计算，很难用连续内存储存。

+ 我在考虑是否先采用类似 Skyline 的储存方式，先划分一段较大的空间，但是这样就失去了CSR的优点，跟Skyline没差别了，至少峰值没差别。但是如果不这样，就要用类似链表的结构储存，效率可能会很低，对每一行是`n^2`的效率，`n`代表非零元个数。这时需要`NEQ`个链表，先暂时储存`columns`。

+ 我觉得好像还可以……直接用`std::set`来实现储存`columns`，同时计算`rowIndex`。可以考虑先这样实现。反正原来的代码这里也不是特别搞笑。用`sed::set`的话，效率可能差不了多少，因为本来就是有序链表。而且更重要的是，这样可以避免内存碎片化。
