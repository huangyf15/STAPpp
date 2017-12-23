#pragma once

template <typename T> class SparseMatrix
{
private:
    SparseMatrix(SparseMatrix<T>&);

public:
    SparseMatrix(unsigned size) {}

    virtual T& operator()(unsigned i, unsigned j) = 0;
};
