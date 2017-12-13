#pragma once

template <typename T> class SparseMatrix
{
public:
    SparseMatrix(unsigned size) {}

    virtual T& operator()(unsigned i, unsigned j) = 0;
};
