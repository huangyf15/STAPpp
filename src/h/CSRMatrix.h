#pragma once
#include <iomanip>
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>

#include "SparseMatrix.h"

#define CSR_USE_VECTOR

#ifdef CSR_USE_VECTOR
typedef std::vector<int> STL_t;
#define CSR_OPT push_back
#else
typedef std::set<int> STL_t;
#define CSR_OPT insert
#endif

template <typename T> class CSRMatrix : public SparseMatrix<T>
{
private:
    STL_t* _tempColumns;

public:
    int size;
    int elementCount;
    T* values;
    int* columns;
    int* rowIndexs;

    CSRMatrix(int m) : SparseMatrix<T>(m)
    {
        size = m;
        values = nullptr;
        columns = nullptr;
        rowIndexs = new int[size + 1];
        _tempColumns = nullptr;
    }

    // call beginPositionMark() will allocate _tempColumns,
    // so that markPostion() can be called.
    void beginPostionMark()
    {
        // _tempColumns = std::set[size]

        _tempColumns = new STL_t[size] {};

    }

    void markPosition(int row, int column)
    {
        // insert column
        _tempColumns[row - 1].CSR_OPT(column);
    }

    void allocate()
    {
        #ifdef CSR_USE_VECTOR
        for (int row = 0; row < size; ++row)
        {
            std::vector<int>& v = _tempColumns[row];
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
        }
        #endif
        // first, calculate row indexs. Notice row index starts from 1
        rowIndexs[0] = 1;
        for (int row = 0; row < size; ++row)
        {
            rowIndexs[row + 1] = rowIndexs[row] + int(_tempColumns[row].size());
        }

        // allocate columns
        elementCount = rowIndexs[size] - rowIndexs[0];
        columns = new int[elementCount];

        // write columns
        int count = 0;
        for (int row = 0; row < size; ++row)
        {
            // column already sorted in _tempColumns
            for (const auto& column : _tempColumns[row])
            {
                columns[count] = column;
                count++;
            }
        }

        // free _tempColumns
        delete[] _tempColumns;
        _tempColumns = nullptr;

        // allocate values last to save memory
        values = new T[elementCount];
        for (int i = 0; i < elementCount; ++i)
            values[i] = T(0);
    }

    // get item at (row, column)
    T& operator()(unsigned row, unsigned column)
    {
        if (row > column)
        {
            return this->operator()(column, row);
        }
        int offset1 = rowIndexs[row - 1] - 1;
        int offset2 = rowIndexs[row] - 1;
        // index lies in [offset1, offset2)
        // find index by bisection

        while (offset2 != (offset1 + 1))
        {
            int offset = (offset1 + offset2) / 2;
            if (columns[offset] > int(column))
            {
                offset2 = offset;
            }
            else
            {
                offset1 = offset;
            }
        }

        return values[offset1];
    }

    int dim() const { return size; }

    ~CSRMatrix()
    {
        delete[] rowIndexs;
        delete[] values;
        delete[] columns;
    }
};

template <typename T> std::ostream& operator<<(std::ostream& out, CSRMatrix<T> mat)
{
    out << "CSR Matrix, size = " << mat.size << std::endl;
    out << "values = " << std::endl << "(";
    for (int i = 0; i < mat.elementCount; ++i)
    {
        out << std::setw(14) << mat.values[i];
    }
    out << ")\ncolumns = \n(";
    for (std::size_t i = 0; i < mat.elementCount; ++i)
    {
        out << std::setw(14) << mat.columns[i];
    }
    out << ")\nrowIndexs = \n(";
    for (std::size_t i = 0; i <= mat.size; ++i)
    {
        out << std::setw(14) << mat.rowIndexs[i];
    }
    out << ")";
    return out;
}
