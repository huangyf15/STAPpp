#pragma once
#include <set>
#include <vector>

#include "SparseMatrix.h"

template <typename T> class CSRMatrix : public SparseMatrix<T>
{
private:
    unsigned size;
    unsigned elementCount;
    T* values;
    unsigned* columns;
    unsigned* rowIndexs;
    std::set<unsigned>* _tempColumns;

public:
    CSRMatrix(unsigned m) : SparseMatrix<T>(m)
    {
        size = m;
        values = nullptr;
        columns = nullptr;
        rowIndexs = new unsigned[size + 1];
        _tempColumns = nullptr;
    }

    // call beginPositionMark() will allocate _tempColumns,
    // so that markPostion() can be called.
    void beginPostionMark()
    {
        // _tempColumns = std::set[size]
        _tempColumns = new std::set<unsigned>[size] {};
    }

    void markPosition(unsigned row, unsigned column)
    {
        // insert column
        _tempColumns[row - 1].insert(column);
    }

    void allocate()
    {
        // first, calculate row indexs.
        // notice row index starts from 1
        rowIndexs[0] = 1;
        for (unsigned row = 0; row < size; ++row)
        {
            rowIndexs[row + 1] = rowIndexs[row] + _tempColumns[row].size();
        }

        // allocate columns
        unsigned elementCount = rowIndexs[size] - rowIndexs[0];
        columns = new unsigned[elementCount];

        // write columns
        unsigned count = 0;
        for (unsigned row = 0; row < size; ++row)
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
    }

    // get item at (row, column)
    T& operator()(unsigned row, unsigned column)
    {
        unsigned offset1 = rowIndexs[row - 1] - 1;
        unsigned offset2 = rowIndexs[row] - 1;
        // index lies in [offset1, offset2)
        // find index by bisection

        while (offset2 != (offset1 + 1))
        {
            unsigned offset = (offset1 + offset2) / 2;
            if (columns[offset] > column)
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
};
