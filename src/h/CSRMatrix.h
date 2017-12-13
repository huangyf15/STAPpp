#pragma once
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>


#include "SparseMatrix.h"

template <typename T> class CSRMatrix : public SparseMatrix<T>
{
private:
    std::set<std::size_t>* _tempColumns;

public:
    std::size_t size;
    std::size_t elementCount;
    T* values;
    std::size_t* columns;
    std::size_t* rowIndexs;

    CSRMatrix(std::size_t m) : SparseMatrix<T>(m)
    {
        size = m;
        values = nullptr;
        columns = nullptr;
        rowIndexs = new std::size_t[size + 1];
        _tempColumns = nullptr;
    }

    // call beginPositionMark() will allocate _tempColumns,
    // so that markPostion() can be called.
    void beginPostionMark()
    {
        // _tempColumns = std::set[size]
        _tempColumns = new std::set<std::size_t>[size] {};
    }

    void markPosition(std::size_t row, std::size_t column)
    {
        // insert column
        _tempColumns[row - 1].insert(column);
    }

    void allocate()
    {
        // first, calculate row indexs.
        // notice row index starts from 1
        rowIndexs[0] = 1;
        for (std::size_t row = 0; row < size; ++row)
        {
            rowIndexs[row + 1] = rowIndexs[row] + _tempColumns[row].size();
        }

        // allocate columns
        elementCount = rowIndexs[size] - rowIndexs[0];
        columns = new std::size_t[elementCount];

        // write columns
        std::size_t count = 0;
        for (std::size_t row = 0; row < size; ++row)
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
        for (std::size_t i = 0; i < elementCount; ++i)
            values[i] = T(0);
    }

    // get item at (row, column)
    T& operator()(unsigned row, unsigned column)
    {
        if (row > column)
        {
            return this->operator()(column, row);
        }
        std::size_t offset1 = rowIndexs[row - 1] - 1;
        std::size_t offset2 = rowIndexs[row] - 1;
        // index lies in [offset1, offset2)
        // find index by bisection

        while (offset2 != (offset1 + 1))
        {
            std::size_t offset = (offset1 + offset2) / 2;
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

    std::size_t dim() const { return size; }
};

template <typename T> std::ostream& operator<<(std::ostream& out, CSRMatrix<T> mat)
{
    out << "CSR Matrix, size = " << mat.size << std::endl;
    out << "values = " << std::endl << "(";
    for (unsigned i = 0; i < mat.elementCount; ++i)
    {
        out << std::setw(14) << mat.values[i];
    }
    out << ")\ncolumns = \n(";
    for (unsigned i = 0; i < mat.elementCount; ++i)
    {
        out << std::setw(14) << mat.columns[i];
    }
    out << ")\nrowIndexs = \n(";
    for (unsigned i = 0; i <= mat.size; ++i)
    {
        out << std::setw(14) << mat.rowIndexs[i];
    }
    out << ")";
    return out;
}
