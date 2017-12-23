#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <functional>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using std::size_t;

template <typename T> class Matrix
{
public:
    typedef size_t Pos_t;
    // typedef T (*Func_t)(Pos_t, Pos_t);
    typedef std::function<T(Pos_t, Pos_t)> Func_t;

    Matrix()
    {
        _elements = nullptr;
        setSize(0, 0);
    };

    Matrix(Pos_t rows, Pos_t columns)
    {
        setSize(rows, columns);
        allocate();
    };

    Matrix(const Matrix<T>& m) { this->operator=(m); }

    Matrix(const std::vector<T>& v) { this->operator=(v); }

    Matrix(Pos_t rows, Pos_t columns, const T* elements)
    {
        setSize(rows, columns);
        allocate();
        std::memcpy(_elements, elements, rows * columns * sizeof(T));
    }

    Matrix(Pos_t rows, Pos_t columns, Func_t func)
    {
        setSize(rows, columns);
        allocate();
        for (Pos_t j = 1; j <= _columns; ++j)
        {
            for (Pos_t i = 1; i <= _rows; ++i)
            {
                at(i, j) = func(i, j);
            }
        }
    }

    std::string c_str(const char* singleElementFormat = " %10.3g") const
    {
        if (!_elements)
        {
            return std::string("(NULL)");
        }
        std::string res = "";

        for (Pos_t i = 1; i <= _rows; ++i)
        {
            std::string thisLineStr;
            char startCode = 0, endCode = 0;
            if (1 == _rows)
            {
                startCode = '(';
                endCode = ')';
            }
            else if (1 == i)
            { // first row
                startCode = '/';
                endCode = '\\';
            }
            else if (i == _rows)
            { // last row
                startCode = '\\';
                endCode = '/';
            }
            else
            {
                startCode = endCode = '|';
            }

            thisLineStr = startCode;

            for (Pos_t j = 1; j <= _columns; ++j)
            {
                char temp[128]{0};
                sprintf(temp, singleElementFormat, c_at(i, j));
                thisLineStr += temp;
            }

            thisLineStr += endCode;
            if (i != _rows)
            {
                thisLineStr += '\n';
            }
            res += thisLineStr;
        }
        return res;
    }

    const T& c_at(Pos_t row, Pos_t column) const { return at(row, column); }

    void set(Pos_t row, Pos_t column, T value) { at(row, column) = value; }

    Matrix<T> operator*(T val) const
    {
        Matrix<T> res(*this);
        for (Pos_t i = 0; i < _rows * _columns; i++)
        {
            res._elements[i] *= val;
        }
        return res;
    }

    T& at(Pos_t row, Pos_t column) const
    {
#ifdef _MATRIX_DEBUG_
        if (row <= 0 || row > _rows || column <= 0 || column > _columns)
        {
            char buff[200];
            sprintf(buff, "Index error happened at method at() with row=%d, column=%d, "
                          "this->rows=%d, this->column=%d",
                    row, column, _rows, _columns);
            throw std::runtime_error(buff);
        }
#endif
        return _elements[row - 1 + (column - 1) * _rows];
    }

    Matrix<T> transpose() const
    {
        Matrix<T> result(_columns, _rows);
        for (Pos_t i = 1; i <= _rows; ++i)
        {
            for (Pos_t j = 1; j <= _columns; ++j)
            {
                result.at(j, i) = c_at(i, j);
            }
        }
        return result;
    }

    Matrix<T> inverseByRow() const
    {
        if (_rows != _columns)
            throw std::logic_error("applying inverse to a non-square matrix.");

        // copy self
        Pos_t n = _rows;
        Matrix<T> ori(*this);
        Matrix<T> res(n, n);
        for (Pos_t i = 1; i <= n; ++i)
            res.at(i, i) = T(1);

        for (Pos_t i = 1; i <= n; ++i)
        {
            Pos_t mark = i;
            T mval = std::abs(ori.c_at(i, i));
            for (Pos_t j = i + 1; j <= n; ++j)
            { // find max element in column
                if (std::abs(ori.c_at(j, i)) > mval)
                {
                    mval = std::abs(ori.c_at(j, i));
                    mark = j;
                }
            }

            if (i != mark)
            {
                ori.swapRow(i, mark);
                res.swapRow(i, mark);
            }
            res.rowTimes(i, T(1) / ori.c_at(i, i));
            ori.rowTimes(i, T(1) / ori.c_at(i, i));

            for (Pos_t j = 1; j <= n; ++j)
            {
                if (j != i)
                {
                    res.addRowTo(j, i, -ori.c_at(j, i));
                    ori.addRowTo(j, i, -ori.c_at(j, i));
                }
            }
        }

        return res;
    }

    Matrix<T>& swapRow(Pos_t i, Pos_t j)
    {
        T buff;
        for (Pos_t n = 1; n <= _columns; ++n)
        {
            buff = c_at(i, n);
            at(i, n) = c_at(j, n);
            at(j, n) = buff;
        }
        return *this;
    }

    Matrix<T>& rowTimes(Pos_t i, T times)
    {
        for (Pos_t n = 1; n <= _columns; ++n)
        {
            at(i, n) *= times;
        }
        return *this;
    }

    Matrix<T>& addRowTo(Pos_t dest, Pos_t src, T times)
    {
        for (Pos_t n = 1; n <= _columns; ++n)
        {
            at(dest, n) += c_at(src, n) * times;
        }
        return *this;
    }

    Matrix<T> inverse() const
    {
        if (_rows != _columns)
            throw std::logic_error("applying inverse to a non-square matrix.");

        Pos_t n = _rows;
        Matrix<T> ori(*this);
        Matrix<T> res(n, n);
        for (Pos_t i = 1; i <= n; ++i)
            res.at(i, i) = T(1);

        for (Pos_t i = 1; i <= n; ++i)
        {
            Pos_t mark = i;
            T mval = std::abs(ori.c_at(i, i));
            for (Pos_t j = i + 1; j <= n; ++j)
            {
                if (std::abs(ori.c_at(i, j)) > mval)
                {
                    mval = std::abs(ori.c_at(i, j));
                    mark = j; // mark column j
                }
            }

            if (i != mark)
            {
                ori.swapColumn(i, mark);
                res.swapColumn(i, mark);
            }
            res.columnTimes(i, T(1) / ori.c_at(i, i));
            ori.columnTimes(i, T(1) / ori.c_at(i, i));

            for (Pos_t j = 1; j <= n; ++j)
            {
                if (j != i)
                {
                    res.addColumnTo(j, i, -ori.c_at(i, j));
                    ori.addColumnTo(j, i, -ori.c_at(i, j));
                }
            }
        }
        return res;
    }

    Matrix<T>& swapColumn(Pos_t i, Pos_t j)
    {
        // swap row i and row j
        T* buff = new T[_rows];
        std::memcpy(buff, _elements + (i - 1) * _rows, _rows * sizeof(T));
        std::memcpy(_elements + (i - 1) * _rows, _elements + (j - 1) * _rows, _rows * sizeof(T));
        std::memcpy(_elements + (j - 1) * _rows, buff, _rows * sizeof(T));
        delete[] buff;
        return *this;
    }

    Matrix<T>& columnTimes(Pos_t i, T times)
    {
        for (Pos_t n = 1; n <= _rows; ++n)
        {
            at(n, i) *= times;
        }
        return *this;
    }

    Matrix<T>& addColumnTo(Pos_t dest, Pos_t src, T times)
    {
        for (Pos_t n = 1; n <= _rows; ++n)
        {
            at(n, dest) += c_at(n, src) * times;
        }
        return *this;
    }

    Matrix<T> partition(Pos_t startRow, Pos_t endRow, Pos_t startColumn, Pos_t endColumn)
    {
        Pos_t rows = endRow - startRow + 1;
        Pos_t columns = endColumn - startColumn + 1;
        Matrix<T> result(rows, columns);

        for (Pos_t i = 1; i <= columns; ++i)
        {
            std::memcpy(result._elements + (i - 1) * rows,
                        _elements + (startColumn + i - 2) * _rows + startRow - 1, sizeof(T) * rows);
        }

        return result;
    }

    Matrix<T> operator*(const Matrix<T>& m)
    {
        if (_columns != m._rows)
        {
            throw std::runtime_error("size don't match for operator *");
        }
        Matrix<T> res(_rows, m._columns);
        for (Pos_t i = 1; i <= res._rows; ++i)
        {
            for (Pos_t j = 1; j <= res._columns; ++j)
            {
                T sum = 0;
                for (Pos_t n = 1; n <= _columns; ++n)
                {
                    // std::cout << "adding this["<<i<<"]["<<n<<"] * m["<<n<<"]["<<j<<"] = "
                    //           <<c_at(i, n) * m.c_at(n, j)<<" to sum="<<sum << std::endl;
                    sum += c_at(i, n) * m.c_at(n, j);
                }
                // std::cout << "res["<<i<<"]["<<j<<"] = "<<sum<<std::endl;
                res.set(i, j, sum);
            }
        }
        return res;
    }

    Matrix<T> operator-(const Matrix<T>& m)
    {
        if (_rows != m._rows || _columns != m._columns)
        {
            throw std::runtime_error("in Matrix<T>::operator-: matrix size don't match.");
        }
        Matrix<T> res(_rows, _columns);

        for (Pos_t i = 0; i < _rows * _columns; ++i)
        {
            res._elements[i] = _elements[i] - m._elements[i];
        }

        return res;
    }

    Matrix<T> operator+(const Matrix<T>& m)
    {
        if (_rows != m._rows || _columns != m._columns)
        {
            throw std::runtime_error("in Matrix<T>::operator+: matrix size don't match.");
        }
        Matrix<T> res(_rows, _columns);

        for (Pos_t i = 0; i < _rows * _columns; ++i)
        {
            res._elements[i] = _elements[i] + m._elements[i];
        }

        return res;
    }

    Matrix<T>& operator+=(const Matrix<T>& m)
    {
        if (_rows != m._rows || _columns != m._columns)
        {
            throw std::runtime_error("in Matrix<T>::operator+=: matrix size don't match.");
        }

        for (Pos_t i = 0; i < _rows * _columns; ++i)
        {
            _elements[i] += m._elements[i];
        }

        return *this;
    }

    Matrix<T>& operator=(const Matrix<T>& m)
    {
        setSize(m._rows, m._columns);
        if (_rows * _columns)
        {
            allocate();
            std::memcpy(_elements, m._elements, _rows * _columns * sizeof(T));
        }
        else
        {
            free();
        }
        return *this;
    }

    Matrix<T>& operator=(const std::vector<T>& v)
    {
        setSize(v.size(), 1);
        if (_rows)
        {
            allocate();
            for (Pos_t i = 1; i <= _rows; ++i)
            {
                at(i, 1) = v[i - 1];
            }
        }
        else
        {
            free();
        }
        return *this;
    }

    ~Matrix()
    {
#ifdef _MATRIX_DEBUG_
        std::clog << "~Matrix()" << std::endl;
#endif
        if (_elements)
        {
            delete[] _elements;
        }
    };

private:
    T* _elements = nullptr;
    Pos_t _rows = 0;
    Pos_t _columns = 0;

    void setSize(Pos_t rows, Pos_t columns)
    {
        _rows = rows;
        _columns = columns;
    }
    void free()
    {
        if (_elements)
        {
            delete[] _elements;
        }
    }
    void allocate()
    {
#ifdef _MATRIX_DEBUG_
        std::clog << "allocating" << std::endl;
#endif
        if (!(_rows * _columns))
        {
            free();
        }
        else
        {
            if (_elements)
                free();
            _elements = new T[_rows * _columns]{T(0)};
#ifdef _MATRIX_DEBUG_
            std::cout << "_elements allocated, size=(" << _rows << ", " << _columns << ")"
                      << std::endl;
#endif
        }
    }
};

template <typename T> std::ostream& operator<<(std::ostream& out, const Matrix<T>& m)
{
    out << m.c_str();
    return out;
}

typedef Matrix<double> DoubleMatrix;
