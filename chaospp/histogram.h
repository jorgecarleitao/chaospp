#ifndef chaospp_histogram_h
#define chaospp_histogram_h

#include <vector>
#include <assert.h>
#include <string>
#include <math.h>

#include "io.h"

namespace histogram {

template <typename T>
class Histogram {
protected:
    T _lowerBound;
    T _upperBound;

    unsigned int _bins;  // number of bins of the histogram
    unsigned int _count; // number of measured samples
    std::vector<unsigned int> _histogram;  // histogram of samples over bins

    inline virtual T v(T value) const {
        return value;
    }
    inline virtual T iv(T value) const {
        return value;
    }

public:
    Histogram(T lowerBound, T upperBound, unsigned int bins) :
            _lowerBound(v(lowerBound)), _upperBound(v(upperBound)), _bins(bins), _histogram(_bins + 1) {
        reset();
    }

    inline unsigned int bins() const {return _bins;}

    inline unsigned int count() const {return _count;}

    unsigned int const& operator[](unsigned int idx) const {
        assert(idx < _histogram.size());
        return _histogram[idx];
    }

    bool invalid_value(T const& value) const {
        return v(value) <= _lowerBound || v(value) >= _upperBound;
    }

    T get_invalid_value() const {
        return iv(_lowerBound) - 1;
    }

    unsigned int bin(T value) const {
        value = v(value);
        if (value <= _lowerBound)
            return 0;
        if (value >= _upperBound)
            return _bins;

        unsigned int bin = (value - _lowerBound)*_bins/(_upperBound - _lowerBound);
        assert(bin < _bins);
        return bin;
    }

    inline double h() const {
        return (_upperBound - _lowerBound)/_bins;
    }

    // inverse of `bin`
    T value(unsigned int bin) const {
        if (bin >= _bins)
            return iv(_upperBound);
        if (bin <= 0)
            return iv(_lowerBound);
        return iv(_lowerBound + (_upperBound - _lowerBound)*bin/_bins);
    }

    virtual void reset() {
        for (unsigned int bin = 0; bin <= _bins; bin++)
            _histogram[bin] = 0;
        _count = 0;
    }

    virtual void add(T value) {
        _histogram[bin(value)]++;
        _count++;
    }

    void print() const {
        unsigned int sum = 0;
        for (unsigned int bin = 0; bin <= _bins; bin++) {
            sum += _histogram[bin];
            if (_histogram[bin] > 0)
                printf("%e %e %d\n", value(bin), _histogram[bin]*1./_count, bin);
        }
        assert(sum == _count);
    }

    virtual void export_pretty(std::string file_name, std::string directory="") const {
        std::vector<std::vector<double> > data;

        for (unsigned int bin = 0; bin <= _bins; bin++)
        {
            if (_histogram[bin] > 0)
            {
                std::vector<double> row(2);
                row[0] = value(bin);
                row[1] = _histogram[bin]*1./_count;
                data.push_back(row);
            }
        }
        io::save(data, directory + file_name);
    }

    virtual void export_histogram(std::string file_name, std::string directory="") const {
        std::vector<std::vector<double> > data;

        for (unsigned int bin = 0; bin <= _bins; bin++)
        {
            if (_histogram[bin] > 0) {
                std::vector<double> row(2);
                row[0] = bin;
                row[1] = _histogram[bin]*1./_count;
                data.push_back(row);
            }
        }
        io::save(data, directory + file_name);
    }
};


template <typename T>
class Log2Histogram : public Histogram<T> {
protected:
    virtual T v(T value) const {
        return log2(value);
    }
    virtual T iv(T value) const {
        return pow(2, value);
    }
public:
    Log2Histogram(T lowerBound,
                  T upperBound,
                  unsigned int bins) : Histogram<T>(v(lowerBound), v(upperBound), bins) {}
};

}

#endif
