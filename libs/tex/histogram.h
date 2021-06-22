#ifndef TEX_HISTOGRAM_HEADER
#define TEX_HISTOGRAM_HEADER

#include <vector>
#include <string>

// 一个直方图，包含固定数量的bin, 
// 计算近似permiles
class Histogram{
    private:
        std::vector<unsigned int> bins;  // 一定范围内的quality会被分到一个bin里
        float min;
        float max;
        int num_values;  // value个数

    public:
        Histogram(float _min, float _max, std::size_t num_bins);

        // 值会被控制在[_min, _max] 之间
        void add_value(float value);

        float get_approx_percentile(float percentile) const;
};

#endif /* TEX_HISTOGRAM_HEADER */
