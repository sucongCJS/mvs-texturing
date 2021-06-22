#include <algorithm>
#include <cassert>
#include <fstream>
#include <cstring>
#include <cerrno>
#include <cmath>

#include <util/file_system.h>
#include <util/exception.h>

#include "histogram.h"

Histogram::Histogram(float _min, float _max, std::size_t num_bins)
    :min(_min), max(_max), num_values(0) {
    bins.resize(num_bins);
}

void Histogram::add_value(float value) {
    float clamped_value = std::max(min, std::min(max, value));  // 让value在min, max之间
    std::size_t index = floor(((clamped_value - min) / (max - min)) * (bins.size() - 1));
    assert(index < bins.size());
    bins[index]++;
    ++num_values;
}

/*
* 返回前百分之percentile占到了min~max之间的多大比重, 如果质量都比较低的话, 返回upper_bound会比较小  
* */
float Histogram::get_approx_percentile(float percentile) const {
    assert(percentile >= 0.0f && percentile <= 1.0f);

    int num = 0;
    float upper_bound = min;
    for (std::size_t i = 0; i < bins.size(); ++i) {
        if (static_cast<float>(num) / num_values > percentile)  // 统计的num超过一定百分比, 返回
            return upper_bound;

        num += bins[i];  // bin[i] 不一定为1
        upper_bound = (static_cast<float>(i) / (bins.size() - 1)) * (max - min) + min;  // 
    }
    return max;
}
