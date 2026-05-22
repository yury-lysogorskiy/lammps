//
// Created by Yury Lysogorskiy on 30.12.25.
//

#ifndef LAMMPS_GRACE_UTILS_H
#define LAMMPS_GRACE_UTILS_H

#include <string>
#include <vector>
#include "lammps.h"

namespace GRACE {
    void log_perf_stats(LAMMPS_NS::LAMMPS *lmp, const std::string &style_name,
                        double total_atoms, long long int total_calls,
                        const std::vector<std::pair<std::string, double>> &timers);

    // Emits a one-shot warning on rank 0 when padding_fraction >= reduction_threshold,
    // a config that would otherwise cause per-step graph recompiles (reduction is
    // suppressed by GracePaddingDimension::update() in that regime).
    void warn_padding_reduction_disabled(LAMMPS_NS::LAMMPS *lmp,
                                         double padding_fraction,
                                         double reduction_threshold);
}

#ifndef NO_GRACE_TF

#include <set>
#include <map>
#include <string>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <cppflow/tensor.h>
#include <tensorflow/c/c_api.h>

namespace GRACE {
    using namespace LAMMPS_NS;

    class GracePaddingDimension {
        public:
            int current_padded_size = 0;
            int num_of_reductions = 0;

            // Settings
            double padding_fraction = 0.01;
            double reduction_threshold_fraction = 0.2;
            int max_reductions = 10;
            bool enabled = true;
            bool verbose = false;
            int n_fake = 0;

            GracePaddingDimension() = default;

            int update(int real_count) {
                if (!enabled) {
                    current_padded_size = real_count;
                    return current_padded_size;
                }

                auto it = padding_history.upper_bound(real_count);

                if (it == padding_history.end()) {
                    current_padded_size = padded_for(real_count);
                    padding_history.insert(current_padded_size);
                    was_updated = true;
                    was_reduced = false;
                } else {
                    current_padded_size = *it;
                    was_updated = false;
                    was_reduced = false;

                    // Refuse to shrink unless padding_fraction is strictly smaller than the
                    // reduction threshold. Otherwise the post-reduction waste (= padding_fraction)
                    // would still exceed the threshold and the next tiny drop in real_count
                    // would re-trigger the reduction, causing per-step graph recompiles.
                    if (it == padding_history.begin() && can_reduce()
                        && padding_fraction < reduction_threshold_fraction) {
                        double waste = (real_count > 0) ?
                                       static_cast<double>(current_padded_size - real_count) / real_count : 0;

                        if (waste > reduction_threshold_fraction) {
                            int new_size = padded_for(real_count);
                            if (new_size < current_padded_size) {
                                current_padded_size = new_size;
                                padding_history.insert(current_padded_size);
                                num_of_reductions++;
                                was_updated = true;
                                was_reduced = true;
                            }
                        }
                    }
                }
                n_fake = current_padded_size - real_count;
                return current_padded_size;
            }

            bool last_update_triggered_resize() const { return was_updated; }
            bool last_update_was_reduction() const { return was_reduced; }
            const char *action_str() const { return was_reduced ? "reducing" : "extending"; }

            void reset() {
                padding_history.clear();
                num_of_reductions = 0;
                current_padded_size = 0;
            }

        private:
            std::set<int> padding_history;
            bool was_updated = false;
            bool was_reduced = false;
            bool can_reduce() const {
                return (max_reductions == -1 || num_of_reductions < max_reductions);
            }
            int padded_for(int real_count) const {
                int extra = static_cast<int>(std::round(real_count * padding_fraction));
                return real_count + std::max(extra, 1);
            }
    };

    void print_tf_inputs(const std::vector<std::tuple<std::string, cppflow::tensor>>& inputs,
                            int me, LAMMPS_NS::LAMMPS *lmp, bool python_ready = false);

    void print_f_data(const double *f_data, int me, LAMMPS_NS::LAMMPS *lmp,
      const std::vector<int>& ind_i_vector, const std::vector<int>& ind_j_vector,
      int * tag);

    /**
     * @brief Get a tensor from a pool or create it if it doesn't exist or size mismatch.
     *
     * If the tensor exists and has the correct byte size, it memcpys the data into the existing buffer.
     * This avoids heap allocations in the hot loop of chunked model calls.
     *
     * IMPORTANT: This function uses direct memcpy into TensorFlow tensor buffers for performance.
     * This is safe ONLY for CPU-resident input tensors. The TensorFlow C API guarantees that
     * TF_TensorData() returns valid writable memory for CPU tensors. Do NOT use this function
     * for GPU tensors or tensors that may be relocated to GPU memory - in such cases, the
     * memcpy would write to invalid memory causing silent corruption or crashes.
     */
    template<typename T>
    cppflow::tensor& get_or_create_tensor(
        std::map<std::string, cppflow::tensor>& pool,
        const std::string& key,
        const std::vector<T>& values,
        const std::vector<int64_t>& shape,
        bool force_recreate = false) 
    {
        size_t byte_size = values.size() * sizeof(T);
        auto it = pool.find(key);

        bool needs_create = force_recreate || (it == pool.end());

        if (!needs_create) {
            // Check byte size of existing tensor
            auto existing_tensor = it->second.get_tensor();
            if (TF_TensorByteSize(existing_tensor.get()) != byte_size) {
                needs_create = true;
            }
        }

        if (needs_create) {
            pool[key] = cppflow::tensor(values, shape);
            return pool[key];
        } else {
            // Safe to memcpy into CPU-resident input tensors
            std::memcpy(TF_TensorData(it->second.get_tensor().get()), values.data(), byte_size);
            return it->second;
        }
    }

    /**
     * @brief Create a tensor with specific dtype from a double vector.
     *
     * This is a simpler version of get_or_create_tensor_dtype that doesn't
     * use a pool. It creates a new tensor each time.
     */
    inline cppflow::tensor create_tensor_dtype(
        const std::vector<double>& values,
        const std::vector<int64_t>& shape,
        cppflow::datatype target_dtype)
    {
        if (target_dtype == TF_DOUBLE) {
            return cppflow::tensor(values, shape);
        }
        // Convert to float
        std::vector<float> f_values(values.begin(), values.end());
        return cppflow::tensor(f_values, shape);
    }

    /**
     * @brief Get a tensor from a pool or create it with specific dtype.
     * 
     * Handles conversion from double to float if target_dtype is TF_FLOAT.
     */
    inline cppflow::tensor& get_or_create_tensor_dtype(
        std::map<std::string, cppflow::tensor>& pool,
        const std::string& key,
        const std::vector<double>& values,
        const std::vector<int64_t>& shape,
        cppflow::datatype target_dtype,
        bool force_recreate = false) 
    {
        if (target_dtype == TF_DOUBLE) {
            return get_or_create_tensor(pool, key, values, shape, force_recreate);
        }
        
        // Convert to float
        std::vector<float> f_values(values.begin(), values.end());
        size_t byte_size = f_values.size() * sizeof(float);
        auto it = pool.find(key);
        bool needs_create = force_recreate || (it == pool.end());

        if (!needs_create) {
            auto existing_tensor = it->second.get_tensor();
            if (TF_TensorByteSize(existing_tensor.get()) != byte_size ||
                TF_TensorType(existing_tensor.get()) != TF_FLOAT) {
                needs_create = true;
            }
        }

        if (needs_create) {
            pool[key] = cppflow::tensor(f_values, shape);
            return pool[key];
        } else {
            std::memcpy(TF_TensorData(it->second.get_tensor().get()), f_values.data(), byte_size);
            return it->second;
        }
    }

    /**
     * @brief Copy TF tensor data into a double vector with widening if necessary.
     */
    inline void copy_tensor_to_vector(const cppflow::tensor& t, std::vector<double>& dest, size_t src_offset, size_t dest_offset, size_t elements) {
        auto tens = t.get_tensor();
        auto dtype = TF_TensorType(tens.get());
        
        if (dtype == TF_DOUBLE) {
            const double *data = static_cast<const double *>(TF_TensorData(tens.get()));
            std::copy_n(data + src_offset, elements, &dest[dest_offset]);
        } else if (dtype == TF_FLOAT) {
            const float *data = static_cast<const float *>(TF_TensorData(tens.get()));
            for (size_t i = 0; i < elements; ++i) dest[dest_offset + i] = static_cast<double>(data[src_offset + i]);
        }
    }
}
#endif
#endif