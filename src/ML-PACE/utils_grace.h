//
// Created by Yury Lysogorskiy on 30.12.25.
//

#ifndef LAMMPS_GRACE_UTILS_H
#define LAMMPS_GRACE_UTILS_H

#ifndef NO_GRACE_TF

#include <vector>
#include <set>
#include <map>
#include <string>
#include <cstring>
#include <cmath>
#include <cppflow/tensor.h>
#include <tensorflow/c/c_api.h>

#include "lammps.h"

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
                    int extra = static_cast<int>(std::round(real_count * padding_fraction));
                    current_padded_size = real_count +  std::max(extra, 1);
                    padding_history.insert(current_padded_size);
                    was_updated = true;
                } else {
                    current_padded_size = *it;
                    was_updated = false;

                    if (it == padding_history.begin() && can_reduce()) {
                        double waste = (real_count > 0) ?
                                       static_cast<double>(current_padded_size - real_count) / real_count : 0;

                        if (waste > reduction_threshold_fraction) {
                            current_padded_size = real_count;
                            padding_history.insert(current_padded_size);
                            num_of_reductions++;
                            was_updated = true;
                        }
                    }
                }
                n_fake = current_padded_size - real_count;
                return current_padded_size;
            }

            bool last_update_triggered_resize() const { return was_updated; }

            void reset() {
                padding_history.clear();
                num_of_reductions = 0;
                current_padded_size = 0;
            }

        private:
            std::set<int> padding_history;
            bool was_updated = false;
            bool can_reduce() const {
                return (max_reductions == -1 || num_of_reductions < max_reductions);
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
}
#endif
#endif