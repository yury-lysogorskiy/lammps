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
#include <string>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <cppflow/tensor.h>

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
}
#endif
#endif