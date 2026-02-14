#ifndef NO_GRACE_TF

//
// Created by Yury Lysogorskiy on 30.12.25.
//
#include "utils_grace.h"

#include <vector>
#include <cppflow/tensor.h>
#include <tensorflow/c/c_api.h>
#include "lammps.h"
#include <iomanip>
#include <sstream>
#include "utils.h"

using namespace LAMMPS_NS;

namespace GRACE {

void print_tf_inputs(const std::vector<std::tuple<std::string, cppflow::tensor>>& inputs,
                            int me, LAMMPS_NS::LAMMPS *lmp, bool python_ready) {
    std::stringstream ss;
    ss << "\n" << std::string(100, '=') << "\n";
    ss << "[GRACE-DEBUG, Proc #" << me << "]: TensorFlow Graph Inputs & Content\n";
    ss << std::string(100, '-') << "\n";

    ss << std::left << std::setw(40) << "Input Name (:port)"
       << " | " << std::setw(10) << "DataType"
       << " | " << std::setw(15) << "Shape" << "\n";
    ss << std::string(70, '-') << "\n";

    for (const auto& item : inputs) {
        const std::string& name = std::get<0>(item);
        const cppflow::tensor& tensor = std::get<1>(item);

        auto shape = tensor.shape().get_data<int64_t>();
        std::string shape_str = "[";
        int64_t total_elements = 1;
        for (size_t i = 0; i < shape.size(); ++i) {
            shape_str += std::to_string(shape[i]) + (i == shape.size() - 1 ? "" : ", ");
            total_elements *= shape[i];
        }
        shape_str += "]";

        TF_DataType dtype_code = TF_TensorType(tensor.get_tensor().get());
        std::string dtype_str;
        if (dtype_code == TF_FLOAT) dtype_str = "float32";
        else if (dtype_code == TF_DOUBLE) dtype_str = "float64";
        else if (dtype_code == TF_INT32) dtype_str = "int32";
        else if (dtype_code == TF_INT64) dtype_str = "int64";
        else dtype_str = "type_" + std::to_string(dtype_code);

        ss << std::left << std::setw(40) << name << " | "
           << std::left << std::setw(10) << dtype_str << " | "
           << std::left << std::setw(15) << shape_str << "\n";
    }

    // --- Python Dictionary Export ---
    if (python_ready) {
        ss << "\n# Proc #"<<me <<" Python-ready dictionary (copy-paste for testing)\n";
        ss << "import numpy as np\n";
        ss << "grace_inputs = {\n";

        for (const auto& item : inputs) {
            std::string key = std::get<0>(item);
            const cppflow::tensor& tensor = std::get<1>(item);

            // Strip prefix "parallel_compute_" and ":0" suffix
            size_t prefix_pos = key.find("parallel_compute_");
            if (prefix_pos != std::string::npos) key.erase(prefix_pos, 17);
            size_t port_pos = key.find(":0");
            if (port_pos != std::string::npos) key.erase(port_pos, 2);

            TF_DataType dtype_code = TF_TensorType(tensor.get_tensor().get());
            void* data = TF_TensorData(tensor.get_tensor().get());

            // Map TF type to NumPy dtype string
            std::string np_dtype = "";
            if (dtype_code == TF_INT32) np_dtype = ", dtype=np.int32";
            else if (dtype_code == TF_INT64) np_dtype = ", dtype=np.int64";
            else if (dtype_code == TF_FLOAT) np_dtype = ", dtype=np.float32";
            else if (dtype_code == TF_DOUBLE) np_dtype = ", dtype=np.float64";

            ss << "    '" << key << "': np.array([";

            int64_t total_elements = 1;
            auto shape = tensor.shape().get_data<int64_t>();
            for (auto s : shape) total_elements *= s;

            // Export full content
            for (int64_t i = 0; i < total_elements; ++i) {
                if (dtype_code == TF_INT32) ss << static_cast<int32_t*>(data)[i];
                else if (dtype_code == TF_DOUBLE) ss << std::fixed << std::setprecision(14) << static_cast<double*>(data)[i];
                else if (dtype_code == TF_FLOAT) ss << std::fixed << std::setprecision(8) << static_cast<float*>(data)[i];
                else if (dtype_code == TF_INT64) ss << static_cast<int64_t*>(data)[i];
                if (i < total_elements - 1) ss << ", ";
            }

            // Construct Python shape tuple
            std::string py_shape = "(";
            for (size_t i = 0; i < shape.size(); ++i) {
                py_shape += std::to_string(shape[i]) + (i == shape.size() - 1 && shape.size() == 1 ? "," : "");
                if (i < shape.size() - 1) py_shape += ",";
            }
            py_shape += ")";

            if (key=="bond_vector")
                ss << "]" << np_dtype << ").reshape(" << py_shape << "),\n";
            else
                ss << "]" << np_dtype << "),\n";
        }
        ss << "}\n";
    }

    ss << std::string(100, '=') << "\n";
    LAMMPS_NS::utils::logmesg(lmp, "{}\n", ss.str());
}



void print_f_data(const double *f_data, int me, LAMMPS *lmp,
    const std::vector<int>& ind_i_vector, const std::vector<int>& ind_j_vector,
    int* tag) {
    int i;
    int j;
    //f-data is [nbonds,3]; ind_i_vector and ind_j_vector have corresponding i-j indices, print all this info pretty
    // ----------------------------------------------------------------------
    // DEBUG: Pretty Print Raw TensorFlow Pair Forces
    // ----------------------------------------------------------------------
    // We iterate up to n_real_neighbours to skip the padding/fake bonds
    // Ensure output from different processors doesn't get garbled
    // (Simple serializing via sleep; for strict ordering use MPI barriers)
    auto atom = lmp->atom;

    if (me == 0) {
        utils::logmesg(lmp, "\n[GRACE-DEBUG] Raw TF Output Tensor (f_data) Preview:\n");
        utils::logmesg(lmp, "Idx  | Proc | Atom I (Tag) -> Atom J (Tag) | Raw Force (fx, fy, fz)\n");
        utils::logmesg(lmp, "-----|------|------------------------------|--------------------------\n");
    }


    for (int k = 0; k < ind_i_vector.size(); ++k) {
        i = ind_i_vector[k];
        j = ind_j_vector[k];

        // Safety check: ensure indices are within local/ghost range
        //if (i < 0 || i >= nall || j < 0 || j >= nall) continue;

        tagint tag_i = tag[i];
        tagint tag_j = tag[j];

        double raw_fx = f_data[3*k + 0];
        double raw_fy = f_data[3*k + 1];
        double raw_fz = f_data[3*k + 2];

        // Optional: Filter out near-zero forces to reduce noise
        //if (raw_fx*raw_fx + raw_fy*raw_fy + raw_fz*raw_fz > 1e-20)
        {
            utils::logmesg(lmp, "{:<4} | {:<4} | {:>8} -> {:<8} | [{: .6f}, {: .6f}, {: .6f}]\n",
                           k, me, tag_i, tag_j, raw_fx, raw_fy, raw_fz);
        }
    }

    if (me == 0) utils::logmesg(lmp, "------------------------------------------------------------------\n");
    // ----------------------------------------------------------------------
}

} // namespace GRACE
#endif
