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

void print_tf_inputs(const std::vector<std::tuple<std::string, cppflow::tensor>>& inputs,
                            int me, LAMMPS_NS::LAMMPS *lmp, bool python_ready) {
    std::stringstream ss;
    ss << "\n" << std::string(100, '=') << "\n";
    ss << "[GRACE-DEBUG, Proc #" << me << "]: TensorFlow Graph Inputs & Content\n";
    ss << std::string(100, '-') << "\n";

    ss << std::left << std::setw(40) << "Input Name (:port)"
       << " | " << std::setw(10) << "DataType"
       << " | " << std::setw(15) << "Shape"
       << " | " << "Content Preview (First 10)" << "\n";
    ss << std::string(100, '-') << "\n";

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
           << std::left << std::setw(15) << shape_str << " | ";

        void* data = TF_TensorData(tensor.get_tensor().get());
        int preview_count = std::min((int64_t)10, total_elements);

        for (int i = 0; i < preview_count; ++i) {
            if (dtype_code == TF_INT32) ss << static_cast<int32_t*>(data)[i];
            else if (dtype_code == TF_DOUBLE) ss << std::fixed << std::setprecision(6) << static_cast<double*>(data)[i];
            else if (dtype_code == TF_FLOAT) ss << std::fixed << std::setprecision(6) << static_cast<float*>(data)[i];
            else if (dtype_code == TF_INT64) ss << static_cast<int64_t*>(data)[i];
            if (i < preview_count - 1) ss << ", ";
        }
        if (total_elements > 10) ss << "...";
        ss << "\n";
    }

    // --- Python Dictionary Export ---
    if (python_ready) {
        ss << "\n# Python-ready dictionary (copy-paste for testing)\n";
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

            ss << "]" << np_dtype << ").reshape(" << py_shape << "),\n";
        }
        ss << "}\n";
    }

    ss << std::string(100, '=') << "\n";
    LAMMPS_NS::utils::logmesg(lmp, "{}\n", ss.str());
}
