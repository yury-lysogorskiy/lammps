#include <iostream>
#include <fstream>
#include <string>

class CSVWriter {
public:
    CSVWriter(const std::string& filename) : filename_(filename) {}

    void writeRow(const std::string& data) {
        std::ofstream file;
        // Use the append mode to add data to the end of the file if it exists
        file.open(filename_, std::ios::out | std::ios::app);

        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename_ << std::endl;
            return;
        }

        file << data << std::endl;
        file.close();
    }

private:
    std::string filename_;
};