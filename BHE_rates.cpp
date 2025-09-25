/*
 * BHE_rates.cpp
 *
 * Utility to estimate minimal encoding lengths (and resulting rates) for
 * bounded-homopolymer DNA encodings while reusing the reference encoder.
 */

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#define BOUNDED_HOMOPOLYMER_NO_MAIN
#include "BoundedHomopolymerEncoding.cpp"
#undef BOUNDED_HOMOPOLYMER_NO_MAIN

class CoutSilencer {
public:
    CoutSilencer() : old_buf_(std::cout.rdbuf(buffer_.rdbuf())) {}
    ~CoutSilencer() { std::cout.rdbuf(old_buf_); }

private:
    std::ostringstream buffer_;
    std::streambuf* old_buf_;
};

int capacity_for(int k, int encoding_length) {
    CoutSilencer silence;
#define BOUNDED_HOMOPOLYMER_SILENT
    BoundedHomopolymerEncoder encoder(k, encoding_length, encoding_length);
#undef BOUNDED_HOMOPOLYMER_SILENT
    return encoder.max_data_length();
}

int parse_positive_int(const std::string& value, const std::string& label) {
    try {
        int parsed = std::stoi(value);
        if (parsed <= 0) {
            throw std::invalid_argument("non-positive");
        }
        return parsed;
    } catch (const std::exception&) {
        throw std::invalid_argument(std::string("Invalid ") + label + ": " + value);
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <encoding_length>" << std::endl;
        return 1;
    }

    int encoding_length;
    try {
        encoding_length = parse_positive_int(argv[1], "encoding length");
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    std::cout << "Encoding length: " << encoding_length << std::endl << std::endl;

    std::cout << std::left << std::setw(6) << "k"
              << std::setw(18) << "max_input_bits"
              << std::setw(18) << "rate(bits/base)" << std::endl;
    std::cout << std::string(42, '-') << std::endl;

    for (int k = 1; k <= 5; ++k) {
        int capacity = capacity_for(k, encoding_length);
        std::ostringstream rate_stream;
        rate_stream << std::fixed << std::setprecision(6)
                    << static_cast<double>(capacity) / static_cast<double>(encoding_length);

        std::cout << std::setw(6) << k
                  << std::setw(18) << capacity
                  << std::setw(18) << rate_stream.str() << std::endl;
    }

    return 0;
}
