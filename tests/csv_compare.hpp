#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

inline std::vector<std::string> splitCsvLine(const std::string& line)
{
    std::vector<std::string> cells;
    std::stringstream ss(line);
    std::string cell;

    while (std::getline(ss, cell, ',')) {
        cells.push_back(cell);
    }

    return cells;
}

inline bool tryParseDouble(const std::string& s, double& value)
{
    try {
        size_t pos = 0;
        value = std::stod(s, &pos);
        return pos == s.size();
    } catch (...) {
        return false;
    }
}

inline bool nearRelative(double a, double b, double relTol, double absTol)
{
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= absTol + relTol * scale;
}

inline bool compareCsvFiles(
    const std::string& computedPath,
    const std::string& referencePath,
    double relTol = 1e-5,
    double absTol = 1e-12
)
{
    std::ifstream computed(computedPath);
    std::ifstream reference(referencePath);

    if (!computed.is_open()) {
        std::cerr << "Cannot open computed file: " << computedPath << "\n";
        return false;
    }

    if (!reference.is_open()) {
        std::cerr << "Cannot open reference file: " << referencePath << "\n";
        return false;
    }

    std::string computedLine;
    std::string referenceLine;
    std::size_t lineNumber = 0;

    while (true) {
        const bool hasComputed = static_cast<bool>(std::getline(computed, computedLine));
        const bool hasReference = static_cast<bool>(std::getline(reference, referenceLine));

        if (!hasComputed && !hasReference) {
            return true;
        }

        ++lineNumber;

        if (hasComputed != hasReference) {
            std::cerr << "Different number of lines in:\n"
                      << "computed:  " << computedPath << "\n"
                      << "reference: " << referencePath << "\n";
            return false;
        }

        const auto computedCells = splitCsvLine(computedLine);
        const auto referenceCells = splitCsvLine(referenceLine);

        if (computedCells.size() != referenceCells.size()) {
            std::cerr << "Different number of columns at line " << lineNumber << "\n"
                      << "computed:  " << computedPath << "\n"
                      << "reference: " << referencePath << "\n";
            return false;
        }

        for (std::size_t j = 0; j < computedCells.size(); ++j) {
            double a = 0.0;
            double b = 0.0;

            const bool aIsNumber = tryParseDouble(computedCells[j], a);
            const bool bIsNumber = tryParseDouble(referenceCells[j], b);

            if (aIsNumber && bIsNumber) {
                if (!nearRelative(a, b, relTol, absTol)) {
                    std::cerr << "CSV numerical mismatch\n"
                              << "file:      " << computedPath << "\n"
                              << "line:      " << lineNumber << "\n"
                              << "column:    " << j + 1 << "\n"
                              << "computed:  " << a << "\n"
                              << "reference: " << b << "\n";
                    return false;
                }
            } else {
                if (computedCells[j] != referenceCells[j]) {
                    std::cerr << "CSV text mismatch\n"
                              << "file:      " << computedPath << "\n"
                              << "line:      " << lineNumber << "\n"
                              << "column:    " << j + 1 << "\n"
                              << "computed:  " << computedCells[j] << "\n"
                              << "reference: " << referenceCells[j] << "\n";
                    return false;
                }
            }
        }
    }
}