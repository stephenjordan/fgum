/**
 * QAOA Objective Verifier
 *
 * Reads JSON result files (results_k*_D*.json) and verifies each objective
 * (satisfaction fraction) by recomputing it via the Basso/Villalonga branch
 * tensor contraction (double-double precision by default).
 *
 * Companion to arXiv:2604.24633
 *
 * Usage:
 *   ./verify results_k3_D4.json
 *   ./verify results_k3_D4.json results_k2_D3.json --precision float64
 *   ./verify qaoa_data/ --p 1-8
 *
 * JSON format:
 *   {"k": 3, "D": 4, "convention": "e^{-ig Z^k}",
 *    "results": {"1": {"gammas": [...], "betas": [...], "objective": ...}, ...}}
 *
 * Angle convention: exp(-i * gamma * Z^k), mixer Rx(2*beta)
 */

#include "contract.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <dirent.h>
#include <sys/stat.h>

// ─────────────────────────────────────────────────────────────────────────────
// Minimal JSON parser (handles the specific results_k*_D*.json structure)
// ─────────────────────────────────────────────────────────────────────────────

static std::string read_file(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) return "";
    return std::string((std::istreambuf_iterator<char>(f)),
                        std::istreambuf_iterator<char>());
}

static size_t skip_ws(const std::string& s, size_t i) {
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t' || s[i] == '\n' || s[i] == '\r'))
        i++;
    return i;
}

static std::string parse_string(const std::string& s, size_t& i) {
    i = skip_ws(s, i);
    if (i >= s.size() || s[i] != '"') return "";
    i++;  // skip opening quote
    std::string result;
    while (i < s.size() && s[i] != '"') {
        if (s[i] == '\\') { i++; if (i < s.size()) result += s[i]; }
        else result += s[i];
        i++;
    }
    i++;  // skip closing quote
    return result;
}

static double parse_number(const std::string& s, size_t& i) {
    i = skip_ws(s, i);
    size_t start = i;
    if (i < s.size() && (s[i] == '-' || s[i] == '+')) i++;
    while (i < s.size() && (isdigit(s[i]) || s[i] == '.' || s[i] == 'e' || s[i] == 'E' || s[i] == '+' || s[i] == '-')) {
        if ((s[i] == '+' || s[i] == '-') && i > start && s[i-1] != 'e' && s[i-1] != 'E') break;
        i++;
    }
    return std::stod(s.substr(start, i - start));
}

static std::vector<double> parse_double_array(const std::string& s, size_t& i) {
    std::vector<double> arr;
    i = skip_ws(s, i);
    if (i >= s.size() || s[i] != '[') return arr;
    i++;  // skip [
    while (true) {
        i = skip_ws(s, i);
        if (i >= s.size() || s[i] == ']') { i++; break; }
        arr.push_back(parse_number(s, i));
        i = skip_ws(s, i);
        if (i < s.size() && s[i] == ',') i++;
    }
    return arr;
}

// Skip a JSON value (string, number, array, object, bool, null)
static void skip_value(const std::string& s, size_t& i) {
    i = skip_ws(s, i);
    if (i >= s.size()) return;
    if (s[i] == '"') { parse_string(s, i); return; }
    if (s[i] == '[') {
        i++; int depth = 1;
        while (i < s.size() && depth > 0) {
            if (s[i] == '[') depth++;
            else if (s[i] == ']') depth--;
            else if (s[i] == '"') { parse_string(s, i); continue; }
            i++;
        }
        return;
    }
    if (s[i] == '{') {
        i++; int depth = 1;
        while (i < s.size() && depth > 0) {
            if (s[i] == '{') depth++;
            else if (s[i] == '}') depth--;
            else if (s[i] == '"') { parse_string(s, i); continue; }
            i++;
        }
        return;
    }
    // number, bool, null
    while (i < s.size() && s[i] != ',' && s[i] != '}' && s[i] != ']' &&
           s[i] != ' ' && s[i] != '\n' && s[i] != '\r' && s[i] != '\t') i++;
}

// ─────────────────────────────────────────────────────────────────────────────
// Verification data structures and logic
// ─────────────────────────────────────────────────────────────────────────────

struct VerificationRow {
    int k, D, p;
    double objective_claimed;
    std::vector<double> gammas;
    std::vector<double> betas;
};

static std::vector<VerificationRow> parse_results_json(const std::string& filepath) {
    std::vector<VerificationRow> rows;
    std::string content = read_file(filepath);
    if (content.empty()) {
        fprintf(stderr, "Error: cannot read '%s'\n", filepath.c_str());
        return rows;
    }

    // Parse top-level k, D
    int k = 0, D = 0;
    size_t i = skip_ws(content, 0);
    if (i >= content.size() || content[i] != '{') return rows;
    i++;  // skip {

    while (i < content.size() && content[i] != '}') {
        i = skip_ws(content, i);
        std::string key = parse_string(content, i);
        i = skip_ws(content, i);
        if (i < content.size() && content[i] == ':') i++;

        if (key == "k") { k = (int)parse_number(content, i); }
        else if (key == "D") { D = (int)parse_number(content, i); }
        else if (key == "results") {
            // Parse results object
            i = skip_ws(content, i);
            if (i >= content.size() || content[i] != '{') break;
            i++;  // skip {

            while (i < content.size() && content[i] != '}') {
                i = skip_ws(content, i);
                std::string p_str = parse_string(content, i);
                int p = std::stoi(p_str);
                i = skip_ws(content, i);
                if (i < content.size() && content[i] == ':') i++;
                i = skip_ws(content, i);

                // Parse result entry object
                if (i >= content.size() || content[i] != '{') break;
                i++;  // skip {

                VerificationRow row;
                row.k = k; row.D = D; row.p = p;
                row.objective_claimed = 0.0;

                while (i < content.size() && content[i] != '}') {
                    i = skip_ws(content, i);
                    std::string field = parse_string(content, i);
                    i = skip_ws(content, i);
                    if (i < content.size() && content[i] == ':') i++;
                    i = skip_ws(content, i);

                    if (field == "gammas") { row.gammas = parse_double_array(content, i); }
                    else if (field == "betas") { row.betas = parse_double_array(content, i); }
                    else if (field == "objective") { row.objective_claimed = parse_number(content, i); }
                    else { skip_value(content, i); }

                    i = skip_ws(content, i);
                    if (i < content.size() && content[i] == ',') i++;
                }
                if (i < content.size()) i++;  // skip }

                if ((int)row.gammas.size() == p && (int)row.betas.size() == p)
                    rows.push_back(row);

                i = skip_ws(content, i);
                if (i < content.size() && content[i] == ',') i++;
            }
            if (i < content.size()) i++;  // skip }
        } else {
            skip_value(content, i);
        }

        i = skip_ws(content, i);
        if (i < content.size() && content[i] == ',') i++;
    }

    return rows;
}

static std::vector<std::string> find_json_files(const std::string& path) {
    std::vector<std::string> files;
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return files;

    if (S_ISDIR(st.st_mode)) {
        DIR* dir = opendir(path.c_str());
        if (!dir) return files;
        struct dirent* entry;
        while ((entry = readdir(dir)) != nullptr) {
            std::string name = entry->d_name;
            if (name.size() > 5 && name.substr(name.size() - 5) == ".json" &&
                name.find("results_k") == 0) {
                files.push_back(path + "/" + name);
            }
        }
        closedir(dir);
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}

static double compute_objective(const std::vector<double>& gammas,
                                const std::vector<double>& betas,
                                int p, int D, int k, bool use_dd) {
    double expectation = contract_symmetric_tree(gammas.data(), betas.data(),
                                                  p, D, k, use_dd, false);
    return (1.0 - expectation) / 2.0;
}

static const char* problem_label(int k, int D) {
    static char buf[48];
    if (k == 2)
        snprintf(buf, sizeof(buf), "MaxCut on %d-regular", D);
    else
        snprintf(buf, sizeof(buf), "MAX-%d-XORSAT on %d-regular", k, D);
    return buf;
}

static void print_separator(int width) {
    for (int i = 0; i < width; i++) putchar('-');
    putchar('\n');
}

static int run_verify(const std::vector<std::string>& paths, bool use_dd,
                      int p_min, int p_max) {
    printf("========================================================================\n");
    printf("  QAOA Objective Verifier\n");
    printf("  Companion to arXiv:2604.24633\n");
    printf("========================================================================\n\n");
    printf("Precision: %s\n", use_dd ? "double-double" : "float64");

    // Collect all rows from all files
    std::vector<VerificationRow> all_rows;
    for (const auto& path : paths) {
        auto files = find_json_files(path);
        for (const auto& file : files) {
            auto rows = parse_results_json(file);
            all_rows.insert(all_rows.end(), rows.begin(), rows.end());
        }
    }

    if (all_rows.empty()) {
        printf("No valid configurations found.\n");
        return 1;
    }

    // Filter by p range
    std::vector<VerificationRow> filtered;
    for (const auto& row : all_rows) {
        if (row.p >= p_min && row.p <= p_max)
            filtered.push_back(row);
    }

    if (p_min > 0 || p_max < 9999)
        printf("Filtering to p in [%d, %d]: %zu of %zu configurations.\n\n",
               p_min, p_max, filtered.size(), all_rows.size());
    else
        printf("Found %zu configurations to verify.\n\n", filtered.size());

    print_separator(100);

    int n_ok = 0, n_mismatch = 0, n_fail = 0;

    for (const auto& row : filtered) {
        auto t0 = std::chrono::high_resolution_clock::now();
        double obj = compute_objective(row.gammas, row.betas, row.p, row.D, row.k, use_dd);
        auto t1 = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(t1 - t0).count();

        if (std::isnan(obj) || !std::isfinite(obj)) {
            printf("  (k=%d, D=%d, p=%2d) [%-26s]  obj=NaN/Inf                       FAIL     [%.1fs]\n",
                   row.k, row.D, row.p, problem_label(row.k, row.D), elapsed);
            n_fail++;
            continue;
        }

        double diff = std::abs(obj - row.objective_claimed);
        const char* status = diff < 1e-6 ? "OK" : "MISMATCH";

        printf("  (k=%d, D=%d, p=%2d) [%-26s]  obj=%.12f  claimed=%.12f  delta=%.2e  %s  [%.1fs]\n",
               row.k, row.D, row.p, problem_label(row.k, row.D),
               obj, row.objective_claimed, diff, status, elapsed);

        if (diff < 1e-6) n_ok++;
        else n_mismatch++;
    }

    print_separator(100);
    printf("\nResults: %d OK, %d MISMATCH, %d FAIL (out of %zu total)\n",
           n_ok, n_mismatch, n_fail, filtered.size());

    if (n_fail == 0 && n_mismatch == 0) {
        printf("\n  All %zu evaluations verified successfully.\n\n", filtered.size());
        return 0;
    }
    return 1;
}

static void usage(const char* prog) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  %s <file.json|directory> [...] [--p P|P_MIN-P_MAX] [--precision dd|float64]\n\n", prog);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  --p 5         Only verify depth p=5\n");
    fprintf(stderr, "  --p 3-8       Only verify depths p=3 through p=8\n");
    fprintf(stderr, "  --precision   dd (default) or float64\n\n");
    fprintf(stderr, "Angle convention: exp(-i*gamma*Z^k), mixer Rx(2*beta)\n\n");
    fprintf(stderr, "Input: results_k*_D*.json files or a directory containing them.\n");
    fprintf(stderr, "JSON format:\n");
    fprintf(stderr, "  {\"k\": K, \"D\": D, \"convention\": \"e^{-ig Z^k}\",\n");
    fprintf(stderr, "   \"results\": {\"P\": {\"gammas\": [...], \"betas\": [...], \"objective\": ...}, ...}}\n");
}

int main(int argc, char** argv) {
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    std::vector<std::string> paths;
    int p_min = 0, p_max = 9999;
    bool use_dd = true;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--p") == 0 && i + 1 < argc) {
            std::string parg = argv[++i];
            auto dash = parg.find('-');
            if (dash != std::string::npos) {
                p_min = std::stoi(parg.substr(0, dash));
                p_max = std::stoi(parg.substr(dash + 1));
            } else {
                p_min = p_max = std::stoi(parg);
            }
        } else if (strcmp(argv[i], "--precision") == 0 && i + 1 < argc) {
            std::string prec = argv[++i];
            use_dd = (prec == "dd");
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
            return 0;
        } else if (argv[i][0] != '-') {
            paths.push_back(argv[i]);
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        }
    }

    if (paths.empty()) {
        // Default to qaoa_data/ relative to the binary location
        std::string bin_path = argv[0];
        auto slash = bin_path.rfind('/');
        std::string bin_dir = (slash != std::string::npos) ? bin_path.substr(0, slash) : ".";
        // binary is in cpp/build/, qaoa_data/ is at ../../qaoa_data/
        std::string default_dir = bin_dir + "/../../qaoa_data";
        struct stat st;
        if (stat(default_dir.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
            paths.push_back(default_dir);
        } else {
            // Try relative to cwd
            if (stat("qaoa_data", &st) == 0 && S_ISDIR(st.st_mode)) {
                paths.push_back("qaoa_data");
            } else {
                fprintf(stderr, "Error: no input files specified and qaoa_data/ not found.\n");
                usage(argv[0]);
                return 1;
            }
        }
    }

    return run_verify(paths, use_dd, p_min, p_max);
}
