#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <thread>
#include <mutex>
#include <atomic>
#include <stdexcept> // For exceptions
#include <limits>    // For numeric_limits
#include <cstdlib>   // For system()
#include <cstdio>    // For tmpnam (or use mkstemp)
#include <memory>    // For unique_ptr if needed
#include <filesystem> // For checking file existence (C++17)

// --- KMC API Headers ---
// Assuming headers are directly in 'include/kmc_api/' relative to this file
// or that the include path is set correctly during compilation (-I include)
#include "kmc_api/kmc_file.h" // Includes kmer_api.h, kmer_defs.h indirectly

// --- Configuration ---
// KMER_SIZE will be determined *from the KMC database* if one is provided,
// otherwise, it must be specified for the KMC run.
// Let's make it configurable via command line if running KMC.
int KMER_SIZE = -1; // Indicate K needs to be set or read
const unsigned int NUM_THREADS = std::min(128u, std::thread::hardware_concurrency()); // Use available cores, capped

// KMC Memory limit (e.g., in GB) - adjust based on your system
const int KMC_MEM_GB = 16;
// Minimum count threshold to consider when calculating median frequency for a read
const uint32_t MIN_KMER_COUNT_THRESHOLD_FOR_MEDIAN = 2;
// Minimum count threshold for KMC itself (optional, -ciX parameter)
// Only relevant if running KMC de novo
uint32_t KMC_MIN_COUNT_THRESHOLD = 2; // Changed from const to be configurable via command line

const size_t READ_BATCH_SIZE = 10000; // Process reads in batches for Phase 2

// --- Data Structures ---
struct FastqRecord {
    std::string header;
    std::string sequence;
    // std::string quality; // We don't need quality for this analysis
};

struct ReadStats {
    std::string read_id;
    double gc_content;
    uint64_t median_kmer_freq; // Use uint64_t as median can reflect large counts
};

// --- Global Variables ---
std::mutex results_mutex; // Mutex to protect shared results vector
std::atomic<uint64_t> reads_processed_count(0); // Atomic counter for progress (Phase 2)
std::atomic<uint64_t> total_reads_in_files(0); // Estimate total reads for progress

// --- Helper Functions ---

// Basic FASTQ Parser (simplified)
template <typename Stream>
bool read_fastq_record(Stream& input, FastqRecord& record) {
    std::string line1, line3, line4;
    // Skip empty lines at the beginning or between records
    while (input.peek() == '\n' || input.peek() == '\r') {
        input.ignore(1);
    }
    if (!std::getline(input, line1) ||
        !std::getline(input, record.sequence) ||
        !std::getline(input, line3) ||
        !std::getline(input, line4)) {
        return false; // Reached EOF or incomplete record at end
    }

    if (line1.empty() || line1[0] != '@' ||
        line3.empty() || line3[0] != '+' ||
        record.sequence.length() != line4.length()) {
         std::cerr << "\nWarning: Skipping potentially malformed FASTQ record near header: " << line1 << std::endl;
         // Attempt to find the next '@'
         while(input.good() && input.peek() != '@') {
             if (!input.ignore(std::numeric_limits<std::streamsize>::max(), '\n')) {
                 return false; // EOF reached while skipping
             }
         }
         if (!input.good()) return false;
         // We are positioned at the start of a potential next record
         return read_fastq_record(input, record); // Recursive call to retry
    }

    // Extract read ID before first space/tab
    size_t first_delim = line1.find_first_of(" \t");
    record.header = line1.substr(1, first_delim != std::string::npos ? first_delim - 1 : std::string::npos);
    return true;
}

// Simple function to estimate total reads (for progress)
uint64_t estimate_total_reads(const std::vector<std::string>& filenames) {
    uint64_t total_lines = 0;
    for (const auto& filename : filenames) {
        // TODO: Handle compressed files (e.g., pipe through zcat | wc -l)
        std::ifstream counter_stream(filename);
        if (!counter_stream) {
            std::cerr << "Warning: Could not open file " << filename << " for read estimation." << std::endl;
            continue;
        }
        total_lines += std::count(std::istreambuf_iterator<char>(counter_stream),
                                  std::istreambuf_iterator<char>(), '\n');
        counter_stream.close();
    }
    // Each FASTQ record is 4 lines
    return (total_lines + 2) / 4; // Add 2 for integer division robustness
}


double calculate_gc(const std::string& seq) {
    if (seq.empty()) return 0.0;
    size_t gc_count = 0;
    size_t valid_bases = 0;
    for (char base : seq) {
        char upper_base = std::toupper(base);
        if (upper_base == 'G' || upper_base == 'C') { gc_count++; valid_bases++; }
        else if (upper_base == 'A' || upper_base == 'T') { valid_bases++; }
        // Ignore Ns and other characters for GC calculation
    }
    return (valid_bases == 0) ? 0.0 : static_cast<double>(gc_count) / valid_bases;
}

// Calculates median using nth_element. Operates on uint64_t now.
uint64_t calculate_median(std::vector<uint64_t>& values) {
     if (values.empty()) return 0; // Return 0 if no k-mers met the threshold
     size_t mid = values.size() / 2;
     std::nth_element(values.begin(), values.begin() + mid, values.end());
     if (values.size() % 2 != 0) { // Odd number of elements
         return values[mid];
     } else { // Even number of elements - use lower of the two middle elements for simplicity
         // uint64_t mid_val1 = values[mid - 1]; // nth_element guarantees this is <= values[mid]
         // For large datasets, the exact median definition (avg of middle 2) might not matter much
         // Let's return values[mid] which is the upper of the two middle ones after nth_element
         // To get the average, we'd need another nth_element on the first half.
         // Simpler approach: Return the element at index `mid` (the upper middle value).
         // If you need the precise average:
         // std::nth_element(values.begin(), values.begin() + mid - 1, values.begin() + mid); // Find the lower middle in the first half
         // return (values[mid - 1] + values[mid]) / 2;
         return values[mid]; // Return upper-middle element for efficiency
     }
}


// --- Worker function for processing a batch of reads (Phase 2 using KMC API) ---
void process_reads_batch_kmc_api(
    const std::vector<FastqRecord>& batch,
    CKMCFile& kmc_reader, // Pass by reference
    int kmer_len_from_db,  // Pass k-mer length read from DB
    std::vector<ReadStats>& local_results,
    uint64_t total_reads_estimate) // For progress reporting
{
    local_results.reserve(batch.size());
    std::vector<uint32_t> kmer_counts_buffer; // Buffer for GetCountersForRead

    for (const auto& record : batch) {
        if (record.sequence.length() < static_cast<size_t>(kmer_len_from_db)) {
            // Update progress even for skipped reads
             reads_processed_count.fetch_add(1);
            continue; // Skip reads shorter than k
        }

        double gc = calculate_gc(record.sequence);
        std::vector<uint64_t> freqs_for_median; // Store frequencies meeting threshold

        // Use KMC API to get counts for all k-mers in the read
        if (kmc_reader.GetCountersForRead(record.sequence, kmer_counts_buffer)) {
            freqs_for_median.reserve(kmer_counts_buffer.size()); // Estimate size
            for (uint32_t count : kmer_counts_buffer) {
                // GetCountersForRead returns 0 for k-mers with 'N' or below DB min_count
                if (count >= MIN_KMER_COUNT_THRESHOLD_FOR_MEDIAN) {
                    freqs_for_median.push_back(static_cast<uint64_t>(count));
                }
            }
        } else {
             // This might happen if the read is shorter than k, but we checked above.
             // Could indicate other issues.
             std::cerr << "\nWarning: GetCountersForRead failed for read: " << record.header << std::endl;
        }

        uint64_t median_freq = calculate_median(freqs_for_median);
        local_results.push_back({record.header, gc, median_freq});

        // Update progress (atomic operation)
        uint64_t current_processed = reads_processed_count.fetch_add(1) + 1;
        if (total_reads_estimate > 0 && current_processed % 100000 == 0) {
             double progress = static_cast<double>(current_processed) / total_reads_estimate * 100.0;
             printf("\rProgress: %.2f%% (%llu / %llu reads)", progress, (long long unsigned)current_processed, (long long unsigned)total_reads_estimate);
            fflush(stdout);
        } else if (current_processed % 100000 == 0){
            printf("\rProgress: %llu reads", (long long unsigned)current_processed);
            fflush(stdout);
        }
    } // end record iteration
}


// --- Main Function ---
int main(int argc, char* argv[]) {
    std::string output_prefix;
    std::vector<std::string> input_files_str;
    std::string existing_kmc_db_prefix = "";
    bool kmc_run_needed = true;

    // --- Basic Argument Parsing ---
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--use-kmc-db" && i + 1 < argc) {
            existing_kmc_db_prefix = argv[++i];
            kmc_run_needed = false;
        } else if (arg == "-k" && i + 1 < argc) {
             try {
                KMER_SIZE = std::stoi(argv[++i]);
                if (KMER_SIZE <= 0) throw std::invalid_argument("K must be positive");
             } catch (const std::exception& e) {
                 std::cerr << "Error: Invalid k-mer size provided with -k: " << argv[i] << std::endl;
                 return 1;
             }
        } else if (arg == "--min-count" && i + 1 < argc) {
             try {
                KMC_MIN_COUNT_THRESHOLD = std::stoi(argv[++i]);
                if (KMC_MIN_COUNT_THRESHOLD < 1) throw std::invalid_argument("Min count must be positive");
             } catch (const std::exception& e) {
                 std::cerr << "Error: Invalid min count threshold provided with --min-count: " << argv[i] << std::endl;
                 return 1;
             }
        } else if (arg == "--help" || arg == "-h") {
             std::cerr << "Usage: " << argv[0] << " <output_prefix> [-k <kmer_size>] [--min-count <min_count>] [--use-kmc-db <existing_kmc_prefix>] <input_fastq1> [input_fastq2 ...]" << std::endl;
             std::cerr << "Options:" << std::endl;
             std::cerr << "  <output_prefix>         : Prefix for output files (_gc_abundance.tsv, _kmc_db*, etc.)" << std::endl;
             std::cerr << "  -k <kmer_size>          : K-mer size (required if KMC needs to run, ignored if --use-kmc-db is used)." << std::endl;
             std::cerr << "  --min-count <min_count> : Minimum kmer count threshold for KMC (default: 2)" << std::endl;
             std::cerr << "  --use-kmc-db <prefix> : Use existing KMC database files (<prefix>.kmc_pre, <prefix>.kmc_suf)." << std::endl;
             std::cerr << "                          If provided, the KMC counting step is skipped." << std::endl;
             std::cerr << "  <input_fastq...>       : One or more input FASTQ files (gzip not supported directly)." << std::endl;
             return 0;
        } else if (output_prefix.empty()) {
            output_prefix = arg;
        } else {
            input_files_str.push_back(arg);
        }
    }

    if (output_prefix.empty() || input_files_str.empty()) {
        std::cerr << "Error: Output prefix and at least one input FASTQ file are required." << std::endl;
        std::cerr << "Usage: " << argv[0] << " <output_prefix> [-k <kmer_size>] [--min-count <min_count>] [--use-kmc-db <existing_kmc_prefix>] <input_fastq1> [input_fastq2 ...]" << std::endl;
        return 1;
    }

    if (kmc_run_needed && KMER_SIZE <= 0) {
        std::cerr << "Error: K-mer size (-k) must be specified when running KMC de novo (i.e., when --use-kmc-db is not used)." << std::endl;
         return 1;
    }
     if (KMER_SIZE > 256 && kmc_run_needed) { // KMC limitation
         std::cerr << "Error: KMC maximum k-mer length is 256." << std::endl;
         return 1;
     }

    std::string kmc_db_prefix_to_use = kmc_run_needed ? (output_prefix + "_kmc_db") : existing_kmc_db_prefix;
    std::string kmc_tmp_dir = output_prefix + "_kmc_tmp"; // Needed if running KMC
    std::string kmc_input_list = output_prefix + "_kmc_input.list"; // Needed if running KMC

    CKMCFile kmc_reader; // KMC database reader object
    CKMCFileInfo kmc_info; // To store DB info

    try {
        std::cout << "=== Virtual Centrifuge Analysis (KMC API Optimized) ===" << std::endl;
        std::cout << "Parameters:" << std::endl;
        std::cout << "  Output Prefix: " << output_prefix << std::endl;
        std::cout << "  Input Files: ";
        for(const auto& f : input_files_str) std::cout << f << " ";
        std::cout << std::endl;
        std::cout << "  Threads: " << NUM_THREADS << std::endl;
        std::cout << "  Min K-mer Freq for Median Calc: " << MIN_KMER_COUNT_THRESHOLD_FOR_MEDIAN << std::endl;

        if (kmc_run_needed) {
            // --- Phase 1: K-mer Counting with KMC (External Call) ---
            std::cout << "\n--- Phase 1: Counting K-mers using KMC ---" << std::endl;
            std::cout << "  K-mer size (K): " << KMER_SIZE << std::endl;
            std::cout << "  KMC Memory: " << KMC_MEM_GB << " GB" << std::endl;
            std::cout << "  KMC Min Count Threshold (-ci): " << KMC_MIN_COUNT_THRESHOLD << std::endl;

            // Create input file list for KMC
            { // Scope to ensure file stream is closed
                std::ofstream input_list_stream(kmc_input_list);
                if (!input_list_stream) {
                    throw std::runtime_error("Failed to create KMC input file list: " + kmc_input_list);
                }
                for (const auto& file : input_files_str) {
                    input_list_stream << file << std::endl;
                }
            } // input_list_stream closed here

            // Construct KMC command
            std::string kmc_command = "kmc";
            kmc_command += " -k" + std::to_string(KMER_SIZE); // K-mer size
            kmc_command += " -t" + std::to_string(NUM_THREADS); // Threads
            kmc_command += " -m" + std::to_string(KMC_MEM_GB); // Memory limit (GB)
            kmc_command += " -ci" + std::to_string(KMC_MIN_COUNT_THRESHOLD); // Min count threshold
            // kmc_command += " -cs <val>"; // Optional: Max count value
            // kmc_command += " -f<a/q/m>"; // Optional: Specify input format if needed (auto usually works)
            // kmc_command += " -b"; // Important: Ensure canonical k-mers are stored (KMC default, but explicit)
            kmc_command += " @" + kmc_input_list; // Input files list marker
            kmc_command += " " + kmc_db_prefix_to_use; // Output database prefix
            kmc_command += " " + kmc_tmp_dir; // Temporary directory path

            std::cout << "Executing KMC: " << kmc_command << std::endl;
             // Create tmp directory if it doesn't exist
            std::filesystem::create_directory(kmc_tmp_dir);
            int kmc_status = system(kmc_command.c_str()); // Execute KMC

            // Basic check: does the output db exist?
            bool db_pre_exists = std::filesystem::exists(kmc_db_prefix_to_use + ".kmc_pre");
            bool db_suf_exists = std::filesystem::exists(kmc_db_prefix_to_use + ".kmc_suf");

            if (kmc_status != 0 || !db_pre_exists || !db_suf_exists) {
                 std::string err_msg = "KMC execution failed or did not produce output files. ";
                 err_msg += "KMC Status: " + std::to_string(kmc_status);
                 err_msg += ", .kmc_pre exists: " + std::string(db_pre_exists ? "yes" : "no");
                 err_msg += ", .kmc_suf exists: " + std::string(db_suf_exists ? "yes" : "no");
                 err_msg += ". Check KMC logs or stderr.";
                throw std::runtime_error(err_msg);
            }
            std::cout << "KMC counting complete." << std::endl;
             // Clean up intermediate files from this phase if desired
             // std::remove(kmc_input_list.c_str());
             // std::filesystem::remove_all(kmc_tmp_dir); // Requires C++17

        } else {
            std::cout << "\n--- Phase 1: Skipped (Using existing KMC database: " << kmc_db_prefix_to_use << ") ---" << std::endl;
             // Check if the provided files exist
             if (!std::filesystem::exists(kmc_db_prefix_to_use + ".kmc_pre") ||
                 !std::filesystem::exists(kmc_db_prefix_to_use + ".kmc_suf")) {
                 throw std::runtime_error("Error: Provided KMC database files not found: " +
                                          kmc_db_prefix_to_use + ".kmc_pre/.kmc_suf");
             }
        }

        // --- Open KMC Database and Get Info ---
        std::cout << "\n--- Opening KMC database for Random Access ---" << std::endl;
        if (!kmc_reader.OpenForRA(kmc_db_prefix_to_use)) {
             throw std::runtime_error("Failed to open KMC database for Random Access: " + kmc_db_prefix_to_use);
        }

        // Get database info
        if (!kmc_reader.Info(kmc_info)) {
             kmc_reader.Close(); // Attempt to close before throwing
             throw std::runtime_error("Failed to get info from KMC database: " + kmc_db_prefix_to_use);
        }

         // Set KMER_SIZE from the database info
         KMER_SIZE = kmc_info.kmer_length;

        std::cout << "KMC DB Info: k=" << kmc_info.kmer_length
                  << ", Total Kmers (in DB): " << kmc_info.total_kmers
                  << ", Min Count Threshold (in DB): " << kmc_info.min_count // This is the -ci value used
                  << ", Max Count: " << kmc_info.max_count
                  << ", Strands Mode: " << (kmc_info.both_strands ? "Canonical (Both Strands)" : "Single Strand")
                  << std::endl;

        // Crucial check for virtual centrifuge concept
        if (!kmc_info.both_strands) {
            std::cerr << "Warning: KMC database was created without canonical k-mer mode (-b option was likely off)." << std::endl;
            std::cerr << "         Median frequencies will be based on single-strand counts, which might not reflect abundance correctly." << std::endl;
            // Proceeding, but the user should be aware.
        }

        // Verify KMC k-mer length matches user expectation if KMC was run de novo
        // (If DB was provided, KMER_SIZE is now set from the DB)
        // No check needed here as KMER_SIZE is now aligned with the DB


        // --- Phase 2: Process Reads for GC and Median Frequency using KMC API ---
        std::cout << "\n--- Phase 2: Calculating GC and Median K-mer Frequency (using KMC API) ---" << std::endl;
        reads_processed_count = 0; // Reset counter

        // Estimate total reads for progress bar
        total_reads_in_files = estimate_total_reads(input_files_str);
        std::cout << "Estimating total reads... approx. " << total_reads_in_files << std::endl;


        std::vector<ReadStats> all_results;
        // Pre-allocate based on estimate if available, otherwise guess
        all_results.reserve(total_reads_in_files.load() > 0 ? total_reads_in_files.load() : 1000000);
        std::vector<std::thread> threads;
        std::vector<FastqRecord> current_batch;
        current_batch.reserve(READ_BATCH_SIZE);

        // Re-iterate through input files for Phase 2
        for (const auto& filename : input_files_str) {
            // TODO: Handle gzipped files if necessary (e.g., pipe through zcat or use a library like zlib)
            std::ifstream fastq_stream(filename);
            if (!fastq_stream) {
                std::cerr << "\nError: Cannot open input file for Phase 2: " << filename << ". Skipping." << std::endl;
                continue; // Skip this file
            }
            std::cout << "\nProcessing file for Phase 2: " << filename << std::endl;

            FastqRecord record;
            while (read_fastq_record(fastq_stream, record)) {
                current_batch.push_back(std::move(record)); // Move record into batch
                if (current_batch.size() >= READ_BATCH_SIZE) {
                    // Wait for a thread slot if necessary
                    while (threads.size() >= NUM_THREADS) {
                        // Look for a finished thread
                        bool joined = false;
                        for (size_t i = 0; i < threads.size(); ++i) {
                            // Check if joinable and try a timed join or check future status if using futures
                            // Simple approach: Join the first one if joinable
                            if (threads[i].joinable()) { // Basic check, might block if first isn't done
                                threads[i].join();
                                threads.erase(threads.begin() + i);
                                joined = true;
                                break; // Exit inner loop once one thread is joined
                            }
                        }
                         if (!joined) {
                             // If no threads were immediately joinable, wait a bit
                             std::this_thread::sleep_for(std::chrono::milliseconds(10));
                         }
                    }

                    // Launch new thread, passing reference to the single KMC reader
                    threads.emplace_back([batch = std::move(current_batch), &kmc_reader, k_len = KMER_SIZE, &all_results, total_estimate = total_reads_in_files.load()]() mutable {
                        std::vector<ReadStats> local_results;
                        // Process the batch using the KMC API
                        process_reads_batch_kmc_api(batch, kmc_reader, k_len, local_results, total_estimate);
                        // Lock and merge results into the main vector
                        std::lock_guard<std::mutex> lock(results_mutex);
                        all_results.insert(all_results.end(),
                                           std::make_move_iterator(local_results.begin()),
                                           std::make_move_iterator(local_results.end()));
                    });
                    // current_batch was moved, clear it and reserve for next batch
                    current_batch.clear(); // Ensure it's empty after move
                    current_batch.reserve(READ_BATCH_SIZE);
                }
            } // End while reading records
            fastq_stream.close();
             printf("\rProgress: 100.00%% (%llu / %llu reads) - File %s Done.\n", (long long unsigned)reads_processed_count.load(), (long long unsigned)total_reads_in_files, filename.c_str()); fflush(stdout);

        } // End file loop for Phase 2

        // Process any remaining reads in the last batch
        if (!current_batch.empty()) {
             // Wait for a thread slot if necessary (copy logic from above)
            while (threads.size() >= NUM_THREADS) {
                 bool joined = false;
                 for (size_t i = 0; i < threads.size(); ++i) {
                     if (threads[i].joinable()) {
                         threads[i].join();
                         threads.erase(threads.begin() + i);
                         joined = true;
                         break;
                     }
                 }
                  if (!joined) std::this_thread::sleep_for(std::chrono::milliseconds(10));
             }

            threads.emplace_back([batch = std::move(current_batch), &kmc_reader, k_len = KMER_SIZE, &all_results, total_estimate = total_reads_in_files.load()]() mutable {
                std::vector<ReadStats> local_results;
                process_reads_batch_kmc_api(batch, kmc_reader, k_len, local_results, total_estimate);
                std::lock_guard<std::mutex> lock(results_mutex);
                all_results.insert(all_results.end(),
                                   std::make_move_iterator(local_results.begin()),
                                   std::make_move_iterator(local_results.end()));
            });
        }

        // Wait for all remaining threads to complete
        std::cout << "\nWaiting for final threads to complete..." << std::endl;
        for (auto& t : threads) {
            if (t.joinable()) {
                t.join();
            }
        }
        printf("\rProgress: 100.00%% (%llu / %llu reads) - All Processing Done.        \n", (long long unsigned)reads_processed_count.load(), (long long unsigned)total_reads_in_files); fflush(stdout);
        std::cout << "Read processing complete. Total reads analyzed in Phase 2: " << reads_processed_count << std::endl;

        // Close the KMC database reader
        kmc_reader.Close();
        std::cout << "KMC database closed." << std::endl;


        // --- Phase 3: Output Results ---
        std::cout << "\n--- Phase 3: Writing Results ---" << std::endl;
        std::string output_tsv_file = output_prefix + "_gc_abundance.tsv";
        std::ofstream out_stream(output_tsv_file);
        if (!out_stream) {
             // No need to close KMC reader here, already closed
             throw std::runtime_error("Error: Cannot open output file: " + output_tsv_file);
        }

        // Write header
        out_stream << "ReadID\tGC_Content\tMedian_Kmer_Frequency\n";
        // Write data rows
        uint64_t rows_written = 0;
        for (const auto& stats : all_results) {
            out_stream << stats.read_id << "\t"
                       << stats.gc_content << "\t"
                       << stats.median_kmer_freq << "\n";
            rows_written++;
        }
        out_stream.close();
        std::cout << "Results for " << rows_written << " reads written to: " << output_tsv_file << std::endl;

        // --- Generate Python plotting script ---
        std::string plot_script_file = output_prefix + "_plot.py";
        std::ofstream py_script(plot_script_file);
        if (py_script) {
            std::cout << "\n--- Generating Python Plotting Script ---" << std::endl;
            py_script << "#!/usr/bin/env python3\n";
            py_script << "# Auto-generated plotting script for Virtual Centrifuge output\n";
            py_script << "import pandas as pd\n";
            py_script << "import matplotlib.pyplot as plt\n";
            py_script << "import numpy as np\n";
            py_script << "import sys\n\n";
            
            py_script << "# Input file - can be overridden with command line argument\n";
            py_script << "tsv_file = '" << output_tsv_file << "'\n";
            py_script << "output_prefix = '" << output_prefix << "'\n";
            py_script << "kmer_size = " << KMER_SIZE << "\n\n";
            
            py_script << "# Allow command line overrides\n";
            py_script << "if len(sys.argv) > 1:\n";
            py_script << "    tsv_file = sys.argv[1]\n";
            py_script << "if len(sys.argv) > 2:\n";
            py_script << "    output_prefix = sys.argv[2]\n\n";
            
            py_script << "print(f\"Reading data from: {tsv_file}\")\n";
            py_script << "# Load the data\n";
            py_script << "df = pd.read_csv(tsv_file, sep='\\t')\n";
            py_script << "print(f\"Loaded {len(df)} data points\")\n\n";
            
            py_script << "# Filter out zeros for log scale\n";
            py_script << "df_plot = df[df['Median_Kmer_Frequency'] > 0].copy()\n";
            py_script << "print(f\"Plotting {len(df_plot)} data points (excluding zero frequencies)\")\n\n";
            
            py_script << "# Create hexbin plot (better for large datasets)\n";
            py_script << "plt.figure(figsize=(10, 8))\n";
            py_script << "hb = plt.hexbin(df_plot['GC_Content'], df_plot['Median_Kmer_Frequency'],\n";
            py_script << "               gridsize=100, cmap='viridis', bins='log',\n";
            py_script << "               xscale='linear', yscale='log')\n";
            py_script << "plt.xlabel('GC Content (Density Analogue)')\n";
            py_script << "plt.ylabel('Median K-mer Frequency (Abundance Analogue - Log Scale)')\n";
            py_script << "plt.title(f'Virtual Ultracentrifugation Plot (k={kmer_size})')\n";
            py_script << "cb = plt.colorbar(hb, label='Log(Number of Reads in Bin)')\n";
            py_script << "plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)\n";
            py_script << "plt.ylim(bottom=max(0.1, df_plot['Median_Kmer_Frequency'].min()))\n";
            py_script << "plt.tight_layout()\n";
            py_script << "hexbin_file = f\"{output_prefix}_hexbin_plot.png\"\n";
            py_script << "plt.savefig(hexbin_file, dpi=150)\n";
            py_script << "print(f\"Hexbin plot saved to: {hexbin_file}\")\n\n";
            
            py_script << "# Create scatter plot (for smaller datasets)\n";
            py_script << "plt.figure(figsize=(10, 8))\n";
            py_script << "plt.scatter(df_plot['GC_Content'], df_plot['Median_Kmer_Frequency'],\n";
            py_script << "           alpha=0.1, s=1, c='blue')\n";
            py_script << "plt.xlabel('GC Content (Density Analogue)')\n";
            py_script << "plt.ylabel('Median K-mer Frequency (Abundance Analogue)')\n";
            py_script << "plt.title(f'Virtual Ultracentrifugation Plot - Scatter (k={kmer_size})')\n";
            py_script << "plt.yscale('log')  # Log scale for frequency\n";
            py_script << "plt.grid(True, alpha=0.3)\n";
            py_script << "plt.tight_layout()\n";
            py_script << "scatter_file = f\"{output_prefix}_scatter_plot.png\"\n";
            py_script << "plt.savefig(scatter_file, dpi=150)\n";
            py_script << "print(f\"Scatter plot saved to: {scatter_file}\")\n";
            py_script << "print(\"Plotting complete!\")\n";
            
            py_script.close();
            
            // Make the script executable
#ifndef _WIN32
            std::string chmod_cmd = "chmod +x " + plot_script_file;
            system(chmod_cmd.c_str());
#endif
            
            std::cout << "Python plotting script generated: " << plot_script_file << std::endl;
            std::cout << "To generate plots, run:" << std::endl;
            std::cout << "  python " << plot_script_file << std::endl;
            std::cout << "  or" << std::endl;
            std::cout << "  ./" << plot_script_file << std::endl;
        }

        // --- Plotting Guidance (same as before) ---
        std::cout << "\n--- Plotting ---" << std::endl;
        std::cout << "A plotting script has been generated for you. However, if you want to customize the plots," << std::endl;
        std::cout << "here's a guide to the data format and plotting approaches:" << std::endl;
        std::cout << "To visualize the 'virtual centrifuge' results, plot the data from:" << std::endl;
        std::cout << "  '" << output_tsv_file << "'" << std::endl;
        std::cout << "Use GC_Content on the X-axis and Median_Kmer_Frequency on the Y-axis." << std::endl;
        std::cout << "Recommended tools: Python (with Matplotlib/Seaborn), R (with ggplot2), or Gnuplot." << std::endl;
        std::cout << "Example (Python/Matplotlib):" << std::endl;
        std::cout << "  import pandas as pd" << std::endl;
        std::cout << "  import matplotlib.pyplot as plt" << std::endl;
        std::cout << "  import numpy as np # for log scale " << std::endl;
        std::cout << "  df = pd.read_csv('" << output_tsv_file << "', sep='\\t')" << std::endl;
        std::cout << "  # Filter out potential zeros before log scale if necessary" << std::endl;
        std::cout << "  df_plot = df[df['Median_Kmer_Frequency'] > 0].copy() # Use .copy() to avoid SettingWithCopyWarning" << std::endl;
        std::cout << "  # Optional: Add small value before log transform if zeros are meaningful but problematic for log" << std::endl;
        std::cout << "  # df_plot['Median_Kmer_Frequency'] = df_plot['Median_Kmer_Frequency'] + 0.1 " << std::endl;
        std::cout << "  plt.figure(figsize=(10, 8))" << std::endl;
        std::cout << "  # Use hexbin for large datasets to show density" << std::endl;
        std::cout << "  hb = plt.hexbin(df_plot['GC_Content'], df_plot['Median_Kmer_Frequency'], gridsize=100, cmap='viridis', bins='log', xscale='linear', yscale='log')" << std::endl;
        std::cout << "  plt.xlabel('GC Content (Density Analogue)')" << std::endl;
        std::cout << "  plt.ylabel('Median K-mer Frequency (Abundance Analogue - Log Scale)')" << std::endl;
        std::cout << "  plt.title('Virtual Ultracentrifugation Plot (k=" << KMER_SIZE << ")')" << std::endl;
        std::cout << "  cb = plt.colorbar(hb, label='Log(Number of Reads in Bin)')" << std::endl;
        std::cout << "  plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)" << std::endl;
        std::cout << "  plt.ylim(bottom=max(0.1, df_plot['Median_Kmer_Frequency'].min())) # Adjust y-axis bottom limit" << std::endl;
        std::cout << "  plt.tight_layout()" << std::endl;
        std::cout << "  plt.savefig('" << output_prefix << "_gc_abundance_plot.png', dpi=150)" << std::endl;
        std::cout << "  # plt.show()" << std::endl;


        std::cout << "\nAnalysis finished successfully." << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "\n\nAn error occurred: " << e.what() << std::endl;
        // Attempt to close KMC reader if it might be open
        if (kmc_reader.Close()) { // Close returns true if it was opened and now closed
             std::cerr << "KMC database reader closed after error." << std::endl;
        }
        // Consider adding cleanup for temporary files/dirs here if needed
        return 1;
    } catch (...) {
        std::cerr << "\n\nAn unknown error occurred." << std::endl;
         if (kmc_reader.Close()) {
             std::cerr << "KMC database reader closed after error." << std::endl;
         }
        // Consider adding cleanup for temporary files/dirs here if needed
        return 1;
    }

    // Optional: Clean up KMC temporary files and database if they were created here
    // Careful not to delete user-provided DBs!
    // if (kmc_run_needed) {
    //     std::remove(kmc_input_list.c_str());
    //     // std::remove((kmc_db_prefix_to_use + ".kmc_pre").c_str()); // Keep the DB by default
    //     // std::remove((kmc_db_prefix_to_use + ".kmc_suf").c_str()); // Keep the DB by default
    //     if (std::filesystem::exists(kmc_tmp_dir)) { // Check before removing
    //         std::filesystem::remove_all(kmc_tmp_dir); // C++17 way to remove directory
    //     }
    // }

    return 0;
}