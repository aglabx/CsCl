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
#include <unordered_map> // For storing k-mer counts
#include <memory>    // For unique_ptr if needed

// --- KMC API Headers ---
// Assuming headers are in a directory structure like "include/KMC/kmc_api/"
// Adjust include paths if KMC API is installed elsewhere
#include "KMC/kmc_api/kmc_file.h" // Includes kmer_api.h, kmer_defs.h indirectly
#include "KMC/kmc_api/kmer_api.h" // Include explicitly just in case

// --- Configuration ---
const int KMER_SIZE = 21; // Recommended k-mer size
const int NUM_THREADS = 100; // std::min(128, std::thread::hardware_concurrency()); // Use available cores



// KMC Memory limit (e.g., in GB) - adjust based on your system
const int KMC_MEM_GB = 16;
// Minimum count threshold to consider when calculating median frequency for a read
const uint32_t MIN_KMER_COUNT_THRESHOLD_FOR_MEDIAN = 2;
// Minimum count threshold for KMC itself (optional, -ciX parameter)
const uint32_t KMC_MIN_COUNT_THRESHOLD = 2; // Example: Filter k-mers counted only once during KMC run

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
    uint64_t median_kmer_freq;
};

// --- Global Variables ---
std::mutex results_mutex; // Mutex to protect shared results vector
std::atomic<uint64_t> reads_processed_count(0); // Atomic counter for progress (Phase 2)

// --- Helper Functions ---

// Basic FASTQ Parser (simplified) - Needed for Phase 2
template <typename Stream>
bool read_fastq_record(Stream& input, FastqRecord& record) {
    std::string line1, line3, line4;
    if (!std::getline(input, record.header) ||
        !std::getline(input, record.sequence) ||
        !std::getline(input, line3) ||
        !std::getline(input, line4)) {
        return false;
    }
    if (record.header.empty() || record.header[0] != '@' ||
        line3.empty() || line3[0] != '+' ||
        record.sequence.length() != line4.length()) {
         std::cerr << "\nWarning: Skipping potentially malformed FASTQ record near header: " << record.header << std::endl;
         while(input.good() && input.peek() != '@') {
             input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
         }
         if (!input.good()) return false;
         // Attempt to re-read starting from the '@' found
         return read_fastq_record(input, record);
    }
    // Extract read ID before first space
    size_t first_space = record.header.find(' ');
    record.header = record.header.substr(1, first_space != std::string::npos ? first_space - 1 : std::string::npos);
    return true;
}


double calculate_gc(const std::string& seq) {
    if (seq.empty()) return 0.0;
    size_t gc_count = 0;
    size_t valid_bases = 0;
    for (char base : seq) {
        char upper_base = std::toupper(base);
        if (upper_base == 'G' || upper_base == 'C') { gc_count++; valid_bases++; }
        else if (upper_base == 'A' || upper_base == 'T') { valid_bases++; }
    }
    return (valid_bases == 0) ? 0.0 : static_cast<double>(gc_count) / valid_bases;
}

uint64_t calculate_median(std::vector<uint64_t>& values) {
     if (values.empty()) return 0; // Return 0 if no k-mers met the threshold
     size_t mid = values.size() / 2;
     if (values.size() % 2 != 0) { // Odd number of elements
         std::nth_element(values.begin(), values.begin() + mid, values.end());
         return values[mid];
     } else { // Even number of elements - average middle two
         std::nth_element(values.begin(), values.begin() + mid - 1, values.end());
         uint64_t mid_val1 = values[mid - 1];
         // Need to find the element at 'mid' position as well
         std::nth_element(values.begin() + mid, values.begin() + mid, values.end());
         uint64_t mid_val2 = values[mid];
         // Return the average (integer division is fine here)
         return (mid_val1 + mid_val2) / 2;
     }
}

// --- K-mer Specific Helpers (Manual implementation needed for Phase 2) ---

// Simple DNA character to 2-bit code (A=00, C=01, G=10, T=11)
// Returns -1 for invalid characters. Uses CKmerAPI::num_codes initialization.
inline int dna_to_code(char base) {
    return CKmerAPI::num_codes[static_cast<unsigned char>(base)];
}

// 2-bit code to complement code
inline int complement_code(int code) {
    return code >= 0 ? 3 - code : -1; // A<->T (0<->3), C<->G (1<->2)
}

// Calculate reverse complement of a k-mer stored in a uint64_t
// Assumes k <= 32
uint64_t reverse_complement_kmer(uint64_t kmer, int k) {
    uint64_t rc_kmer = 0;
    uint64_t mask = 3ULL; // 0b11
    for (int i = 0; i < k; ++i) {
        int base_code = (kmer >> (2 * i)) & mask;
        int complement = complement_code(base_code);
        if (complement < 0) {
             // This should ideally not happen if input kmer is valid DNA
             // Return original or handle error appropriately
             // For simplicity, let's return original, assuming Ns were filtered
             return kmer;
        }
        // Place complemented base at the reversed position
        rc_kmer |= static_cast<uint64_t>(complement) << (2 * (k - 1 - i));
    }
    return rc_kmer;
}

// Get canonical k-mer (lexicographically smaller of k-mer and its reverse complement)
// Assumes k <= 32
uint64_t get_canonical_kmer(uint64_t kmer_val, int k) {
    uint64_t rc_kmer = reverse_complement_kmer(kmer_val, k);
    return std::min(kmer_val, rc_kmer);
}

// --- Worker function for processing a batch of reads (Phase 2) ---
void process_reads_batch(
    const std::vector<FastqRecord>& batch,
    const std::unordered_map<uint64_t, uint32_t>& kmer_counts, // Pass the loaded counts map
    std::vector<ReadStats>& local_results)
{
    local_results.reserve(batch.size());
    // Calculate mask for k-mers (works correctly for k <= 32)
    const uint64_t kmer_mask = (KMER_SIZE == 32) ? UINT64_MAX : (1ULL << (KMER_SIZE * 2)) - 1;

    for (const auto& record : batch) {
        if (record.sequence.length() < KMER_SIZE) {
            continue; // Skip reads shorter than k
        }

        double gc = calculate_gc(record.sequence);
        std::vector<uint64_t> freqs; // Store frequencies of k-mers in this read
        freqs.reserve(record.sequence.length()); // Rough estimate

        uint64_t current_kmer = 0;
        int valid_bases_in_kmer = 0;

        for (char base : record.sequence) {
            int code = dna_to_code(base); // Use helper using CKmerAPI::num_codes
            if (code < 0) { // Handle 'N' or invalid base
                valid_bases_in_kmer = 0; // Reset count of valid bases in the window
                current_kmer = 0;        // Reset k-mer value (optional, helps debugging)
                continue;                // Skip this base and potentially break the k-mer chain
            }

            // Shift previous k-mer left, add new base code, apply mask
            current_kmer = ((current_kmer << 2) | code) & kmer_mask;

            if (++valid_bases_in_kmer >= KMER_SIZE) {
                // We have a full k-mer of valid bases
                uint64_t canonical_kmer = get_canonical_kmer(current_kmer, KMER_SIZE);

                // Look up the canonical k-mer in the map
                auto it = kmer_counts.find(canonical_kmer);
                if (it != kmer_counts.end()) { // Found k-mer
                    // Check if count meets the threshold for median calculation
                    if (it->second >= MIN_KMER_COUNT_THRESHOLD_FOR_MEDIAN) {
                        freqs.push_back(it->second);
                    }
                }
                // If k-mer not found, it might have been filtered by KMC's -ci threshold,
                // or it's an error k-mer not present above threshold. Don't add to freqs.
            }
        } // end base iteration

        uint64_t median_freq = calculate_median(freqs);
        local_results.push_back({record.header, gc, median_freq});

        // Update progress (atomic operation)
        uint64_t current_processed = reads_processed_count.fetch_add(1) + 1;
        if (current_processed % 100000 == 0) {
            std::cout << "." << std::flush;
        }
    } // end record iteration
}


// --- Main Function ---
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <output_prefix> <input_fastq1> [input_fastq2 ...]" << std::endl;
        return 1;
    }

    std::string output_prefix = argv[1];
    std::vector<std::string> input_files_str; // Keep original strings
    for (int i = 2; i < argc; ++i) {
        input_files_str.push_back(argv[i]);
    }

    if (KMER_SIZE > 32) {
        std::cerr << "Error: This implementation currently only supports K <= 32." << std::endl;
        return 1;
    }


    std::cout << "=== Virtual Centrifuge Analysis (KMC 3.2.4 API version) ===" << std::endl;
    std::cout << "Parameters:" << std::endl;
    std::cout << "  K-mer size (K): " << KMER_SIZE << std::endl;
    std::cout << "  Threads: " << NUM_THREADS << std::endl;
    std::cout << "  KMC Memory: " << KMC_MEM_GB << " GB" << std::endl;
    std::cout << "  KMC Min Count Threshold (-ci): " << KMC_MIN_COUNT_THRESHOLD << std::endl;
    std::cout << "  Min K-mer Freq for Median Calc: " << MIN_KMER_COUNT_THRESHOLD_FOR_MEDIAN << std::endl;
    std::cout << "  Output Prefix: " << output_prefix << std::endl;
    std::cout << "  Input Files: ";
    for(const auto& f : input_files_str) std::cout << f << " ";
    std::cout << std::endl;

    // Generate temporary filenames for KMC
    std::string kmc_db_name = output_prefix + "_kmc_db"; // KMC uses this as prefix
    std::string kmc_tmp_dir = output_prefix + "_kmc_tmp"; // KMC needs a directory path
    std::string kmc_input_list = output_prefix + "_kmc_input.list";

    // Map to hold k-mer counts loaded from KMC DB
    // Key: uint64_t representation of canonical k-mer
    // Value: uint32_t count (KMC API typically uses uint32 for counts unless very high)
    std::unordered_map<uint64_t, uint32_t> kmer_counts_map;

    try {
        // --- Phase 1: K-mer Counting with KMC (External Call) ---
        std::cout << "\n--- Phase 1: Counting K-mers using KMC ---" << std::endl;

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
        // Basic command structure
        std::string kmc_command = "kmc";
        kmc_command += " -k" + std::to_string(KMER_SIZE); // K-mer size
        kmc_command += " -t" + std::to_string(NUM_THREADS); // Threads
        kmc_command += " -m" + std::to_string(KMC_MEM_GB); // Memory limit (GB)
        kmc_command += " -ci" + std::to_string(KMC_MIN_COUNT_THRESHOLD); // Min count threshold
        // kmc_command += " -cs1000000"; // Optional: Max count value if needed
        // kmc_command += " -f" // Optional: Specify input format (fastq/fasta/etc) if needed
        kmc_command += " @" + kmc_input_list; // Input files list marker
        kmc_command += " " + kmc_db_name; // Output database prefix
        kmc_command += " " + kmc_tmp_dir; // Temporary directory path

        std::cout << "Executing KMC: " << kmc_command << std::endl;
        int kmc_status = system(kmc_command.c_str()); // Execute KMC
        if (kmc_status != 0) {
             // system() return value interpretation can be complex,
             // but non-zero usually indicates an error.
            throw std::runtime_error("KMC execution failed. Check KMC logs or stderr. Status: " + std::to_string(kmc_status));
        }
        std::cout << "KMC counting complete." << std::endl;

        // --- Load KMC Counts into memory using KMC API ---
        std::cout << "\n--- Loading KMC counts into memory ---" << std::endl;
        CKMCFile kmc_db; // KMC database reader object
        if (!kmc_db.OpenForListing(kmc_db_name)) {
             throw std::runtime_error("Failed to open KMC database for listing: " + kmc_db_name);
        }

        // Get database info
        CKMCFileInfo kmc_info;
        if (!kmc_db.Info(kmc_info)) {
             kmc_db.Close(); // Attempt to close before throwing
             throw std::runtime_error("Failed to get info from KMC database: " + kmc_db_name);
        }

        std::cout << "KMC DB Info: k=" << kmc_info.kmer_length
                  << ", Total Kmers (in DB): " << kmc_info.total_kmers
                  << ", Min Count: " << kmc_info.min_count
                  << ", Max Count: " << kmc_info.max_count
                  << ", Strands: " << (kmc_info.both_strands ? "Both" : "Single")
                  << std::endl;

        // Verify KMC k-mer length matches configured KMER_SIZE
        if (kmc_info.kmer_length != KMER_SIZE) {
             std::cerr << "Warning: KMC database k (" << kmc_info.kmer_length
                       << ") differs from requested K (" << KMER_SIZE << ")" << std::endl;
             // Decide whether to proceed or exit based on severity
             // For now, we'll proceed but results might be unexpected
        }

        // Allocate CKmerAPI object once with the correct k-mer length from the DB
        CKmerAPI kmer_obj(kmc_info.kmer_length);
        uint64 count64; // ReadNextKmer expects uint64 for count
        std::vector<uint64> kmer_int_buffer; // Buffer for to_long

        uint64_t loaded_count = 0;
        uint64_t report_interval = 1000000;
        std::cout << "Loading k-mers: " << std::flush;

        // Iterate through all k-mers in the KMC database
        while (kmc_db.ReadNextKmer(kmer_obj, count64)) {
            // KMC should store canonical k-mers by default if both_strands is true
            // Convert the k-mer object to its uint64 representation
            kmer_obj.to_long(kmer_int_buffer);

            // Assuming K <= 32, the vector should have size 1
            if (!kmer_int_buffer.empty()) {
                 // Check if count exceeds uint32 max if needed, otherwise cast
                 if (count64 > UINT32_MAX) {
                     static bool count_warn = false;
                     if (!count_warn) {
                          std::cerr << "\nWarning: K-mer count " << count64 << " exceeds uint32_t max. Storing as UINT32_MAX." << std::endl;
                          count_warn = true;
                     }
                     kmer_counts_map[kmer_int_buffer[0]] = UINT32_MAX;
                 } else {
                      kmer_counts_map[kmer_int_buffer[0]] = static_cast<uint32_t>(count64);
                 }

                 loaded_count++;
                 if(loaded_count % report_interval == 0) std::cout << "." << std::flush;

            } else {
                 // This case should ideally not happen for K <= 32
                 std::cerr << "\nWarning: kmer_obj.to_long() returned empty vector for k=" << kmc_info.kmer_length << std::endl;
            }
        }
        kmc_db.Close(); // Close the KMC database file handles
        std::cout << "\nLoaded " << kmer_counts_map.size() << " unique k-mers (above KMC threshold " << KMC_MIN_COUNT_THRESHOLD << ") from KMC database." << std::endl;


        // --- Phase 2: Process Reads for GC and Median Frequency ---
        std::cout << "\n--- Phase 2: Calculating GC and Median K-mer Frequency ---" << std::endl;
        reads_processed_count = 0; // Reset counter

        std::vector<ReadStats> all_results;
        all_results.reserve(1000000); // Pre-allocate some space
        std::vector<std::thread> threads;
        std::vector<FastqRecord> current_batch;
        current_batch.reserve(READ_BATCH_SIZE);

        // Re-iterate through input files for Phase 2 using manual parsing
        for (const auto& filename : input_files_str) {
            // TODO: Handle gzipped files if necessary (e.g., pipe through zcat or use a library like zlib)
            std::ifstream fastq_stream(filename);
            if (!fastq_stream) {
                std::cerr << "\nError: Cannot open input file for Phase 2: " << filename << std::endl;
                continue; // Skip this file
            }
            std::cout << "\nProcessing file for Phase 2: " << filename << std::endl;

            FastqRecord record;
            while (read_fastq_record(fastq_stream, record)) {
                current_batch.push_back(record);
                if (current_batch.size() >= READ_BATCH_SIZE) {
                    if (threads.size() >= NUM_THREADS) {
                        // Wait for the oldest thread to finish
                        if (threads[0].joinable()) threads[0].join();
                        threads.erase(threads.begin()); // Remove finished thread
                    }
                    // Launch new thread, passing the const reference to the map
                    threads.emplace_back([batch = std::move(current_batch), &kmer_counts_map, &all_results]() mutable {
                        std::vector<ReadStats> local_results;
                        // Pass the map by const reference to the processing function
                        process_reads_batch(batch, kmer_counts_map, local_results);
                        // Lock and merge results into the main vector
                        std::lock_guard<std::mutex> lock(results_mutex);
                        all_results.insert(all_results.end(),
                                           std::make_move_iterator(local_results.begin()),
                                           std::make_move_iterator(local_results.end()));
                    });
                    // current_batch was moved, clear it (though move leaves it in valid state)
                    current_batch.clear();
                    current_batch.reserve(READ_BATCH_SIZE); // Reserve for next batch
                }
            }
            fastq_stream.close(); // Close the file stream for this file
        } // End file loop for Phase 2

        // Process any remaining reads in the last batch
        if (!current_batch.empty()) {
            if (threads.size() >= NUM_THREADS) {
                if (threads[0].joinable()) threads[0].join();
                threads.erase(threads.begin());
            }
            threads.emplace_back([batch = std::move(current_batch), &kmer_counts_map, &all_results]() mutable {
                std::vector<ReadStats> local_results;
                process_reads_batch(batch, kmer_counts_map, local_results);
                std::lock_guard<std::mutex> lock(results_mutex);
                all_results.insert(all_results.end(),
                                   std::make_move_iterator(local_results.begin()),
                                   std::make_move_iterator(local_results.end()));
            });
        }

        // Wait for all remaining threads to complete
        for (auto& t : threads) {
            if (t.joinable()) {
                t.join();
            }
        }
        std::cout << "\nRead processing complete. Total reads analyzed in Phase 2: " << reads_processed_count << std::endl;

        // --- Phase 3: Output Results ---
        std::cout << "\n--- Phase 3: Writing Results ---" << std::endl;
        std::string output_tsv_file = output_prefix + "_gc_abundance.tsv";
        std::ofstream out_stream(output_tsv_file);
        if (!out_stream) {
             throw std::runtime_error("Error: Cannot open output file: " + output_tsv_file);
        }

        // Write header
        out_stream << "ReadID\tGC_Content\tMedian_Kmer_Frequency\n";
        // Write data rows
        for (const auto& stats : all_results) {
            out_stream << stats.read_id << "\t"
                       << stats.gc_content << "\t"
                       << stats.median_kmer_freq << "\n";
        }
        out_stream.close(); // Close the output file stream
        std::cout << "Results written to: " << output_tsv_file << std::endl;

        // --- Plotting Guidance (same as before) ---
        std::cout << "\n--- Plotting ---" << std::endl;
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
        std::cout << "  df_plot = df[df['Median_Kmer_Frequency'] > 0]" << std::endl;
        std::cout << "  plt.figure(figsize=(10, 8))" << std::endl;
        std::cout << "  # Use hexbin for large datasets to show density" << std::endl;
        std::cout << "  plt.hexbin(df_plot['GC_Content'], df_plot['Median_Kmer_Frequency'], gridsize=100, cmap='viridis', bins='log', xscale='linear', yscale='log')" << std::endl;
        std::cout << "  plt.xlabel('GC Content (Density Analogue)')" << std::endl;
        std::cout << "  plt.ylabel('Median K-mer Frequency (Abundance Analogue - Log Scale)')" << std::endl;
        std::cout << "  plt.title('Virtual Ultracentrifugation Plot')" << std::endl;
        std::cout << "  plt.colorbar(label='Log(Number of Reads in Bin)')" << std::endl;
        std::cout << "  plt.grid(True, alpha=0.3)" << std::endl;
        std::cout << "  plt.savefig('" << output_prefix << "_gc_abundance_plot.png', dpi=150)" << std::endl;
        std::cout << "  # plt.show()" << std::endl;


        std::cout << "\nAnalysis finished." << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "\nAn error occurred: " << e.what() << std::endl;
        // Consider adding cleanup for temporary files/dirs here if needed
        return 1;
    } catch (...) {
        std::cerr << "\nAn unknown error occurred." << std::endl;
        // Consider adding cleanup for temporary files/dirs here if needed
        return 1;
    }

    // Optional: Clean up KMC temporary files and database
    // Note: Removing directories might require platform-specific code or filesystem library
    // std::remove(kmc_input_list.c_str());
    // std::remove((kmc_db_name + ".kmc_pre").c_str());
    // std::remove((kmc_db_name + ".kmc_suf").c_str());
    // std::filesystem::remove_all(kmc_tmp_dir); // C++17 way to remove directory

    return 0;
}
