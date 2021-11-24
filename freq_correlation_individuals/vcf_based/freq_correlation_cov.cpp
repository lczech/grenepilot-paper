/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2021 Lucas Czech

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lczech@carnegiescience.edu>
    Department of Plant Biology, Carnegie Institution For Science
    260 Panama Street, Stanford, CA 94305, USA
*/

#include "genesis/genesis.hpp"

#include <array>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <limits>
#include <string>
#include <map>
#include <unordered_map>

using namespace genesis;
using namespace genesis::population;
using namespace genesis::utils;

using Position = std::pair<std::string, size_t>;

// struct Position
// {
//     std::string chrom;
//     size_t pos;
//
//     bool operator < ( Position const& rhs )
//     {
//         return chrom < rhs.chrom || ( chrom == rhs.chrom && pos < rhs.pos );
//     }
// };

struct Data
{
    double seed_freq;
    std::array<std::array<double, 16>, 4> flle_freqs;
    bool in_flle = false;
};

int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;
    LOG_INFO << "started";

    // We use data from:
    //     /Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/SRA/2018-06-26-ena-seeds/
    //     /Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/SRA/2018-10-01-ena-flowerpools/
    // The called frequency vcfs are here:
    //     /lustre/scratch/lczech/grenepipe-runs/ath-evo-seeds-rerun/annotated/all.vcf.gz
    //     /lustre/scratch/lczech/grenepipe-runs/ath-evo-flowerpools/annotated/all.vcf.gz
    // Which is what we expect as input
    // Still we take them as args, in case they move.
    if( argc != 3 ) {
        throw std::runtime_error( "need a vcf file for seeds and one for flowers and leaves" );
    }

    // Seed data and flower/leaf (flle) data
    std::string const seed_file = std::string( argv[1] );
    std::string const flle_file = std::string( argv[2] );

    // Compute frequencies of seeds, by combining all counts of the S1..S8 samples.
    // We do it this way, as then frequencies are weighted by number of reads in each sample.
    // We also collect this first and keep it in memory, as the two vcf files contain different
    // positions (I think because the flowers/leaves were not called with the known variants vcf),
    // and so we want an easy way to find matching and ignore missing positions.
    std::map<Position, Data> data;

    // Iterate the seed VCF file and compute frequencies.
    size_t seed_cnt = 0;
    size_t seed_skp = 0;
    for( auto record = VcfInputIterator( seed_file ); record; ++record ) {

        // Basic filtering: Only SNPs, PASS, need allelic depth field, biallelic SNPs only.
        if(
            ! record->is_snp() || ! record->pass_filter() ||
            ! record->has_format( "AD" ) || record->get_alternatives_count() != 1
        ) {
            continue;
        }

        // Add up counts.
        size_t ref = 0;
        size_t tot = 0;
        record->assert_format( "AD" );
        for( auto sample : record->get_format_int( "AD" )) {
            ref += sample.get_value_at(0);
            for( size_t i = 0; i < sample.values_per_sample(); ++i ) {
                tot += sample.get_value_at(i);
            }
        }

        // Quality control: Skip SNPs with low coverage, and with no alternatives called in the seed bank.
        if( tot < 50 || ref == 0 || ref == tot ) {
            ++seed_skp;
            continue;
        }

        // Compute and store allele frequency.
        auto const af = static_cast<double>( ref ) / static_cast<double>( tot );
        if( data.count({ record->get_chromosome(), record->get_position() }) > 0 ) {
            LOG_ERR << "pos existis: " << record->get_chromosome() << ":" << record->get_position();
            return 1;
        }
        data[{ record->get_chromosome(), record->get_position() }].seed_freq = af;

        ++seed_cnt;
    }
    LOG_INFO << "processed " << seed_cnt << " seed SNPs";
    LOG_INFO << "skipped   " << seed_skp << " seed SNPs";


    // Map of flower/leaf samples.
    // Their order is different in the vcf, so we need proper mapping :-(

    // sample  origin
    // S1      Flowerpool1001und2
    // S2      Flowerpool50A
    // S3      Flowerpool50B
    // S4      Flowerpool25A
    // S5      Flowerpool25B
    // S6      Flowerpool10B2
    // S7      Flowerpool5B
    // S8      Plantpool100x
    // S9      Plantpool50a
    // S10     Plantpool50b
    // S11     Plantpool25a
    // S12     Plantpool25b
    // S13     Plantpool10b1
    // S14     Plantpool10b2
    // S15     Plantpool5x
    // S16     PlantpoolB12345


    // Get the sample names in the order they are in the vcf.
    std::vector<std::string> sample_names;

    // Iterate the seed VCF file and compute frequencies.
    size_t flle_cnt = 0;
    size_t skip_cnt = 0;
    for( auto record = VcfInputIterator( flle_file ); record; ++record ) {

        if( sample_names.empty() ) {
            sample_names = record.header().get_sample_names();
            LOG_INFO << "flower and leaf sample names: " << utils::join( sample_names );
        }

        // Basic filtering: Only SNPs, PASS, need allelic depth field, biallelic SNPs only.
        if(
            ! record->is_snp() || ! record->pass_filter() ||
            ! record->has_format( "AD" ) || record->get_alternatives_count() != 1
        ) {
            continue;
        }

        // See if the current pos is in the seeds. If not, skip.
        if( data.count({ record->get_chromosome(), record->get_position() }) == 0 ) {
            ++skip_cnt;
            continue;
        }
        data[{ record->get_chromosome(), record->get_position() }].in_flle = true;

        // Add up counts.
        record->assert_format( "AD" );
        size_t smp_cnt = 0;
        for( auto sample : record->get_format_int( "AD" )) {
            size_t ref = sample.get_value_at(0);
            size_t tot = 0;
            for( size_t i = 0; i < sample.values_per_sample(); ++i ) {
                tot += sample.get_value_at(i);
            }

            // Compute and store allele frequency.
            auto af = static_cast<double>( ref ) / static_cast<double>( tot );

            // Quality control: Skip SNPs with low coverage,
            // and with no alternatives called in the seed bank.
            // We do this here by setting af to nan, so that it is not considered in the correlation.
            if( ref == 0 || ref == tot ) {
                af = std::numeric_limits<double>::quiet_NaN();
            }

            auto& entry = data[{ record->get_chromosome(), record->get_position() }];
            if( tot < 25 ) {
                entry.flle_freqs[0][smp_cnt] = std::numeric_limits<double>::quiet_NaN();
                entry.flle_freqs[1][smp_cnt] = std::numeric_limits<double>::quiet_NaN();
                entry.flle_freqs[2][smp_cnt] = std::numeric_limits<double>::quiet_NaN();
                entry.flle_freqs[3][smp_cnt] = std::numeric_limits<double>::quiet_NaN();
            } else if( tot < 50 ) {
                entry.flle_freqs[0][smp_cnt] = af;
                entry.flle_freqs[1][smp_cnt] = std::numeric_limits<double>::quiet_NaN();
                entry.flle_freqs[2][smp_cnt] = std::numeric_limits<double>::quiet_NaN();
                entry.flle_freqs[3][smp_cnt] = std::numeric_limits<double>::quiet_NaN();
            } else if( tot < 75 ) {
                entry.flle_freqs[0][smp_cnt] = af;
                entry.flle_freqs[1][smp_cnt] = af;
                entry.flle_freqs[2][smp_cnt] = std::numeric_limits<double>::quiet_NaN();
                entry.flle_freqs[3][smp_cnt] = std::numeric_limits<double>::quiet_NaN();
            } else if( tot < 100 ) {
                entry.flle_freqs[0][smp_cnt] = af;
                entry.flle_freqs[1][smp_cnt] = af;
                entry.flle_freqs[2][smp_cnt] = af;
                entry.flle_freqs[3][smp_cnt] = std::numeric_limits<double>::quiet_NaN();
            } else {
                entry.flle_freqs[0][smp_cnt] = af;
                entry.flle_freqs[1][smp_cnt] = af;
                entry.flle_freqs[2][smp_cnt] = af;
                entry.flle_freqs[3][smp_cnt] = af;
            }

            ++smp_cnt;
        }

        ++flle_cnt;
    }
    LOG_INFO << "processed " << flle_cnt << " flle SNPs";
    LOG_INFO << "skipped   " << skip_cnt << " flle SNPs";

    // Count how many positions in the data from the seeds were not found in the flle data,
    // and delete them so that they do not occur in the correlation.
    size_t omit_cnt = 0;
    auto itr = data.begin();
    while( itr != data.end() ) {
        if( ! itr->second.in_flle ){
            itr = data.erase(itr);
            ++omit_cnt;
        } else {
            ++itr;
        }
    }
    LOG_INFO << "skipped   " << omit_cnt << " seed SNPs";

    auto smp_to_name = std::unordered_map<std::string, std::string>{
        { "S1",  "Flowerpool1001und2" },
        { "S2",  "Flowerpool50A" },
        { "S3",  "Flowerpool50B" },
        { "S4",  "Flowerpool25A" },
        { "S5",  "Flowerpool25B" },
        { "S6",  "Flowerpool10B2" },
        { "S7",  "Flowerpool5B" },
        { "S8",  "Plantpool100x" },
        { "S9",  "Plantpool50a" },
        { "S10", "Plantpool50b" },
        { "S11", "Plantpool25a" },
        { "S12", "Plantpool25b" },
        { "S13", "Plantpool10b1" },
        { "S14", "Plantpool10b2" },
        { "S15", "Plantpool5x" },
        { "S16", "PlantpoolB12345" }
    };

    // We use a lambda that returns a tranforming rage to select an entry from the flle data.
    // auto select_entry = []( ForwardIterator begin, ForwardIterator end, size_t index ){
    //     return utils::make_transform_range(
    //         [index]( std::array<double, 16> const& flle_freqs ){
    //             return flle_freqs[index];
    //         },
    //         begin, end
    //     );
    // };

    // Order of entries?!
    // S1      S10     S11     S12     S13     S14     S15     S16     S2
    // S3      S4      S5      S6      S7      S8      S9
    // Better look it up...

    // We are lazy today and just copy everything into vectors... could be solved with transforming
    // iterators for speed and mem saving, but we need to be quick with results today! Fix later?!

    // Seed freqs as vector.
    std::vector<double> seed_freqs;
    seed_freqs.reserve( data.size() );
    for( auto const& entry : data ) {
        seed_freqs.push_back( entry.second.seed_freq );
    }
    if( seed_freqs.size() != data.size() ) {
        LOG_ERR << "diff size " << data.size() << " != " << seed_freqs.size();
    }

    // Prep output
    std::ofstream csv_out( "freq_correlation_cov.csv" );
    csv_out << "LHS\tRHS\tPCC_25\tPCC_50\tPCC_75\tPCC_100\n";

    // Sample freqs as vectors.
    std::array<std::array<std::vector<double>, 16>, 4> flle_freqs;
    for( size_t cov = 0; cov < 4; ++cov ) {
        for( size_t i = 0; i < 16; ++i ) {
            flle_freqs[cov][i].reserve( data.size() );
            for( auto const& entry : data ) {
                flle_freqs[cov][i].push_back( entry.second.flle_freqs[cov][i] );
            }
            if( flle_freqs[cov][i].size() != data.size() ) {
                LOG_ERR << "diff size " << data.size() << " != " << flle_freqs[cov][i].size();
            }
        }
    }

    // Compute correlation between seeds and all flle columns.
    for( size_t i = 0; i < 16; ++i ) {
        csv_out << "seed\t" << smp_to_name[sample_names[i]];

        // Go through all coverage levels
        for( size_t cov = 0; cov < 4; ++cov ) {
            csv_out << "\t" << pearson_correlation_coefficient( seed_freqs, flle_freqs[cov][i] );
        }
        csv_out << "\n";
    }

    // Compute pairwise correlation between flowers and leaves
    // (we do all, to keep this simple here...)
    for( size_t i = 0; i < 16; ++i ) {
        for( size_t j = i+1; j < 16; ++j ) {
            csv_out << smp_to_name[sample_names[i]] << "\t" << smp_to_name[sample_names[j]];

            for( size_t cov = 0; cov < 4; ++cov ) {
                csv_out << "\t" << pearson_correlation_coefficient( flle_freqs[cov][i], flle_freqs[cov][j] );
            }
            csv_out << "\n";
        }
    }

    LOG_INFO << "finished";
    return 0;
}
