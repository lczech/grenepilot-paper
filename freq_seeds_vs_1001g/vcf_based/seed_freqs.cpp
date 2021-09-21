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
#include <cmath>
#include <iomanip>
#include <fstream>
#include <limits>
#include <string>
#include <map>
#include <random>
#include <unordered_map>

using namespace genesis;
using namespace genesis::population;
using namespace genesis::utils;

// using Position = std::pair<std::string, size_t>;
//
// struct Data
// {
//     double seed_freq;
//     std::array<double, 16> flle_freqs;
//     bool in_flle = false;
// };


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
    if( argc != 2 ) {
        // throw std::runtime_error( "need a vcf file for seeds and one for 1001g" );
        throw std::runtime_error( "need a vcf file for seeds" );
    }

    // Seed data and 1001g (one thousand one genomges = otog, because variable names
    // cannot start with numbers) data
    std::string const seed_file = std::string( argv[1] );
    // std::string const otog_file = std::string( argv[2] );

    // Compute frequencies of seeds, by combining all counts of the S1..S8 samples.
    // We do it this way, as then frequencies are weighted by number of reads in each sample.
    // We also collect this first and keep it in memory, as the two vcf files contain different
    // positions (I think because the flowers/leaves were not called with the known variants vcf),
    // and so we want an easy way to find matching and ignore missing positions.
    // std::map<Position, Data> data;

    std::ofstream csv_out( "seed_freqs.csv" );
    csv_out << "CHR\tSNP\tA1\tA2\tMAF\tNCHROBS\n";

    // Iterate the seed VCF file and compute frequencies.
    size_t seed_cnt = 0;
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
        if( tot < 30 || ref == 0 || ref == tot ) {
            continue;
        }

        // Compute and store allele frequency.
        // auto const af = 1.0 - static_cast<double>( ref ) / static_cast<double>( tot );
        auto af = static_cast<double>( ref ) / static_cast<double>( tot );
        if( af > 0.5 ) {
            af = 1.0 - af;
        }

        // if( data.count({ record->get_chromosome(), record->get_position() }) > 0 ) {
        //     LOG_ERR << "pos existis: " << record->get_chromosome() << ":" << record->get_position();
        //     return 1;
        // }
        // data[{ record->get_chromosome(), record->get_position() }].seed_freq = af;

        csv_out << record->get_chromosome() << "\t";
        csv_out << record->get_chromosome() << "_" << record->get_position() << "\t";
        csv_out << record->get_variant( 0 ) << "\t" << record->get_variant( 1 ) << "\t";
        csv_out << af << "\t" << tot << "\n";

        ++seed_cnt;
    }
    LOG_INFO << "processed " << seed_cnt << " seed SNPs";

    LOG_INFO << "finished";
    return 0;
}

