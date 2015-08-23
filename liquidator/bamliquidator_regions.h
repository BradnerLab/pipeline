#ifndef PIPELINE_BAMLIQUIDATORINTERNAL_BAMLIQUIDATOR_REGIONS_H
#define PIPELINE_BAMLIQUIDATORINTERNAL_BAMLIQUIDATOR_REGIONS_H

#include <hdf5.h>
#include <hdf5_hl.h>

#include <boost/algorithm/string.hpp>

#include <map>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

const size_t region_name_length = 64;

// this Region must match exactly the structure in HDF5
// -- see bamliquidator_batch.py function create_regions_table
struct Region
{
  uint32_t bam_file_key;
  char chromosome[64];
  char region_name[region_name_length];
  uint64_t start;
  uint64_t stop;
  char strand;
  uint64_t count;
  double normalized_count;

  bool is_valid(const std::map<std::string, size_t>& chromosome_to_length)
  {
    const auto it = chromosome_to_length.find(chromosome);
    if ( it != chromosome_to_length.end() )
    {
      if ( stop <= it->second )
      {
        return true;
      }
    }

    return false;
  }
};

inline std::ostream& operator<<(std::ostream& os, const Region& r)
{
  os << "bam file key " << r.bam_file_key << ' ' << r.chromosome << ' '
     << r.region_name << ' ' << r.start << " -> " << r.stop << ' '
     << r.strand << ' ' << r.normalized_count;
  return os;
}

// default_strand: optional argument, default _ indicates to use 
//                 gff strand column or . (both) for .bed region file
inline
std::vector<Region> parse_regions(const std::string& region_file_path,
                                  const std::string& region_format,
                                  const unsigned int bam_file_key,
                                  const std::map<std::string, size_t>& chromosome_to_length = std::map<std::string, size_t>(), 
                                  const char default_strand = '_') 
{
  using namespace liquidator;

  unsigned int chromosome_column = 0;
  unsigned int name_column = 0;
  unsigned int start_column = 0;
  unsigned int stop_column = 0;
  unsigned int strand_column = 0;
  unsigned int min_columns = 0;

  if (region_format == "gff")
  {
    chromosome_column = 0;
    name_column = 1;
    start_column = 3;
    stop_column = 4;
    strand_column = 6;
    min_columns = 7;
  }
  else if (region_format == "bed")
  {
    chromosome_column = 0;
    name_column = 3;
    start_column = 1;
    stop_column = 2;
    strand_column = 5;
    min_columns = 3;
  }
  else
  {
    throw std::runtime_error("unsupported region file format (" + region_format + " for " + region_file_path
                             + "), please supply a .gff or .bed file");
  }

  std::ifstream region_file(region_file_path.c_str());
  if (!region_file.is_open())
  {
    throw std::runtime_error("failed to open region_file " + region_file_path);
  }

  std::vector<Region> regions;
  int line_number = 1;
  for(std::string line; std::getline(region_file, line); ++line_number)
  {
    std::vector<std::string> columns;
    boost::split(columns, line, boost::is_any_of("\t"));
    if (columns.size() < min_columns)
    {
      std::stringstream ss;
      ss << "Not enough columns parsing line " << line_number << " '" << line << "' of " << region_file_path;
      throw std::runtime_error(ss.str());
    }
    Region region;
    region.bam_file_key = bam_file_key; 
    copy(region.chromosome,  columns[chromosome_column], sizeof(Region::chromosome));
    if (columns.size() > name_column)
    {
      copy(region.region_name, columns[name_column], sizeof(Region::region_name));
      if (columns[name_column].size() >= sizeof(Region::region_name))
      {
        Logger::warn() << "Truncated region on line " << line_number << " from '" << columns[name_column] << "' to '" << region.region_name << "'";
      }
    }
    else
    {
      copy(region.region_name, "", sizeof(Region::region_name));
    }
    region.start = boost::lexical_cast<uint64_t>(columns[start_column]);
    region.stop  = boost::lexical_cast<uint64_t>(columns[stop_column]);
    if (region.start > region.stop)
    {
      std::swap(region.start, region.stop);
    }

    if (columns.size() > strand_column)
    {
      if (columns[strand_column].size() != 1)
      {
        std::stringstream ss;
        ss << "error parsing strand: '" << columns[strand_column] << "' on line " << line_number;
        throw std::runtime_error(ss.str());
      }
      region.strand = default_strand == '_'
                    ? columns[strand_column][0]
                    : default_strand;
    }
    else
    {
      region.strand = default_strand == '_'
                    ? '.'
                    : default_strand; 
    }
    region.count = 0;
    region.normalized_count = 0.0;

    if (chromosome_to_length.empty() || region.is_valid(chromosome_to_length))
    {
      regions.push_back(region);
    }
    else
    {
      Logger::warn() << "Excluding invalid region on line " << line_number << ": " << region;
    }
  }

  return regions;
}

inline void write(hid_t& file, std::vector<Region>& regions)
{
  const size_t record_size = sizeof(Region);

  size_t record_offset[] = { HOFFSET(Region, bam_file_key),
                             HOFFSET(Region, chromosome),
                             HOFFSET(Region, region_name),
                             HOFFSET(Region, start),
                             HOFFSET(Region, stop),
                             HOFFSET(Region, strand),
                             HOFFSET(Region, count),
                             HOFFSET(Region, normalized_count) };

  size_t field_sizes[] = { sizeof(Region::bam_file_key),
                           sizeof(Region::chromosome),
                           sizeof(Region::region_name),
                           sizeof(Region::start),
                           sizeof(Region::stop),
                           sizeof(Region::strand),
                           sizeof(Region::count),
                           sizeof(Region::normalized_count) };

  herr_t status = H5TBappend_records(file, "region_counts", regions.size(), record_size, record_offset,
                                     field_sizes, regions.data());
  if (status != 0)
  {
    std::stringstream ss;
    ss << "Error appending record, status = " << status;
    throw std::runtime_error(ss.str());
  }
}

#endif
