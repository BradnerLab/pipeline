#ifndef LIQUIDATOR_FASTA_SCORER_H_INCLUDED
#define LIQUIDATOR_FASTA_SCORER_H_INCLUDED

#include "score_matrix.h"

#include <string>
#include <vector>

namespace liquidator
{

void process_fasta(const std::vector<ScoreMatrix>& matrices,
                   const std::string& fasta_file_path,
                   const std::string& output_file_path);

}

#endif
