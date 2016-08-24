#include "fasta_scorer.h"

#include "fasta_reader.h"
#include "fimo_style_printer.h"

#include <fstream>

namespace liquidator
{

void process_fasta(const std::vector<ScoreMatrix>& matrices,
                   const std::string& fasta_file_path,
                   const std::string& output_file_path)
{
    std::ifstream fasta_input(fasta_file_path);
    if (!fasta_input)
    {
        throw std::runtime_error("failed to open " + fasta_file_path);
    }

    std::ofstream output(output_file_path);
    FimoStylePrinter printer(output);

    FastaReader fasta_reader(fasta_input);
    std::string sequence;
    std::string sequence_name;
    while (fasta_reader.next_read(sequence, sequence_name))
    {
        for (const auto& matrix : matrices)
        {
            printer.sequence_name = &sequence_name;
            matrix.score(sequence, printer);
        }
    }
}

}
