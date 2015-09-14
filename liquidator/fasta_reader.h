#ifndef LIQUIDATOR_FASTA_READER_H_INCLUDED
#define LIQUIDATOR_FASTA_READER_H_INCLUDED

#include <istream>

namespace liquidator
{

class FastaReader
{
public:
    inline FastaReader(std::istream& fasta_file);

    inline bool next_read(std::string& sequence, std::string& sequence_name);

private:
    std::istream& m_fasta_file;
};

FastaReader::FastaReader(std::istream &fasta_file)
    : m_fasta_file(fasta_file)
{}

bool FastaReader::next_read(std::string& sequence, std::string& sequence_name)
{
    // this function assumes that every sequence is preceded by a ">" line, with no empty lines allowed.
    // this is an approximation of the fasta format.
    if (!m_fasta_file) return false;
    std::getline(m_fasta_file, sequence_name);
    if (!m_fasta_file) return false;
    if (sequence_name.empty() || sequence_name[0] != '>')
    {

        throw std::runtime_error("fasta sequence description is missing");
    }
    sequence_name.erase(0, 1); // remove the '>' from the name
    std::getline(m_fasta_file, sequence);
    return m_fasta_file;
}

}

#endif
