#ifndef LIQUIDATOR_FASTA_READER_H_INCLUDED
#define LIQUIDATOR_FASTA_READER_H_INCLUDED

#include <istream>

namespace liquidator
{

class FastaReader
{
public:
    inline FastaReader(std::istream& fasta_file);

    inline bool next_read(std::string& sequence, std::string& sequence_name, char& strand);

private:
    std::istream& m_fasta_file;
};

FastaReader::FastaReader(std::istream &fasta_file)
    : m_fasta_file(fasta_file)
{}

bool FastaReader::next_read(std::string &sequence, std::string &sequence_name, char& strand)
{
    if (!m_fasta_file)
    {
        throw std::runtime_error("can't read next sequence because istream is not good");
    }
    return m_fasta_file;
}

}

#endif
