#ifndef LIQUIDATOR_BAM_SCORER_H_INCLUDED
#define LIQUIDATOR_BAM_SCORER_H_INCLUDED

#include "bamliquidator_regions.h"
#include "score_matrix.h"

#include <samtools/bam.h>

#include <boost/filesystem.hpp>
#include <tbb/concurrent_queue.h>

#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <mutex>
#include <deque>

static const int MAX_QUEUED_READS = 200;
static const int MAX_THREAD_CHUNK = 100000;

namespace liquidator
{

inline bool unmapped(const bam1_t& read)
{
    // The 3rd bit being set in the flag means it is unmapped.
    // Binary with 3rd bit set (0b100) is 4.
    static const uint32_t unmapped_bit = 4;
    return read.core.flag & unmapped_bit;
}

inline bool reverse_complemented(const bam1_t& read)
{
    // The 5th bit being set in the flag means it is reverse complemented.
    // Binary with 5th bit set (0b10000) is 16.
    static const uint32_t reverse_complemented_bit = 16;
    return read.core.flag & reverse_complemented_bit;
}

class BamScorer
{   
public:
    enum PrintStyle
    {
        None,
        Fimo,
        MappedFimo
    };

    BamScorer(const std::string& bam_input_file_path,
              const std::vector<ScoreMatrix>& matrices,
              PrintStyle print_style,
              bool only_score_unmapped,
              const std::string& bam_output_file_path,
              const std::string& region_file_path = "")
    :
        m_input(bam_open(bam_input_file_path.c_str(), "r")),
        m_bam_output_file_path(bam_output_file_path),
        m_output(0),
        m_header(bam_header_read(m_input)),
        m_index(bam_index_load(bam_input_file_path.c_str())),
        m_matrices(matrices),
        m_print_style(print_style),
        m_only_score_unmapped(only_score_unmapped),
        m_read(0),
        m_read_count(0),
        m_unmapped_count(0),
        m_read_hit_count(0),
        m_unmapped_hit_count(0),
        m_total_hit_count(0)
    {
        if (m_input == 0 || m_header == 0 || m_index == 0)
        {
            throw std::runtime_error("failed to open " + bam_input_file_path);
        }

        if (!bam_output_file_path.empty())
        {
            m_output = bam_open(bam_output_file_path.c_str(), "w");
            bam_header_write(m_output, m_header);
        }

        if (m_print_style != None)
        {
            std::cout << "#pattern name\tsequence name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched sequence" << std::endl;
        }

        m_reading = true;
        std::thread score_thread1(&BamScorer::score_thread_actual, this);
        //std::thread score_thread2(&BamScorer::score_thread_actual, this);
        //std::thread score_thread3(&BamScorer::score_thread_actual, this);
        //std::thread score_thread4(&BamScorer::score_thread_actual, this);
        if (!region_file_path.empty())
        {
            printf("REGIONS\n");
            score_regions(region_file_path);
        }
        else
        {
            printf("ALL READS\n");
            score_all_reads();
        }
        m_reading = false;
        score_thread1.join();
        //score_thread2.join();
        //score_thread3.join();
        //score_thread4.join();
    }

    ~BamScorer()
    {
        auto print_percent = [](const std::string& upper_label, size_t upper_value, const std::string& lower_label, size_t lower_value) {
            std::cout << "# (" << upper_label << ") / (" << lower_label << ") = " << upper_value << '/' << lower_value << " = " << 100*(double(upper_value)/lower_value) << '%' << std::endl;
        };

        if (!m_only_score_unmapped)
        {
            print_percent("reads hit", m_read_hit_count, "total reads", m_read_count);
            print_percent("mapped hit", m_read_hit_count - m_unmapped_hit_count, "mapped reads", m_read_count - m_unmapped_count);
        }
        print_percent("unmapped hit", m_unmapped_hit_count, "unmapped reads", m_unmapped_count);
        if (!m_only_score_unmapped)
        {
            print_percent("unmapped hit", m_unmapped_hit_count, "total hit", m_read_hit_count);
        }
        print_percent("unmapped reads", m_unmapped_count, "total reads", m_read_count);
        std::cout << "# total hits: " << m_total_hit_count << " (average hits per hit read = " << double(m_total_hit_count)/m_read_hit_count << ")" << std::endl;

        bam_index_destroy(m_index);
        bam_header_destroy(m_header);
        bam_close(m_input);

        if (m_output)
        {
            bam_close(m_output);

            // Suprisingly the bam seems already sorted, even if regions provided
            // in such a way that bam_write1 I expect is called out of order.
            // At some point I expect that I will need to do some kind of explicit sort.
            int index_rc = bam_index_build(m_bam_output_file_path.c_str());
            if (index_rc != 0)
            {
                std::cerr << "Failed to build index for " << m_bam_output_file_path << std::endl;
            }
        }
    }

    void operator()(const std::string& motif_name,
                    size_t start,
                    size_t stop,
                    const ScoreMatrix::Score& score)
    {
        if (score.pvalue() < 0.0001)
        {
            ++m_total_hit_count;
            if (m_print_style != None)
            {
                std::cout << motif_name << '\t';
                if (m_print_style == MappedFimo)
                {
                    const char* chromosome = m_read->core.tid >= 0 ? m_header->target_name[m_read->core.tid] : "*";
                    std::cout << (unmapped(*m_read) ? "un" : "") << "mapped:" << chromosome << ":" << (char*) m_read->data << '\t'
                              << m_read->core.pos + start << '\t'
                              << m_read->core.pos + stop << '\t'
                              << (score.is_reverse_complement() ? '-' : '+') << '\t';
                }
                else // m_print_style == Fimo
                {
                    // fimo reading a fasta has no way of determining if a sequence is reverse or forward mapping,
                    // so if read is reverse complemented then we should reverse the direction to get perfect fimo matching.
                    bool fimo_reverse = score.is_reverse_complement();
                    size_t fimo_start = start;
                    size_t fimo_stop  = stop;
                    if (reverse_complemented(*m_read))
                    {
                        fimo_reverse = !fimo_reverse;
                        fimo_start = m_read->core.l_qseq - stop  + 1;
                        fimo_stop  = m_read->core.l_qseq - start + 1;
                    }
                    std::cout << (char*) m_read->data << '\t'
                              << fimo_start << '\t'
                              << fimo_stop << '\t'
                              << (fimo_reverse ? '-' : '+') << '\t';
                }

                std::cout.precision(6);
                std::cout << score.score() << '\t';
                std::cout.precision(3);

                std::cout << score.pvalue() << '\t'
                          << '\t' // omit q-value for now
                          << score << std::endl;
            }
        }
    }

private:

    struct BamAllocator
    {
        BamAllocator() { }

        ~BamAllocator()
        {
            // note: we are not calling bam_destroy1(), which was just free()ing the
            // bam and the bam data. since the vectors manage the bam memory,
            // we just need to free the bam dat, but if samtools changes this will
            // likely break
            free(bam.data);
        }

        // note: we are not calling bam_init1(), which was just calloc()ing,
        // if that ever changes this will likely break
        bam1_t bam;
    };
    
    typedef std::shared_ptr < std::vector < BamAllocator > > BamVecPtr;
    
    void score_thread_actual()
    {
        BamVecPtr vec;

        int scored_total = 0;
        while(m_reading || m_queued_reads.size() > 0)
        {
            m_queue_mutex.lock();
            if(m_queued_reads.empty())
            {
                m_queue_mutex.unlock();
                std::this_thread::yield();
                continue;
            }

            vec = m_queued_reads.front();
            m_queued_reads.pop_front();
            
            m_queue_mutex.unlock();

            for(size_t i = 0; i < vec->size(); i++)
            {
                scored_total++;
                BamAllocator &raii_read = (*vec)[i];
                score_read(&raii_read.bam);
            }

            //printf("scored: %d\n", scored_total);
        }
    }
    
    void score_all_reads()
    {
        // todo: the unmapped reads seem to all be at the very end of the loop.
        //       to speed up scoring just the unmapped reads, we could probably skip to the last indexed read and start there.
        //       although, that might be relying on undocumented behavior that could change in future releases, so maybe that is a bad idea.
        //       also, there seems to be some mechanism for storing unmapped reads that correspond to a chromosome, so that is probably a doubly bad idea.
        //       see https://www.biostars.org/p/86405/#86439

        BamVecPtr vec(new std::vector < BamAllocator >);
        size_t vec_count = 0;
        vec->resize(MAX_THREAD_CHUNK);
        memset(&((*vec)[0]), 0, MAX_THREAD_CHUNK * sizeof(BamAllocator));

        int total_read = 0;
        int vec_total = 0;
        while(1)
        {
            bool break_main = false;
            
            m_queue_mutex.lock();
            if(m_queued_reads.size() > MAX_QUEUED_READS)
            {
                m_queue_mutex.unlock();
                std::this_thread::yield();
                continue;
            }
            m_queue_mutex.unlock();
            
            while(vec_count < MAX_THREAD_CHUNK)
            {
                int ret;
                ret = bam_read1(m_input, &(*vec)[vec_count].bam);
                if(ret < 0)
                {
                    printf("EOF? %d\n", total_read);
                    break_main = true;
                    break;
                }
                
                vec_count++;
                total_read++;
            }

            // since the vector is allocated in big chunks, we need to remove the
            // tail end of the vector when we reach EOF
            if(vec_count != vec->size())
                vec->erase(vec->begin() + vec_count, vec->end());
            vec_total += vec->size();
            m_queue_mutex.lock();
            m_queued_reads.push_back(vec);
            m_queue_mutex.unlock();
            
            vec.reset(new std::vector < BamAllocator >);
            vec_count = 0;
            vec->resize(MAX_THREAD_CHUNK);
            memset(&((*vec)[0]), 0, MAX_THREAD_CHUNK * sizeof(BamAllocator));

            if(break_main)
            {
                printf("Breaking main loop\n");
                break;
            }            
        }

        printf("vec total: %d\n", vec_total);
    }

    void score_regions(const std::string& region_file_path)
    {
        const std::string region_extension = boost::filesystem::extension(region_file_path).erase(0, 1);
        for (const Region& region : parse_regions(region_file_path, region_extension, 0))
        {
            // todo: don't I just need to parse_region once per chromosome to get the tid? perhaps there is a faster way to do this without parsing a whole region string?
            // todo: consider adding a util function to do this and remove duplicate code in bamliquidator.cpp
            std::stringstream coord;
            coord << region.chromosome << ':' << region.start << '-' << region.stop;

            int ref,beg,end;
            const int region_parse_rc = bam_parse_region(m_header, coord.str().c_str(), &ref, &beg, &end);
            if (region_parse_rc != 0)
            {
                std::stringstream error_msg;
                error_msg << "bam_parse_region failed with return code " << region_parse_rc;
                throw std::runtime_error(error_msg.str());
            }
            if(ref<0)
            {
                // this bam doesn't have this chromosome
                continue;
            }

            const int fetch_rc = bam_fetch(m_input, m_index, ref, beg, end, this, bam_fetch_func);
            if (fetch_rc != 0)
            {
                std::stringstream error_msg;
                error_msg << "bam_fetch failed with return code " << fetch_rc;
                throw std::runtime_error(error_msg.str());
            }
        }
    }

    //void score_read(std::shared_ptr < BamAllocator > raii_read)
    void score_read(const bam1_t *read)
    {
        ++m_read_count;
        if (unmapped(*read))
        {
            ++m_unmapped_count;
        }
        else if (m_only_score_unmapped)
        {
            return;
        }

        const bam1_core_t *c = &read->core;
        uint8_t *s = bam1_seq(read);

        // [s, s+c->l_qseq) is the sequence, with two bases packed into each byte.
        // I bet we could directly search that instead of first copying into a string
        // but lets get something simple working first. An intermediate step could be
        // to search integers without using bam_nt16_rev_table (and I wouldn't have
        // to worry about the packing complexity).

        if (m_sequence.size() != size_t(c->l_qseq))
        {
            // assuming that all reads are uniform length, this will only happen once
            m_sequence = std::string(c->l_qseq, ' ');
        }
        for (int i = 0; i < c->l_qseq; ++i)
        {
            m_sequence[i] = bam_nt16_rev_table[bam1_seqi(s, i)];
        }
        
        const size_t hit_count_before_this_read = m_total_hit_count;
        m_read = read;
        for (const auto& matrix : m_matrices)
        {
            matrix.score(m_sequence, *this);
        }
        if (m_total_hit_count > hit_count_before_this_read)
        {
            ++m_read_hit_count;
            if (unmapped(*read))
            {
                ++m_unmapped_hit_count;
            }
            if (m_output)
            {
                // TODO: if we want to run multiple score threads in parallel, we will need to
                // add an output queue similar to input queue here
                bam_write1(m_output, read);
            }
        }
    }

    void fetch_helper(const bam1_t *read)
    {
        /*std::shared_ptr < BamAllocator > raii_read(new BamAllocator);
        bam_copy1(raii_read->bam, read);

        while(m_queued_reads.size() >= MAX_QUEUED_READS)
            std::this_thread::yield();
        
            m_queued_reads.push_back(raii_read);*/
    }

    static int bam_fetch_func(const bam1_t* read, void* handle)
    {
        BamScorer& scorer = *static_cast<BamScorer*>(handle);
        scorer.score_read(read);
        
        //scorer.fetch_helper(read);        
        
        return 0;
    }

private:
    bamFile m_input;
    const std::string m_bam_output_file_path;
    bamFile m_output;
    bam_header_t* m_header;
    bam_index_t* m_index;
    const std::vector<ScoreMatrix>& m_matrices;
    const PrintStyle m_print_style;
    const bool m_only_score_unmapped;
    const bam1_t* m_read;
    size_t m_read_count;
    size_t m_unmapped_count;
    size_t m_read_hit_count;
    size_t m_unmapped_hit_count;
    size_t m_total_hit_count;
    std::string m_sequence;
    std::deque < BamVecPtr > m_queued_reads;
    std::mutex m_queue_mutex;
    bool m_reading;
};

}

#endif
