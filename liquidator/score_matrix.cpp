#include "score_matrix.h"
#include "detail/score_matrix_detail.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

namespace liquidator
{
constexpr std::array<double, 4> ScoreMatrix::default_acgt_background;

std::array<double, AlphabetSize> ensure_no_zeros_in_background(const std::array<double, AlphabetSize>& original_background)
{
    // I found this behavior in meme 4.11.3.  Comments suggest this is new behavior from Feb 2017.
    // Seems like bad logic, but ensuring bug for bug compatability to ease testing.
    const double small_ammount_to_add = 0.0000005;
    std::array<double, AlphabetSize> bg_no_zeros = original_background;
    double length = 0;
    for (double& bg : bg_no_zeros)
    {
        bg += small_ammount_to_add;
        length += bg;
    }

    // normalize
    for (double& bg : bg_no_zeros)
    {
        bg = bg/length;
    }

    return bg_no_zeros;
}

ScoreMatrix::ScoreMatrix(const std::string& name,
                         const std::array<double, AlphabetSize>& original_background,
                         bool average_background_for_reverse,
                         const std::vector<std::array<double, AlphabetSize>>& pwm,
                         unsigned number_of_sites,
                         bool is_reverse_complement,
                         double pseudo_sites)
    : m_name(name),
      m_is_reverse_complement(is_reverse_complement),
      m_scale(0),
      m_min_before_scaling(0)
{
    // To be honest, I don't understand why the original background would ever be used instead of just using the adjusted background everywhere.
    // This was implemented this way to ensure exact same results as fimo, to ease testing.  This might be bug for bug compatibility.
    const std::array<double, AlphabetSize> adjusted_background = detail::adjust_background(original_background, average_background_for_reverse);
    detail::PWM unscaledPwm {number_of_sites, name, pwm};
    auto min_max = detail::log_adjusted_likelihood_ratio(unscaledPwm, original_background, adjusted_background, pseudo_sites);
    static const unsigned range = 1000;
    auto scaledPWM = detail::scale(unscaledPwm, min_max, range);
    m_matrix = scaledPWM.matrix;
    m_scale = scaledPWM.scale;
    m_min_before_scaling = scaledPWM.min_before_scaling;

    std::vector<double> table = detail::probability_distribution(scaledPWM.matrix, adjusted_background);
    detail::pdf_to_pvalues(table);
    m_pvalues = table;
}

ScoreMatrix::Score
ScoreMatrix::score_sequence(const std::string& sequence, size_t begin, size_t end) const
{
    const unsigned scaled_score = detail::score(m_matrix, sequence, begin, end);
    assert(scaled_score < m_pvalues.size());
    const double pvalue = m_pvalues[scaled_score];
    const double unscaled_score = double(scaled_score)/m_scale + m_matrix.size()*m_min_before_scaling;
    return Score(sequence, m_is_reverse_complement, begin, end, pvalue, unscaled_score);
}

std::vector<ScoreMatrix>
ScoreMatrix::read(std::istream& meme_style_pwm,
                  const std::array<double, AlphabetSize>& acgt_background,
                  const std::string& motif_name,
                  bool include_reverse_complement,
                  double pseudo_sites)
{
    const std::array<double, AlphabetSize>& acgt_background_with_no_zeros = ensure_no_zeros_in_background(acgt_background);

    std::vector<detail::PWM> pwms = detail::read_pwm(meme_style_pwm);
    std::vector<ScoreMatrix> score_matrices;
    for (auto& pwm : pwms)
    {
        if (!motif_name.empty() && pwm.name != motif_name)
        {
            continue;
        }
        score_matrices.push_back(ScoreMatrix(pwm.name, acgt_background_with_no_zeros, include_reverse_complement, pwm.matrix, pwm.number_of_sites, false, pseudo_sites));
        if (include_reverse_complement)
        {
            detail::reverse_complement(pwm.matrix);
            score_matrices.push_back(ScoreMatrix(pwm.name, acgt_background_with_no_zeros, true, pwm.matrix, pwm.number_of_sites, true, pseudo_sites));
        }
    }
    return score_matrices;
}

std::array<double, AlphabetSize>
ScoreMatrix::read_background(std::istream& background_file)
{
    std::array<unsigned, AlphabetSize> counts{ {0, 0, 0, 0} };
    std::array<double, AlphabetSize> background;
    for(std::string line; getline(background_file, line); )
    {
        std::vector<std::string> split;
        boost::split(split, line, boost::is_any_of(" "), boost::token_compress_on);
        if (split.size() != 2 || split[0].size() != 1)
        {
            continue;
        }
        const size_t index = alphabet_index(split[0][0]);
        if (index >= AlphabetSize)
        {
            continue;
        }
        background[index] = boost::lexical_cast<double>(split[1]);
        ++counts[index];
    }
    for (auto& count : counts)
    {
        if (count != 1)
        {
            throw std::runtime_error("failed to read background");
        }
    }
    return background;
}

ScoreMatrix::Score::Score(const std::string &sequence, bool is_reverse_complement, size_t begin, size_t end, double pvalue, double score)
    : m_sequence(sequence),
      m_is_reverse_complement(is_reverse_complement),
      m_begin(begin),
      m_end(end),
      m_pvalue(pvalue),
      m_score(score)
{}

}
