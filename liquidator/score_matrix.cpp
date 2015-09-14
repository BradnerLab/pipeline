#include "score_matrix.h"
#include "detail/score_matrix_detail.h"

namespace liquidator
{

ScoreMatrix::ScoreMatrix(const std::string& name,
                         const std::array<double, AlphabetSize>& background,
                         const std::vector<std::array<double, AlphabetSize>>& pwm,
                         unsigned number_of_sites,
                         bool is_reverse_complement,
                         double pseudo_sites)
    : m_name(name),
      m_is_reverse_complement(is_reverse_complement),
      m_background(background),
      m_scale(0),
      m_min_before_scaling(0)
{
    detail::PWM unscaledPwm {number_of_sites, name, pwm};
    auto min_max = detail::log_adjusted_likelihood_ratio(unscaledPwm, background, pseudo_sites);
    static const unsigned range = 1000;
    auto scaledPWM = detail::scale(unscaledPwm, min_max, range);
    m_matrix = scaledPWM.matrix;
    m_scale = scaledPWM.scale;
    m_min_before_scaling = scaledPWM.min_before_scaling;

    std::vector<double> table = detail::probability_distribution(scaledPWM.matrix, background);
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
                  std::array<double, AlphabetSize> acgt_background,
                  bool include_reverse_complement,
                  double pseudo_sites)
{
    std::vector<detail::PWM> pwms = detail::read_pwm(meme_style_pwm);
    std::vector<ScoreMatrix> score_matrices;
    for (auto& pwm : pwms)
    {
        score_matrices.push_back(ScoreMatrix(pwm.name, acgt_background, pwm.matrix, pwm.number_of_sites, false, pseudo_sites));
        if (include_reverse_complement)
        {
            detail::reverse_complement(pwm.matrix);
            score_matrices.push_back(ScoreMatrix(pwm.name, acgt_background, pwm.matrix, pwm.number_of_sites, true, pseudo_sites));
        }
    }
    return score_matrices;
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
