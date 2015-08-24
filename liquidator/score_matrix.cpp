#include "score_matrix.h"

namespace liquidator
{

ScoreMatrix::ScoreMatrix(double min)
:
    m_min(min)
{}

std::vector<ScoreMatrix> read(std::istream& meme_style_pwm,
                              std::array<double, AlphabetSize> acgt_background,
                              double pseudo_sites)
{
    std::vector<detail::PWM> pwms = detail::read_pwm(meme_style_pwm);
    std::vector<ScoreMatrix> score_matrices;
    for (const auto& pwm : pwms)
    {

    }
    return score_matrices;
}

}
