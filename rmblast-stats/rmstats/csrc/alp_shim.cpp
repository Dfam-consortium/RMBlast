// extern "C" shim over Sls::AlignmentEvaluer (ALP library) for the rmstats
// "alp-fit" feature. All Sls::error exceptions are caught at this boundary
// and translated to a status code plus message.

#include <cstring>
#include <string>

#include "sls_alignment_evaluer.hpp"

extern "C" {

// Mirrors the fields of Sls::ALP_set_of_parameters that rmstats consumes.
// Must match #[repr(C)] struct AlpFitRaw in src/alp.rs.
struct AlpFitRaw {
    double lambda, lambda_error;
    double k, k_error;
    double c, c_error;
    double a_i, a_i_error;
    double a_j, a_j_error;
    double b_i, b_i_error;
    double b_j, b_j_error;
    double alpha_i, alpha_i_error;
    double alpha_j, alpha_j_error;
    double beta_i, beta_i_error;
    double beta_j, beta_j_error;
    double sigma, sigma_error;
    double tau, tau_error;
    double gapless_a, gapless_a_error;
    double gapless_alpha, gapless_alpha_error;
    double calc_time;
};

// Run Sls::AlignmentEvaluer::initGapped and extract the fitted parameters.
//
// scores: row-major alphabet_size x alphabet_size matrix.
// freqs1/freqs2: background frequencies per sequence (length alphabet_size).
// If deterministic != 0, install fixed sample counts via
// set_gapped_computation_parameters_simplified(max_time) and call initGapped
// with max_time <= 0 (LAST-style reproducible mode); otherwise run the plain
// time-budgeted mode with max_time.
//
// Returns 0 on success; 1 on Sls::error (message copied to err_buf); 2 on
// any other exception; 3 if the evaluer is not in a good state afterwards.
int rmstats_alp_init_gapped(
    long alphabet_size,
    const long* scores,
    const double* freqs1,
    const double* freqs2,
    long gap_open1, long gap_epen1,
    long gap_open2, long gap_epen2,
    int insertions_after_deletions,
    double eps_lambda, double eps_k,
    double max_time, double max_mem,
    long rand_seed,
    double temperature,
    int deterministic,
    AlpFitRaw* out,
    char* err_buf, long err_buf_len)
{
    if (err_buf && err_buf_len > 0) err_buf[0] = '\0';
    try {
        std::vector<const long*> rows(alphabet_size);
        for (long i = 0; i < alphabet_size; i++)
            rows[i] = scores + i * alphabet_size;

        Sls::AlignmentEvaluer evaluer;
        double init_max_time = max_time;
        if (deterministic) {
            evaluer.set_gapped_computation_parameters_simplified(max_time);
            init_max_time = 0.0;
        }
        evaluer.initGapped(alphabet_size, rows.data(), freqs1, freqs2,
                           gap_open1, gap_epen1, gap_open2, gap_epen2,
                           insertions_after_deletions != 0,
                           eps_lambda, eps_k,
                           init_max_time, max_mem, rand_seed, temperature);
        if (!evaluer.isGood())
            return 3;

        const Sls::ALP_set_of_parameters& p = evaluer.parameters();
        out->lambda = p.lambda;               out->lambda_error = p.lambda_error;
        out->k = p.K;                         out->k_error = p.K_error;
        out->c = p.C;                         out->c_error = p.C_error;
        out->a_i = p.a_I;                     out->a_i_error = p.a_I_error;
        out->a_j = p.a_J;                     out->a_j_error = p.a_J_error;
        out->b_i = p.b_I;                     out->b_i_error = p.b_I_error;
        out->b_j = p.b_J;                     out->b_j_error = p.b_J_error;
        out->alpha_i = p.alpha_I;             out->alpha_i_error = p.alpha_I_error;
        out->alpha_j = p.alpha_J;             out->alpha_j_error = p.alpha_J_error;
        out->beta_i = p.beta_I;               out->beta_i_error = p.beta_I_error;
        out->beta_j = p.beta_J;               out->beta_j_error = p.beta_J_error;
        out->sigma = p.sigma;                 out->sigma_error = p.sigma_error;
        out->tau = p.tau;                     out->tau_error = p.tau_error;
        out->gapless_a = p.gapless_a;         out->gapless_a_error = p.gapless_a_error;
        out->gapless_alpha = p.gapless_alpha; out->gapless_alpha_error = p.gapless_alpha_error;
        out->calc_time = p.m_CalcTime;
        return 0;
    } catch (const Sls::error& e) {
        if (err_buf && err_buf_len > 0) {
            strncpy(err_buf, e.st.c_str(), (size_t)err_buf_len - 1);
            err_buf[err_buf_len - 1] = '\0';
        }
        return 1;
    } catch (...) {
        if (err_buf && err_buf_len > 0) {
            strncpy(err_buf, "unknown exception", (size_t)err_buf_len - 1);
            err_buf[err_buf_len - 1] = '\0';
        }
        return 2;
    }
}

double rmstats_alp_default_temperature(void)
{
    return Sls::default_importance_sampling_temperature;
}

} // extern "C"
