// alp_fit — Gumbel-parameter fitter for RepeatMasker nucleotide matrices.
//
// Reads a RepeatMasker matrix file (12-letter alphabet A R G C Y T K M S W N X
// with an optional "# FREQS A .. C .. G .. T .." comment line), extracts the
// 4x4 A/C/G/T core, and runs the ALP library (Sls::AlignmentEvaluer::initGapped)
// to estimate the Gumbel parameters lambda/K plus the finite-size-correction
// parameters (a_I/a_J, b_I/b_J, alpha_I/alpha_J, beta_I/beta_J, sigma, tau).
//
// Output (stdout) is a machine-readable "key<TAB>value[<TAB>error]" block,
// including the 11-column NCBI rmblast_*_values row:
//   {open, extend, INT2_MAX, lambda, K, H, alpha, beta, C, alpha_v, sigma}
// with the RMBlast mapping: alpha = a_J + a_I, H = lambda/alpha,
// beta = b_J + b_I, alpha_v = Alpha_J + Alpha_I, sigma = Sigma.
//
// The full evaluer state (including bootstrap arrays for error propagation)
// can be saved with -out <file> and reloaded via operator>>.
//
// Build: see Makefile (links ../../ALP_1.98/ALP_1.98/ALP_1.98_LIB/cpp/libalp.a)

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

#include "sls_alignment_evaluer.hpp"

static void die(const std::string &msg)
{
    std::cerr << "alp_fit: error: " << msg << "\n";
    exit(1);
}

struct NtMatrix {
    // 4x4 core in A,C,G,T order
    long scores[4][4];
    double freqs[4];       // A,C,G,T background frequencies
    bool have_freqs;
};

// Parse a RepeatMasker matrix file: comment lines (#...), one of which may be
// "# FREQS A 0.295 C 0.205 G 0.205 T 0.295"; then a header line of column
// letters; then one row per letter: "<letter> s1 s2 ...".
static NtMatrix read_rm_matrix(const std::string &path)
{
    std::ifstream in(path.c_str());
    if (!in) die("cannot open matrix file: " + path);

    NtMatrix m;
    m.have_freqs = false;
    for (int i = 0; i < 4; i++) m.freqs[i] = 0.25;

    std::vector<char> col_letters;
    std::map<char, std::vector<long> > rows;
    std::string line;
    while (std::getline(in, line)) {
        // strip leading whitespace
        size_t p = line.find_first_not_of(" \t\r");
        if (p == std::string::npos) continue;
        if (line[p] == '#') {
            std::istringstream ss(line.substr(p + 1));
            std::string tok;
            ss >> tok;
            if (tok == "FREQS") {
                std::string letter;
                double f;
                while (ss >> letter >> f) {
                    switch (toupper(letter[0])) {
                    case 'A': m.freqs[0] = f; break;
                    case 'C': m.freqs[1] = f; break;
                    case 'G': m.freqs[2] = f; break;
                    case 'T': m.freqs[3] = f; break;
                    default: break; // ignore other letters
                    }
                }
                m.have_freqs = true;
            }
            continue;
        }
        std::istringstream ss(line.substr(p));
        if (col_letters.empty()) {
            // header line of single-letter column names
            std::string tok;
            bool ok = true;
            std::vector<char> letters;
            while (ss >> tok) {
                if (tok.size() != 1 || !isalpha((unsigned char)tok[0])) { ok = false; break; }
                letters.push_back((char)toupper(tok[0]));
            }
            if (!ok || letters.empty())
                die("expected header line of column letters in " + path);
            col_letters = letters;
        } else {
            std::string tok;
            ss >> tok;
            if (tok.size() != 1 || !isalpha((unsigned char)tok[0]))
                die("expected row letter at start of line: " + line);
            char rl = (char)toupper(tok[0]);
            std::vector<long> vals;
            long v;
            while (ss >> v) vals.push_back(v);
            if (vals.size() != col_letters.size())
                die("row " + std::string(1, rl) + " has wrong number of scores");
            rows[rl] = vals;
        }
    }

    static const char core[4] = {'A', 'C', 'G', 'T'};
    // column index of each core letter
    int col_idx[4];
    for (int j = 0; j < 4; j++) {
        col_idx[j] = -1;
        for (size_t c = 0; c < col_letters.size(); c++)
            if (col_letters[c] == core[j]) { col_idx[j] = (int)c; break; }
        if (col_idx[j] < 0)
            die(std::string("matrix missing column for letter ") + core[j]);
    }
    for (int i = 0; i < 4; i++) {
        std::map<char, std::vector<long> >::const_iterator it = rows.find(core[i]);
        if (it == rows.end())
            die(std::string("matrix missing row for letter ") + core[i]);
        for (int j = 0; j < 4; j++)
            m.scores[i][j] = it->second[col_idx[j]];
    }
    return m;
}

// Plain ALP-format inputs (first token = alphabet size), as used by alp.exe.
static NtMatrix read_plain(const std::string &matpath, const std::string &freqpath)
{
    NtMatrix m;
    m.have_freqs = false;
    for (int i = 0; i < 4; i++) m.freqs[i] = 0.25;

    std::ifstream in(matpath.c_str());
    if (!in) die("cannot open matrix file: " + matpath);
    long n;
    if (!(in >> n) || n != 4) die("plain matrix must be 4x4 (first token 4)");
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            if (!(in >> m.scores[i][j])) die("bad matrix entry");

    if (!freqpath.empty()) {
        std::ifstream fin(freqpath.c_str());
        if (!fin) die("cannot open freqs file: " + freqpath);
        if (!(fin >> n) || n != 4) die("freqs file must start with 4");
        for (int i = 0; i < 4; i++)
            if (!(fin >> m.freqs[i])) die("bad freqs entry");
        m.have_freqs = true;
    }
    return m;
}

int main(int argc, char **argv)
{
    std::string rm_matrix, plain_matrix, plain_freqs, out_file;
    long gap_open = -1, gap_extend = -1, gap_open2 = -1, gap_extend2 = -1;
    double eps_lambda = 0.01, eps_K = 0.05;
    double max_time = 60.0, max_mem = 2000.0;
    long seed = 1;
    bool iad = false;              // insertions_after_deletions
    bool quiet = false;

    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        std::string v = (i + 1 < argc) ? argv[i + 1] : "";
        if (a == "-matrix")            { rm_matrix = v; i++; }
        else if (a == "-scoremat")     { plain_matrix = v; i++; }
        else if (a == "-freqs1")       { plain_freqs = v; i++; }
        else if (a == "-gapopen")      { gap_open = atol(v.c_str()); i++; }
        else if (a == "-gapextend")    { gap_extend = atol(v.c_str()); i++; }
        else if (a == "-gapopen2")     { gap_open2 = atol(v.c_str()); i++; }
        else if (a == "-gapextend2")   { gap_extend2 = atol(v.c_str()); i++; }
        else if (a == "-eps_lambda")   { eps_lambda = atof(v.c_str()); i++; }
        else if (a == "-eps_K")        { eps_K = atof(v.c_str()); i++; }
        else if (a == "-max_time")     { max_time = atof(v.c_str()); i++; }
        else if (a == "-max_mem")      { max_mem = atof(v.c_str()); i++; }
        else if (a == "-rand")         { seed = atol(v.c_str()); i++; }
        else if (a == "-iad")          { iad = (v == "true" || v == "1"); i++; }
        else if (a == "-out")          { out_file = v; i++; }
        else if (a == "-quiet")        { quiet = true; }
        else die("unknown argument: " + a + "\nusage: alp_fit -matrix <rm.matrix> "
                 "-gapopen N -gapextend N [-gapopen2 N -gapextend2 N] "
                 "[-eps_lambda F] [-eps_K F] [-max_time S] [-max_mem MB] "
                 "[-rand SEED] [-iad true|false] [-out file.par] [-quiet]\n"
                 "   or: alp_fit -scoremat <plain4x4> -freqs1 <plain4> ...");
    }
    if (rm_matrix.empty() && plain_matrix.empty()) die("need -matrix or -scoremat");
    if (gap_open < 0 || gap_extend < 0) die("need -gapopen and -gapextend");
    if (gap_open2 < 0) gap_open2 = gap_open;
    if (gap_extend2 < 0) gap_extend2 = gap_extend;

    NtMatrix m = rm_matrix.empty() ? read_plain(plain_matrix, plain_freqs)
                                   : read_rm_matrix(rm_matrix);
    if (!m.have_freqs)
        std::cerr << "alp_fit: warning: no background frequencies found; using 0.25 each\n";

    const long *rowptr[4];
    for (int i = 0; i < 4; i++) rowptr[i] = m.scores[i];

    Sls::AlignmentEvaluer evaluer;
    // Time-budgeted computation, LAST-style temperature retry loop.
    int attempt = 0;
    for (;; ++attempt) {
        double t = Sls::default_importance_sampling_temperature + 0.01 * attempt;
        try {
            evaluer.initGapped(4, rowptr, m.freqs, m.freqs,
                               gap_open, gap_extend, gap_open2, gap_extend2,
                               iad, eps_lambda, eps_K,
                               max_time, max_mem, seed, t);
            break;
        } catch (const Sls::error &e) {
            std::cerr << "alp_fit: attempt " << attempt << " (temperature " << t
                      << ") failed: " << e.st << " (code " << e.error_code << ")\n";
            if (attempt == 20) die("initGapped failed after 21 temperature attempts");
        }
    }
    if (!evaluer.isGood()) die("evaluer not in good state after initGapped");

    const Sls::ALP_set_of_parameters &p = evaluer.parameters();

    // Derived NCBI rmblast_*_values row (see blast_stat.c -RMH- comment):
    double alpha = p.a_J + p.a_I;
    double H = (alpha != 0.0) ? p.lambda / alpha : 0.0;
    double beta = p.b_J + p.b_I;
    double alpha_v = p.alpha_J + p.alpha_I;

    std::cout.precision(10);
    std::cout << std::fixed;
    std::cout << "matrix\t" << (rm_matrix.empty() ? plain_matrix : rm_matrix) << "\n";
    std::cout << "gap_open\t" << gap_open << "\n";
    std::cout << "gap_extend\t" << gap_extend << "\n";
    std::cout << "gap_open2\t" << gap_open2 << "\n";
    std::cout << "gap_extend2\t" << gap_extend2 << "\n";
    std::cout << "insertions_after_deletions\t" << (iad ? "true" : "false") << "\n";
    std::cout << "seed\t" << seed << "\n";
    std::cout << "temperature_attempts\t" << (attempt + 1) << "\n";
    std::cout << "freqs_ACGT\t" << m.freqs[0] << "\t" << m.freqs[1] << "\t"
              << m.freqs[2] << "\t" << m.freqs[3] << "\n";

    std::cout << "lambda\t" << p.lambda << "\t" << p.lambda_error << "\n";
    std::cout << "K\t" << p.K << "\t" << p.K_error << "\n";
    std::cout << "C\t" << p.C << "\t" << p.C_error << "\n";
    std::cout << "a_I\t" << p.a_I << "\t" << p.a_I_error << "\n";
    std::cout << "a_J\t" << p.a_J << "\t" << p.a_J_error << "\n";
    std::cout << "b_I\t" << p.b_I << "\t" << p.b_I_error << "\n";
    std::cout << "b_J\t" << p.b_J << "\t" << p.b_J_error << "\n";
    std::cout << "alpha_I\t" << p.alpha_I << "\t" << p.alpha_I_error << "\n";
    std::cout << "alpha_J\t" << p.alpha_J << "\t" << p.alpha_J_error << "\n";
    std::cout << "beta_I\t" << p.beta_I << "\t" << p.beta_I_error << "\n";
    std::cout << "beta_J\t" << p.beta_J << "\t" << p.beta_J_error << "\n";
    std::cout << "sigma\t" << p.sigma << "\t" << p.sigma_error << "\n";
    std::cout << "tau\t" << p.tau << "\t" << p.tau_error << "\n";
    std::cout << "gapless_a\t" << p.gapless_a << "\t" << p.gapless_a_error << "\n";
    std::cout << "gapless_alpha\t" << p.gapless_alpha << "\t" << p.gapless_alpha_error << "\n";
    std::cout << "calc_time\t" << p.m_CalcTime << "\n";

    // NCBI 11-column row, ready to paste
    std::cout << "ncbi_row\t" << gap_open << "\t" << gap_extend << "\tINT2_MAX\t"
              << p.lambda << "\t" << p.K << "\t" << H << "\t"
              << alpha << "\t" << beta << "\t"
              << p.C << "\t" << alpha_v << "\t" << p.sigma << "\n";

    if (!out_file.empty()) {
        std::ofstream out(out_file.c_str());
        if (!out) die("cannot open output file: " + out_file);
        out << evaluer;
        if (!quiet)
            std::cerr << "alp_fit: full evaluer state written to " << out_file << "\n";
    }
    return 0;
}
