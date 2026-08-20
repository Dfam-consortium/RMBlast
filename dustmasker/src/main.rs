//! `dustmasker` — minimal, drop-in replacement for NCBI's dustmasker application,
//! built on the DUST implementation in `rmblast-lib` (a faithful port of NCBI's
//! `CSymDustMasker`).
//!
//! SCOPE: this covers the FASTA-in / text-out subset that legacy RepeatMasker- and
//! RepeatScout-era tooling actually calls, e.g.
//!
//!   dustmasker -in seqs.fa -outfmt fasta | ...
//!
//! Supported: `-in -out -window -level -linker -infmt fasta
//!            -outfmt {interval,fasta,acclist} -hard_masking`.
//!
//! Deliberately NOT ported (each exits non-zero with a clear message rather than
//! silently producing something subtly different):
//!   * `-infmt blastdb`
//!   * `-outfmt seqloc_*` / `maskinfo_*` (ASN.1 / XML)
//!   * `-parse_seqids`
//!
//! Sequence identifiers are echoed back exactly as they appeared in the input
//! defline.  NCBI instead re-generates them through its object manager, which for
//! plain FASTA without `-parse_seqids` yields the same string; with
//! `-parse_seqids` it emits a normalised Seq-id (`lcl|foo `), which is why that
//! flag is refused here.

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::process::ExitCode;

use anyhow::{bail, Context, Result};
use clap::Parser;

use rmblast_lib::encoding::BLASTNA_TO_IUPAC;
use rmblast_lib::filter::dust::{dustmasker_intervals, DUST_LEVEL, DUST_LINKER, DUST_WINDOW};
use rmblast_lib::seq::FastaReader;

// ──────────────────────────────────────────────────────────────────────────────
// Argument-normalization shim for NCBI single-dash compatibility
// ──────────────────────────────────────────────────────────────────────────────
//
// NCBI's C++ Toolkit accepts long options with a single leading dash
// (`-outfmt fasta`); clap requires two (`--outfmt fasta`).  Rewrite
// `-<letter><rest>` to `--<letter><rest>` before clap sees argv.  Bare `-`
// (stdin), `--`, and single-character options pass through unchanged.
//
// Same shim as rmblastn/src/main.rs; single-character flags cannot be bundled.
fn normalize_args() -> Vec<std::ffi::OsString> {
    std::env::args_os()
        .map(|arg| {
            let s = match arg.to_str() {
                Some(s) => s,
                None => return arg,
            };
            if s.starts_with("--") || s == "-" {
                return arg;
            }
            if let Some(rest) = s.strip_prefix('-') {
                let mut chars = rest.chars();
                if let Some(first) = chars.next() {
                    if first.is_ascii_alphabetic() && chars.next().is_some() {
                        return format!("--{}", rest).into();
                    }
                }
            }
            arg
        })
        .collect()
}

#[derive(Parser, Debug)]
#[command(
    name = "dustmasker",
    version,
    about = "Low complexity region masker based on Symmetric DUST algorithm",
    after_help = "\
SCOPE:
  This is the FASTA-in / text-out subset of NCBI dustmasker.  Unsupported
  options (-infmt blastdb, -outfmt seqloc_*/maskinfo_*, -parse_seqids) exit
  with an error rather than producing near-miss output.

  Both single-dash (NCBI style) and double-dash long options are accepted.
  Sequence identifiers are echoed verbatim from the input defline.
"
)]
struct Args {
    /// input file name
    #[arg(long, default_value = "-")]
    r#in: String,

    /// output file name
    #[arg(long, default_value = "-")]
    out: String,

    /// DUST window length
    #[arg(long, default_value_t = DUST_WINDOW)]
    window: usize,

    /// DUST level (score threshold for subwindows)
    #[arg(long, default_value_t = DUST_LEVEL)]
    level: u32,

    /// DUST linker (how close masked intervals should be to get merged together)
    #[arg(long, default_value_t = DUST_LINKER)]
    linker: usize,

    /// input format (possible values: fasta)
    #[arg(long, default_value = "fasta")]
    infmt: String,

    /// output format (possible values: interval, fasta, acclist)
    #[arg(long, default_value = "interval")]
    outfmt: String,

    /// Use hard masking (N) for fasta outfmt instead of lowercase
    #[arg(long = "hard_masking", alias = "hard-masking", default_value_t = false)]
    hard_masking: bool,

    /// Parse Seq-ids in FASTA input (NOT PORTED)
    #[arg(long = "parse_seqids", alias = "parse-seqids", default_value_t = false)]
    parse_seqids: bool,
}

#[derive(Clone, Copy, PartialEq)]
enum OutFmt {
    Interval,
    Fasta,
    Acclist,
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("dustmasker: {:#}", e);
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<()> {
    let args = Args::parse_from(normalize_args());

    if args.parse_seqids {
        bail!(
            "-parse_seqids is not ported.  This dustmasker echoes the input \
             defline verbatim; it does not reproduce NCBI's Seq-id \
             normalisation (e.g. 'lcl|foo ').  Re-run without -parse_seqids."
        );
    }

    match args.infmt.as_str() {
        "fasta" => {}
        "blastdb" => bail!("-infmt blastdb is not ported; supply the sequences as FASTA"),
        other => bail!("unknown input format '{}' (supported: fasta)", other),
    }

    let outfmt = match args.outfmt.as_str() {
        "interval" => OutFmt::Interval,
        "fasta" => OutFmt::Fasta,
        "acclist" => OutFmt::Acclist,
        fmt @ ("seqloc_asn1_bin" | "seqloc_asn1_text" | "seqloc_xml" | "maskinfo_asn1_bin"
        | "maskinfo_asn1_text" | "maskinfo_xml") => bail!(
            "-outfmt {} is not ported (ASN.1/XML serialisation); \
             supported formats are interval, fasta, acclist",
            fmt
        ),
        other => bail!(
            "unknown output format '{}' (supported: interval, fasta, acclist)",
            other
        ),
    };

    // NCBI raises this same error (dust_mask_app.cpp:x_GetWriter).
    if args.hard_masking && outfmt != OutFmt::Fasta {
        bail!("Hard masking can only be applied for fasta output");
    }

    let input: Box<dyn BufRead> = if args.r#in == "-" {
        Box::new(BufReader::new(io::stdin()))
    } else {
        Box::new(BufReader::new(
            File::open(&args.r#in).with_context(|| format!("opening input '{}'", args.r#in))?,
        ))
    };
    let mut output: Box<dyn Write> = if args.out == "-" {
        Box::new(BufWriter::new(io::stdout()))
    } else {
        Box::new(BufWriter::new(
            File::create(&args.out).with_context(|| format!("creating output '{}'", args.out))?,
        ))
    };

    let mut reader = FastaReader::new(input);
    while let Some(rec) = reader.next_record()? {
        // NCBI skips zero-length bioseqs.
        if rec.is_empty() {
            continue;
        }
        let bases = rec.bases();
        let masks = dustmasker_intervals(bases, args.window, args.level, args.linker);

        match outfmt {
            OutFmt::Interval => {
                writeln!(output, ">{}", rec.defline)?;
                for &(s, e) in &masks {
                    writeln!(output, "{} - {}", s, e)?;
                }
            }
            OutFmt::Acclist => {
                for &(s, e) in &masks {
                    writeln!(output, ">{}\t{}\t{}", rec.defline, s, e)?;
                }
            }
            OutFmt::Fasta => {
                writeln!(output, ">{}", rec.defline)?;
                write_masked_fasta(&mut output, bases, &masks, args.hard_masking)?;
            }
        }
    }

    output.flush()?;
    Ok(())
}

/// Mirror of `CMaskWriterFasta::Print`: IUPAC sequence in 60-column lines with
/// masked positions either lowercased (soft, the default) or replaced by `N`.
fn write_masked_fasta(
    out: &mut dyn Write,
    bases: &[u8],
    masks: &[(usize, usize)],
    hard_masking: bool,
) -> Result<()> {
    let mut line = Vec::with_capacity(60);
    let mut imask = 0usize;

    for (i, &b) in bases.iter().enumerate() {
        // masks are sorted and non-overlapping; advance past any that ended.
        while imask < masks.len() && i > masks[imask].1 {
            imask += 1;
        }
        let masked = imask < masks.len() && i >= masks[imask].0;

        let letter = BLASTNA_TO_IUPAC[(b & 15) as usize];
        line.push(if masked {
            if hard_masking { b'N' } else { letter.to_ascii_lowercase() }
        } else {
            letter
        });

        if line.len() == 60 {
            out.write_all(&line)?;
            out.write_all(b"\n")?;
            line.clear();
        }
    }
    if !line.is_empty() {
        out.write_all(&line)?;
        out.write_all(b"\n")?;
    }
    Ok(())
}
