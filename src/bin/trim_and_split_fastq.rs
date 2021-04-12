//! Trim sequences from FASTQ files and split into chunks in FASTQ and/or FASTA
//! format.

use anyhow::{ensure, Context, Result};
use log::info;
use mga2::fastq::{FastqReader, FastqRecord, FastqWriter};
use std::path::PathBuf;
use structopt::StructOpt;

/// Configuration parameters specified as command-line options.
#[derive(StructOpt)]
#[structopt(
    about = "Trim sequences from FASTQ files and split into chunks in FASTQ and/or FASTA format."
)]
struct Config {
    /// Input FASTQ files containing sequences to be trimmed and split into
    /// chunks.
    #[structopt(parse(from_os_str))]
    fastq_files: Vec<PathBuf>,

    /// The prefix for output file names.
    #[structopt(short, long)]
    output_prefix: String,

    /// The chunk size for splitting records into separate files.
    #[structopt(short = "n", long, default_value = "1000000")]
    chunk_size: u32,

    /// The start position in reads (preceeding bases will be trimmed).
    #[structopt(short, long)]
    start: Option<usize>,

    /// The maximum length of trimmed sequences.
    #[structopt(short, long)]
    length: Option<usize>,

    /// Output FASTA files in addition to FASTQ chunks
    /// (note that these will not be trimmed).
    #[structopt(long)]
    output_fasta: bool,
}

fn main() -> Result<()> {
    env_logger::init();

    let config = Config::from_args();

    ensure!(
        !config.fastq_files.is_empty(),
        "No input FASTQ files specified"
    );

    ensure!(config.chunk_size > 0, "Invalid chunk size");

    if let Some(start) = config.start {
        ensure!(
            start > 0,
            "Invalid start position for trimming - numbering is 1-based"
        );
    }

    if let Some(length) = config.length {
        ensure!(length > 0, "Invalid length for trimmed sequences");
    }

    trim_and_split(
        &config.fastq_files,
        &config.output_prefix,
        config.chunk_size,
        config.start,
        config.length,
        config.output_fasta,
    )?;

    Ok(())
}

/// Read records from a set of FASTQ files, trimming the sequence and quality
/// strings and writing the resulting trimmed FASTQ records in chunks to files.
fn trim_and_split(
    fastq_files: &[PathBuf],
    output_prefix: &str,
    chunk_size: u32,
    start: Option<usize>,
    length: Option<usize>,
    output_fasta: bool,
) -> Result<()> {
    let mut count = 0;

    let mut output_file_count = 1;
    info!("Creating chunk {}", output_file_count);

    let output_fastq_path = Some(PathBuf::from(format!(
        "{}.{}.fq",
        output_prefix, output_file_count
    )));
    let mut fastq_writer = FastqWriter::to_file(&output_fastq_path)?;

    let mut fasta_writer = if output_fasta {
        let output_fasta_path = Some(PathBuf::from(format!(
            "{}.{}.fa",
            output_prefix, output_file_count
        )));
        Some(FastqWriter::to_file(&output_fasta_path)?)
    } else {
        None
    };

    let mut chunk_count = 0;

    for fastq_file in fastq_files {
        let filename = fastq_file
            .to_str()
            .context("Error obtaining FASTQ file name")?;
        info!("Reading {}", filename);

        let mut fastq_reader = FastqReader::from_file(fastq_file)?;

        let mut record = FastqRecord::new();

        while fastq_reader.read_next_into(&mut record)? {
            // check whether we need to create a new chunk file
            if chunk_count > 0 && chunk_count % chunk_size == 0 {
                fastq_writer.flush()?;
                if let Some(ref mut fasta_writer) = fasta_writer {
                    fasta_writer.flush()?
                }

                output_file_count += 1;
                info!("Creating chunk {}", output_file_count);

                let output_fastq_path = Some(PathBuf::from(format!(
                    "{}.{}.fq",
                    output_prefix, output_file_count
                )));
                fastq_writer = FastqWriter::to_file(&output_fastq_path)?;

                fasta_writer = if output_fasta {
                    let output_fasta_path = Some(PathBuf::from(format!(
                        "{}.{}.fa",
                        output_prefix, output_file_count
                    )));
                    Some(FastqWriter::to_file(&output_fasta_path)?)
                } else {
                    None
                };
            }

            let mut trimmed_record = record.clone();

            match length {
                Some(length) => match start {
                    Some(start) => trimmed_record.trim(start, Some(start + length - 1))?,
                    None => trimmed_record.trim_to_length(length)?,
                },
                None => {
                    if let Some(start) = start {
                        trimmed_record.trim(start, None)?;
                    }
                }
            }

            chunk_count += 1;
            fastq_writer.write_fastq(&trimmed_record)?;
            if let Some(ref mut fasta_writer) = fasta_writer {
                fasta_writer.write_fasta(&record)?;
            }

            count += 1;
            if count % 1_000_000 == 0 {
                info!("{} million records read", count / 1_000_000);
            }
        }
    }

    info!("{} records read", count);

    fastq_writer.flush()?;
    if let Some(ref mut fasta_writer) = fasta_writer {
        fasta_writer.flush()?;
    }

    Ok(())
}
