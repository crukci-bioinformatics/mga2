use anyhow::{bail, Result};
use log::{info, warn};
use mga2::fastq::{create_fastq_reader, create_fastq_writer, FastqRecord};
use std::path::PathBuf;
use structopt::StructOpt;

/// Trim sequences from FASTQ files and split into chunks in FASTQ and/or FASTA
/// format
#[derive(StructOpt)]
struct Config {
    /// FASTQ file to sample records from
    #[structopt(parse(from_os_str))]
    fastq_files: Vec<PathBuf>,

    /// The prefix for output file names
    #[structopt(short, long)]
    output_prefix: String,

    /// The chunk size for splitting records into separate files
    #[structopt(short = "n", long, default_value = "1000000")]
    chunk_size: u32,

    /// The start position in reads (preceeding bases will be trimmed)
    #[structopt(short, long)]
    start: Option<usize>,

    /// The maximum length of trimmed sequences
    #[structopt(short, long)]
    length: Option<usize>,

    /// Output FASTA files in addition to FASTQ chunks
    /// (note that these will not be trimmed)
    #[structopt(long)]
    output_fasta: bool,
}

fn main() -> Result<()> {
    env_logger::init();

    let config = Config::from_args();

    if config.fastq_files.is_empty() {
        bail!("No input FASTQ files specified");
    }

    if config.chunk_size == 0 {
        bail!("Invalid chunk size");
    }

    if let Some(start) = config.start {
        if start == 0 {
            bail!("Invalid start position for trimming - numbering is 1-based");
        }
    }

    if let Some(length) = config.length {
        if length == 0 {
            bail!("Zero length specified for trimmed sequences");
        }
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

fn trim_and_split(
    fastq_files: &[PathBuf],
    output_prefix: &str,
    chunk_size: u32,
    start: Option<usize>,
    length: Option<usize>,
    output_fasta: bool,
) -> Result<()> {
    let mut count = 0;
    let mut filtered_count = 0;

    let mut output_file_count = 1;

    let output_fastq_file = format!("{}.{}.fq", output_prefix, output_file_count);
    let output_fasta_file = format!("{}.{}.fa", output_prefix, output_file_count);

    info!("Creating {}", output_fastq_file);
    let output_fastq_path = Some(PathBuf::from(output_fastq_file));
    let mut fastq_writer = create_fastq_writer(&output_fastq_path)?;

    let mut fasta_writer = if output_fasta {
        info!("Creating {}", output_fasta_file);
        let output_fasta_path = Some(PathBuf::from(output_fasta_file));
        Some(create_fastq_writer(&output_fasta_path)?)
    } else {
        None
    };

    let mut chunk_count = 0;

    for fastq_file in fastq_files {
        let filename = fastq_file.to_str().unwrap();
        info!("Reading {}", filename);

        let mut fastq_reader = create_fastq_reader(fastq_file)?;

        let mut record = FastqRecord::new();

        while fastq_reader.read_next_into(&mut record)? {
            // check whether we need to create a new chunk file
            if chunk_count > 0 && chunk_count % chunk_size == 0 {
                fastq_writer.flush()?;
                if let Some(ref mut fasta_writer) = fasta_writer {
                    fasta_writer.flush()?
                }

                output_file_count += 1;

                let output_fastq_file = format!("{}.{}.fq", output_prefix, output_file_count);
                let output_fasta_file = format!("{}.{}.fa", output_prefix, output_file_count);

                info!("Creating {}", output_fastq_file);
                let output_fastq_path = Some(PathBuf::from(output_fastq_file));
                fastq_writer = create_fastq_writer(&output_fastq_path)?;

                fasta_writer = if output_fasta {
                    info!("Creating {}", output_fasta_file);
                    let output_fasta_path = Some(PathBuf::from(output_fasta_file));
                    Some(create_fastq_writer(&output_fasta_path)?)
                } else {
                    None
                };
            }

            let untrimmed_record = record.clone();

            match length {
                Some(length) => match start {
                    Some(start) => record.trim(start, Some(start + length - 1))?,
                    None => record.trim_to_length(length)?,
                },
                None => {
                    if let Some(start) = start {
                        record.trim(start, None)?;
                    }
                }
            }

            if record.seq.is_empty() {
                filtered_count += 1;
            } else {
                chunk_count += 1;
                fastq_writer.write_fastq(&record)?;
                if let Some(ref mut fasta_writer) = fasta_writer {
                    fasta_writer.write_fasta(&untrimmed_record)?;
                }
            }

            count += 1;
            if count % 1_000_000 == 0 {
                info!("{}", count);
            }
        }
    }

    if count % 1_000_000 != 0 {
        info!("{}", count);
    }

    info!(
        "{} of {} records filtered after trimming",
        filtered_count, count
    );

    fastq_writer.flush()?;

    if filtered_count == count {
        warn!("No sequences remain after trimming and zero-length sequences filtered");
    }

    Ok(())
}
