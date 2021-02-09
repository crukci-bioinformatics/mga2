use anyhow::{bail, Result};
use log::info;
use mga::fastq::{create_fastq_reader, create_fastq_writer, FastqRecord};
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
}

fn main() -> Result<()> {
    env_logger::init();

    let config = Config::from_args();

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
    )?;

    Ok(())
}

fn trim_and_split(
    fastq_files: &[PathBuf],
    output_prefix: &str,
    chunk_size: u32,
    start: Option<usize>,
    length: Option<usize>,
) -> Result<()> {
    let mut count = 0;
    let mut filtered_count = 0;

    let mut output_file_count = 1;
    let output_fastq_file = format!("{}.{}.fq", output_prefix, output_file_count);
    info!("Creating {}", output_fastq_file);
    let output_fastq_path = Some(PathBuf::from(output_fastq_file));
    let mut fastq_writer = create_fastq_writer(&output_fastq_path)?;
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
                output_file_count += 1;
                let output_fastq_file = format!("{}.{}.fq", output_prefix, output_file_count);
                info!("Creating {}", output_fastq_file);
                let output_fastq_path = Some(PathBuf::from(output_fastq_file));
                fastq_writer = create_fastq_writer(&output_fastq_path)?;
            }

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
        bail!("No sequences remain after trimming and zero-length sequences filtered");
    }

    Ok(())
}
