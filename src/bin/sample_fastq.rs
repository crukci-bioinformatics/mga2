//! Sample a specified number of records from a FASTQ file.

use anyhow::{bail, Context, Result};
use log::info;
use mga2::fastq::{create_fastq_reader, create_fastq_writer, FastqRecord};
use rand::{thread_rng, Rng};
use serde::Serialize;
use std::collections::HashSet;
use std::path::{Path, PathBuf};
use structopt::StructOpt;

/// Configuration parameters specified as command-line options.
#[derive(StructOpt)]
#[structopt(about = "Sample a specified number of records from a FASTQ file.")]
struct Config {
    /// Identifier for dataset.
    #[structopt(short, long)]
    id: Option<String>,

    /// FASTQ file to sample records from.
    #[structopt(parse(from_os_str))]
    fastq_files: Vec<PathBuf>,

    /// FASTQ file to which sampled records will be written.
    #[structopt(short, long)]
    output_file: Option<PathBuf>,

    /// The number of records to sample.
    #[structopt(short = "n", long)]
    sample_size: u32,

    /// The maximum number of records to read/sample from.
    #[structopt(short, long)]
    max_number_to_sample_from: Option<u64>,

    /// The minimum length of sequences to include in the sample.
    #[structopt(short = "l", long)]
    min_sequence_length: Option<u32>,

    // Prepend dataset identifier to record identifier.
    #[structopt(short, long)]
    prepend_id: bool,

    // Check record identifiers are unique.
    #[structopt(long)]
    check_unique_record_ids: bool,

    /// Summary of the numbers of records read and sampled.
    #[structopt(short, long)]
    summary_file: Option<PathBuf>,
}

fn main() -> Result<()> {
    env_logger::init();

    let config = Config::from_args();

    if config.fastq_files.is_empty() {
        bail!("No input FASTQ files specified");
    }

    if config.sample_size == 0 {
        bail!("Invalid sample size");
    }

    if let Some(max_number_to_sample_from) = config.max_number_to_sample_from {
        if max_number_to_sample_from == 0 {
            bail!("Invalid maximum number of records to sample from");
        }
    }

    if let Some(min_sequence_length) = config.min_sequence_length {
        if min_sequence_length == 0 {
            bail!("Invalid minimum sequence length");
        }
    }

    let (mut records, total) = sample_fastq(
        &config.fastq_files,
        config.sample_size,
        config.max_number_to_sample_from,
        config.min_sequence_length,
    )?;

    if config.prepend_id {
        if let Some(id) = &config.id {
            for record in records.iter_mut() {
                record.id = format!("{}|{}", id, record.id);
            }
        }
    }

    write_fastq_records(&records, &config.output_file)?;

    if let Some(summary_file) = config.summary_file {
        write_summary(&config.id, total, records.len() as u32, &summary_file)?;
    }

    if config.check_unique_record_ids {
        check_unique_record_ids(&records)?;
    }

    Ok(())
}

/// Read records from a set of FASTQ files and return a randomly selected sample
/// of the specified size.
fn sample_fastq(
    fastq_files: &[PathBuf],
    sample_size: u32,
    max_number_to_sample_from: Option<u64>,
    min_sequence_length: Option<u32>,
) -> Result<(Vec<FastqRecord>, u64)> {
    let mut sampled_records = Vec::new();
    let mut number_of_records_read: u64 = 0;

    let mut rng = thread_rng();

    // loop over and sample from input FASTQ files using reservoir sampling
    'outer: for fastq_file in fastq_files {
        let filename = fastq_file.to_str().unwrap();
        info!("Reading {}", filename);

        let mut reader = create_fastq_reader(fastq_file)?;

        let mut count = 0;
        let mut record = FastqRecord::new();

        while reader.read_next_into(&mut record)? {
            count += 1;
            number_of_records_read += 1;

            if number_of_records_read % 10_000_000 == 0 {
                info!("{}", number_of_records_read);
            }

            if min_sequence_length.is_none()
                || min_sequence_length.unwrap() as usize <= record.seq.len()
            {
                if sampled_records.len() < sample_size as usize {
                    sampled_records.push(record.clone());
                } else {
                    let index = rng.gen_range(0..number_of_records_read);
                    if index < sample_size as u64 {
                        sampled_records[index as usize] = record.clone();
                    }
                }
            }

            if let Some(max_number_to_sample_from) = max_number_to_sample_from {
                if number_of_records_read == max_number_to_sample_from {
                    break 'outer;
                }
            }
        }

        info!("{} records read from {}", count, filename);
    }

    info!("Total records read: {}", number_of_records_read);
    info!("Number of records sampled: {}", sampled_records.len());

    Ok((sampled_records, number_of_records_read))
}

/// Write FASTQ records to a file or stdout.
fn write_fastq_records(records: &[FastqRecord], output_file: &Option<PathBuf>) -> Result<()> {
    match output_file {
        Some(file) => info!(
            "Writing sampled FASTQ records to {}",
            file.to_str().unwrap()
        ),
        None => info!("Writing sampled FASTQ records to stdout"),
    }

    let mut writer = create_fastq_writer(output_file)?;

    // write sampled records
    let mut count = 0;
    for record in records {
        count += 1;
        writer.write_fastq(record)?;
        if count % 10_000_000 == 0 {
            info!("{}", count);
        }
    }

    if count % 10_000_000 != 0 {
        info!("{}", count);
    }

    writer.flush()?;

    Ok(())
}

/// Sampling summary of the number of FASTQ records read and the number sampled.
#[derive(Debug, Serialize)]
struct Summary {
    /// The ID for this sampling run.
    id: String,

    /// The number of FASTQ records read.
    read: u64,

    /// The number of FASTQ records in the sample.
    sampled: u32,
}

/// Write the sampling summary to a CSV file.
fn write_summary(id: &Option<String>, read: u64, sampled: u32, summary_file: &Path) -> Result<()> {
    let summary_filename = summary_file.to_str().unwrap();
    info!("Writing summary to {}", summary_filename);
    let id = match id {
        Some(id) => id.clone(),
        None => String::from(""),
    };
    let summary = Summary { id, read, sampled };
    let mut summary_writer = csv::Writer::from_path(summary_file)
        .with_context(|| format!("Error creating summary file {}", summary_filename))?;
    summary_writer
        .serialize(summary)
        .with_context(|| format!("Error writing summary to {}", summary_filename))?;
    summary_writer.flush().with_context(|| {
        format!(
            "Error writing summary file {} to completion",
            summary_filename
        )
    })?;

    Ok(())
}

/// Check that the given set of FASTQ records have unique identifiers.
fn check_unique_record_ids(records: &[FastqRecord]) -> Result<()> {
    info!("Checking sampled records have unique identifiers");
    let mut ids = HashSet::new();
    for record in records {
        if ids.contains(&record.id) {
            bail!("Duplicate record identifiers found among sampled records");
        }
        ids.insert(record.id.clone());
    }
    Ok(())
}
