//! Functions for reading and writing FASTQ records.
//! 
//! Example of reading a FASTQ file, iterating over each record and counting
//! the cumulative number of bases. This creates a new
//! [`FastqRecord`](struct.FastqRecord.html) struct for each record.
//! 
//! ```
//! # use anyhow::Result;
//! use mga2::fastq::FastqReader;
//! use std::io::{self, BufReader};
//!
//! # fn main() -> Result<()> {
//! let mut reader = FastqReader::new(BufReader::new(io::stdin()));
//!
//! let mut number_of_records = 0;
//! let mut number_of_bases = 0;
//!
//! while let Some(record) = reader.read_next()? {
//!     number_of_records += 1;
//!     number_of_bases += record.seq.len();
//! }
//!
//! println!("Number of reads: {}", number_of_records);
//! println!("Number of bases: {}", number_of_bases);
//! # Ok(())
//! # }
//! ```
//! 
//! A single FastqRecord struct can be reused for each iteration to avoid the
//! cost of creating new instances.
//! 
//! ```
//! # use anyhow::Result;
//! use mga2::fastq::{FastqReader, FastqRecord};
//! use std::io::{self, BufReader};
//! 
//! # fn main() -> Result<()> {
//! let mut reader = FastqReader::new(BufReader::new(io::stdin()));
//! let mut record = FastqRecord::new();
//! 
//! let mut number_of_records = 0;
//! let mut number_of_bases = 0;
//!
//! while reader.read_next_into(&mut record)? {
//!     number_of_records += 1;
//!     number_of_bases += record.seq.len();
//! }
//! # Ok(())
//! # }
//! ```

use anyhow::{bail, Context, Result};
use flate2::bufread::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{stdout, BufRead, BufReader, BufWriter, ErrorKind, Read, Write};
use std::path::{Path, PathBuf};

#[derive(Clone, Default)]
pub struct FastqRecord {
    pub id: String,
    pub desc: Option<String>,
    pub seq: String,
    pub qual: String,
}

impl FastqRecord {
    pub fn new() -> FastqRecord {
        FastqRecord {
            id: String::new(),
            desc: None,
            seq: String::new(),
            qual: String::new(),
        }
    }

    pub fn with_capacity(id_capacity: usize, sequence_capacity: usize) -> FastqRecord {
        FastqRecord {
            id: String::with_capacity(id_capacity),
            desc: None,
            seq: String::with_capacity(sequence_capacity),
            qual: String::with_capacity(sequence_capacity),
        }
    }

    pub fn clear(&mut self) {
        self.id.clear();
        self.desc = None;
        self.seq.clear();
        self.qual.clear();
    }

    pub fn is_empty(&mut self) -> bool {
        self.seq.is_empty()
    }

    pub fn check(&self) -> Result<()> {
        if self.id.is_empty() {
            bail!("Missing identifier for FASTQ record");
        }

        if !self.seq.is_ascii() {
            bail!("Sequence string contains non-ASCII character(s)");
        }

        if !self.qual.is_ascii() {
            bail!("Quality score string contains non-ASCII character(s)");
        }

        if self.seq.len() != self.qual.len() {
            bail!("Sequence and quality strings have different lengths");
        }

        Ok(())
    }

    pub fn trim(&mut self, start: usize, end: Option<usize>) -> Result<()> {
        let length = self.seq.len();
        if self.qual.len() != length {
            bail!("Attempt to trim FASTQ record with differing sequence and quality lengths");
        }

        if start == 0 {
            bail!("Attempt to trim FASTQ record using invalid zero start position (numbering is 1-based)");
        }

        if let Some(end) = end {
            if end < start {
                bail!("Attempt to trim FASTQ record with end < start");
            }
        }

        if start > length {
            self.seq = String::new();
            self.qual = String::new();
        } else {
            match end {
                Some(end) => {
                    self.seq = self.seq.as_str()[(start - 1)..end.min(length)].to_string();
                    self.qual = self.qual.as_str()[(start - 1)..end.min(length)].to_string();
                }
                None => {
                    self.seq = self.seq.as_str()[(start - 1)..].to_string();
                    self.qual = self.qual.as_str()[(start - 1)..].to_string();
                }
            }
        }

        Ok(())
    }

    pub fn trim_to_length(&mut self, length: usize) -> Result<()> {
        if self.seq.len() != self.qual.len() {
            bail!("Attempt to trim FASTQ record with differing sequence and quality lengths");
        }
        if length < self.seq.len() {
            self.seq = self.seq.as_str()[..length].to_string();
            self.qual = self.qual.as_str()[..length].to_string();
        }
        Ok(())
    }
}

pub struct FastqReader<R: BufRead> {
    reader: R,
    name: String,
    id_capacity: usize,
    sequence_capacity: usize,
    test_separator_char: [u8; 1],
    line_count: u32,
}

impl<R: BufRead> FastqReader<R> {
    pub fn new(reader: R) -> Self {
        FastqReader::with_capacity(reader, "unnamed", 50, 160)
    }

    pub fn with_name(reader: R, name: &str) -> Self {
        FastqReader::with_capacity(reader, name, 50, 160)
    }

    pub fn with_capacity(
        reader: R,
        name: &str,
        id_capacity: usize,
        sequence_capacity: usize,
    ) -> Self {
        FastqReader {
            reader,
            name: name.to_string(),
            id_capacity,
            sequence_capacity,
            test_separator_char: [0],
            line_count: 0,
        }
    }

    fn read_next_line(&mut self, buffer: &mut String) -> usize {
        let number_of_bytes = self.reader.read_line(buffer).unwrap();
        self.line_count += 1;
        number_of_bytes
    }

    pub fn read_next(&mut self) -> Result<Option<FastqRecord>> {
        let mut record = FastqRecord::with_capacity(self.id_capacity, self.sequence_capacity);
        if self.read_next_into(&mut record)? {
            Ok(Some(record))
        } else {
            Ok(None)
        }
    }

    pub fn read_next_into(&mut self, record: &mut FastqRecord) -> Result<bool> {
        record.clear();

        loop {
            match self.reader.read_exact(&mut self.test_separator_char) {
                Ok(_) => {
                    if self.test_separator_char[0] == b'@' {
                        break;
                    } else if self.test_separator_char[0] == b'\n' {
                        // skip empty line
                        continue;
                    } else {
                        bail!(
                            "expected '@' character at beginning of line {}, {}",
                            self.line_count + 1,
                            self.name
                        );
                    }
                }
                Err(error) => {
                    if error.kind() == ErrorKind::UnexpectedEof {
                        return Ok(false);
                    } else {
                        bail!(
                            "unexpected problem reading line {}, {}",
                            self.line_count + 1,
                            self.name
                        );
                    }
                }
            }
        }

        self.read_next_line(&mut record.id);

        match record.id.find(' ') {
            Some(length) => {
                record.desc = Some(record.id[length + 1..record.id.len() - 1].to_string());
                record.id.truncate(length);
            }
            None => {
                record.id.pop();
            }
        }

        if record.id.is_empty() {
            bail!(
                "error extracting sequence id for record at line {}, {}",
                self.line_count,
                self.name
            );
        }

        let mut length = 0;

        loop {
            let number_of_bytes = self.read_next_line(&mut record.seq);
            match number_of_bytes {
                0 => bail!(
                    "incomplete record at line {}, {}",
                    self.line_count - 1,
                    self.name
                ),
                1 => {
                    // empty line
                    record.seq.pop();
                }
                _ => {
                    // check for separator character
                    let first_char = &record.seq[length..length + 1];
                    if first_char == "+" {
                        // ignore whatever else was on the '+' line that separates the sequence and quality scores
                        record.seq.truncate(length);
                        break;
                    }
                    // remove newline character
                    record.seq.pop();
                    // check for spaces in sequence
                    if record.seq[length..].find(' ').is_some() {
                        bail!(
                            "sequence contains space(s) for record at line {}, {}",
                            self.line_count,
                            self.name
                        );
                    }
                    length = record.seq.len();
                }
            };
        }

        let mut qual_length = 0;
        while qual_length < length {
            let number_of_bytes = self.read_next_line(&mut record.qual);
            match number_of_bytes {
                0 => bail!(
                    "incomplete record at line {}, {}",
                    self.line_count - 1,
                    self.name
                ),
                _ => {
                    record.qual.pop();
                }
            };
            qual_length = record.qual.len();
        }

        if record.qual.len() != length {
            bail!(
                "sequence and quality lengths differ for record ending at line {}, {}",
                self.line_count,
                self.name
            );
        }

        // check for spaces in quality scores
        if record.qual.find(' ').is_some() {
            bail!(
                "qualities contains space(s) for record ending at line {}, {}",
                self.line_count,
                self.name
            );
        }

        Ok(true)
    }
}

pub struct FastqWriter<W: Write> {
    writer: BufWriter<W>,
}

impl<W: Write> FastqWriter<W> {
    pub fn new(writer: BufWriter<W>) -> Self {
        FastqWriter { writer }
    }

    pub fn write_fastq(&mut self, record: &FastqRecord) -> Result<()> {
        write_fastq(&mut self.writer, record)
            .with_context(|| format!("Error writing FASTA format for record: {}", record.id))?;
        Ok(())
    }

    pub fn write_fasta(&mut self, record: &FastqRecord) -> Result<()> {
        write_fasta(&mut self.writer, record)
            .with_context(|| format!("Error writing FASTA format for record: {}", record.id))?;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()?;
        Ok(())
    }
}

pub fn create_fastq_reader(fastq_file: &Path) -> Result<FastqReader<BufReader<Box<dyn Read>>>> {
    let fastq_file_name = match fastq_file.to_str() {
        Some(name) => String::from(name),
        None => {
            bail!("invalid file name for {:?}", fastq_file)
        }
    };

    let file = File::open(fastq_file)
        .with_context(|| format!("Error opening file {}", fastq_file_name))?;

    let reader: Box<dyn Read> = if fastq_file.to_str().unwrap().ends_with(".gz") {
        Box::new(MultiGzDecoder::new(BufReader::with_capacity(
            64 * 1024,
            file,
        )))
    } else {
        Box::new(file)
    };

    let buffered_reader = BufReader::with_capacity(64 * 1024, reader);
    let fastq_reader = FastqReader::with_name(buffered_reader, fastq_file_name.as_str());

    Ok(fastq_reader)
}

pub fn create_fastq_writer(output_file: &Option<PathBuf>) -> Result<FastqWriter<Box<dyn Write>>> {
    let writer: Box<dyn Write> = match output_file {
        Some(output_file) => {
            let output_filename = output_file.to_str().unwrap();
            let file = File::create(output_file)?;
            if output_filename.ends_with(".gz") {
                Box::new(GzEncoder::new(file, Compression::fast()))
            } else {
                Box::new(file)
            }
        }
        None => Box::new(stdout()),
    };

    let buffered_writer = BufWriter::new(writer);
    let fastq_writer = FastqWriter::new(buffered_writer);

    Ok(fastq_writer)
}

fn write_fastq(writer: &mut dyn Write, record: &FastqRecord) -> Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(record.id.as_bytes())?;
    if let Some(desc) = &record.desc {
        writer.write_all(b" ")?;
        writer.write_all(desc.as_bytes())?;
    }
    writer.write_all(b"\n")?;
    writer.write_all(record.seq.as_bytes())?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(record.qual.as_bytes())?;
    writer.write_all(b"\n")?;
    Ok(())
}

const FASTA_WIDTH: usize = 80;

fn write_fasta(writer: &mut dyn Write, record: &FastqRecord) -> Result<()> {
    writer.write_all(b">")?;
    writer.write_all(record.id.as_bytes())?;
    if let Some(desc) = &record.desc {
        writer.write_all(b" ")?;
        writer.write_all(desc.as_bytes())?;
    }
    writer.write_all(b"\n")?;
    let sequence = record.seq.as_bytes();
    let length = sequence.len();
    let mut position = 0;
    while position < length {
        writer.write_all(&sequence[position..length.min(position + FASTA_WIDTH)])?;
        writer.write_all(b"\n")?;
        position += FASTA_WIDTH;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const RECORD_ID: &str = "MDE123";
    const DESCRIPTION: &str = "a sample read for testing";
    const SEQUENCE: &str = "TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA";
    const QUALITIES: &str = "AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ";

    const EMPTY_RECORD: &[u8] = b"";

    const FASTQ_RECORD: &[u8] = b"@MDE123 a sample read for testing
TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA
+
AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
";

    const INCOMPLETE_RECORD: &[u8] = b"@MDE123 a sample read for testing
TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA
+
";

    fn create_record() -> FastqRecord {
        FastqRecord {
            id: RECORD_ID.to_string(),
            desc: Some(DESCRIPTION.to_string()),
            seq: SEQUENCE.to_string(),
            qual: QUALITIES.to_string(),
        }
    }

    #[test]
    fn valid_record() {
        let mut record = create_record();
        assert!(record.check().is_ok(), "Invalid record");
        record.desc = None;
        assert!(record.check().is_ok(), "Invalid record");
    }

    #[test]
    fn missing_id() {
        let mut record = FastqRecord::new();
        record.desc = Some(DESCRIPTION.to_string());
        record.seq = SEQUENCE.to_string();
        record.qual = QUALITIES.to_string();
        let result = record.check();
        assert!(result.is_err(), "Expecting error for missing identifier");
        let error = result.unwrap_err();
        assert_eq!(error.to_string(), "Missing identifier for FASTQ record");
    }

    #[test]
    fn invalid_sequence() {
        let mut record = FastqRecord::new();
        record.id = RECORD_ID.to_string();
        record.seq = "TGTGACCCAAGAAGTTGTTAAAATéTCCGGAGGTAGCCATTATATACCAA".to_string();
        record.qual = QUALITIES.to_string();
        let result = record.check();
        assert!(
            result.is_err(),
            "Expecting error for non-ASCII character in sequence"
        );
        let error = result.unwrap_err();
        assert_eq!(
            error.to_string(),
            "Sequence string contains non-ASCII character(s)"
        );
    }

    #[test]
    fn invalid_quality_string() {
        let mut record = FastqRecord::new();
        record.id = RECORD_ID.to_string();
        record.seq = SEQUENCE.to_string();
        record.qual = "AAFFFJJJJJJJJJJJJJJJIJJéJJJJJJJJJJJJJJJJJJJJJJJJJJ".to_string();
        let result = record.check();
        assert!(
            result.is_err(),
            "Expecting error for non-ASCII character in quality"
        );
        let error = result.unwrap_err();
        assert_eq!(
            error.to_string(),
            "Quality score string contains non-ASCII character(s)"
        );
    }

    #[test]
    fn differing_sequence_and_quality_strings() {
        let mut record = FastqRecord::new();
        record.id = RECORD_ID.to_string();
        record.seq = SEQUENCE.to_string();
        record.qual = "AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJ".to_string();
        let result = record.check();
        assert!(
            result.is_err(),
            "Expecting error for differing lengths of sequence and quality strings"
        );
        let error = result.unwrap_err();
        assert_eq!(
            error.to_string(),
            "Sequence and quality strings have different lengths"
        );
    }

    #[test]
    fn read_empty_record() {
        let mut reader = FastqReader::new(EMPTY_RECORD);
        let result = reader.read_next();
        assert!(result.is_ok(), "Error reading empty FASTQ record");
        let record = result.unwrap();
        assert!(record.is_none(), "Record found when none expected");
    }

    #[test]
    fn read_single_record() {
        let mut reader = FastqReader::new(FASTQ_RECORD);
        let result = reader.read_next();
        assert!(result.is_ok(), "Error reading FASTQ record");
        let record = result.unwrap();
        assert!(record.is_some(), "No record read");
        let record = record.unwrap();
        assert_eq!(record.id, RECORD_ID.to_string());
        assert_eq!(record.desc, Some(DESCRIPTION.to_string()));
        assert_eq!(record.seq, SEQUENCE.to_string());
        assert_eq!(record.qual, QUALITIES.to_string());
        assert!(record.check().is_ok(), "Invalid record");
        let result = reader.read_next();
        assert!(result.is_ok(), "Error reading FASTQ record");
        let record = result.unwrap();
        assert!(record.is_none(), "Record found when none expected");
    }

    #[test]
    fn read_into_record() {
        let mut reader = FastqReader::new(FASTQ_RECORD);
        let mut record = FastqRecord::new();
        let result = reader.read_next_into(&mut record);
        assert!(result.is_ok(), "Error reading FASTQ record");
        assert_eq!(record.id, RECORD_ID.to_string());
        assert_eq!(record.desc, Some(DESCRIPTION.to_string()));
        assert_eq!(record.seq, SEQUENCE.to_string());
        assert_eq!(record.qual, QUALITIES.to_string());
        assert!(record.check().is_ok(), "Invalid record");
        let result = reader.read_next_into(&mut record);
        assert!(result.is_ok(), "Error reading FASTQ record");
        let record_found = result.unwrap();
        assert!(!record_found, "Record found when none expected");
    }

    #[test]
    fn read_incomplete_record() {
        let mut reader = FastqReader::new(INCOMPLETE_RECORD);
        let mut record = FastqRecord::new();
        let result = reader.read_next_into(&mut record);
        assert!(result.is_err());
        let error = result.unwrap_err();
        assert!(error.to_string().starts_with("incomplete record at line"));
    }

    const INCOMPLETE_SECOND_RECORD: &[u8] = b"@MDE123 a sample read for testing
TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA
+
AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@
";

    #[test]
    fn incomplete_second_record() {
        let mut reader = FastqReader::new(INCOMPLETE_SECOND_RECORD);
        let mut record = FastqRecord::new();
        let result = reader.read_next_into(&mut record);
        assert!(result.is_ok(), "Error reading valid first FASTQ record");
        assert!(record.check().is_ok(), "Expecting valid first record");
        let result = reader.read_next_into(&mut record);
        assert!(
            result.is_err(),
            "Expecting error for incomplete second record"
        );
        let error = result.unwrap_err();
        assert!(error
            .to_string()
            .starts_with("error extracting sequence id for record at line"));
    }

    #[test]
    fn write_fastq_record() {
        let record = create_record();
        let mut writer = FastqWriter::new(BufWriter::new(Vec::new()));
        writer
            .write_fastq(&record)
            .expect("Error writing FASTQ record");
        writer.flush().expect("Error flushing FASTQ writer");
        assert_eq!(writer.writer.get_ref(), &FASTQ_RECORD);
    }

    const FASTA_RECORD: &[u8] = b">MDE123 a sample read for testing
TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA
";

    #[test]
    fn write_fasta_record() {
        let record = create_record();
        let mut writer = FastqWriter::new(BufWriter::new(Vec::new()));
        writer
            .write_fasta(&record)
            .expect("Error writing FASTA record");
        writer.flush().expect("Error flushing FASTQ writer");
        assert_eq!(writer.writer.get_ref(), &FASTA_RECORD);
    }

    #[test]
    fn trim_start_and_end() {
        let mut record = create_record();
        let result = record.trim(11, Some(35));
        assert!(result.is_ok(), "Unexpected error during trimming");
        assert!(record.check().is_ok(), "Invalid record after trimming");
        assert_eq!(record.seq, "GAAGTTGTTAAAATTTCCGGAGGTA".to_string());
    }

    #[test]
    fn trim_start() {
        let mut record = create_record();
        let result = record.trim(11, None);
        assert!(result.is_ok(), "Unexpected error during trimming");
        assert!(record.check().is_ok(), "Invalid record after trimming");
        assert_eq!(record.seq, "GAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA".to_string());
        assert_eq!(record.qual, "JJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ".to_string());
    }

    #[test]
    fn trim_beyond_end() {
        let mut record = create_record();
        let result = record.trim(11, Some(1000));
        assert!(result.is_ok(), "Unexpected error during trimming");
        assert!(record.check().is_ok(), "Invalid record after trimming");
        assert_eq!(record.seq, "GAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA".to_string());
        assert_eq!(record.qual, "JJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ".to_string());
    }
    
    #[test]
    fn trim_to_length() {
        let mut record = create_record();
        let result = record.trim_to_length(36);
        assert!(result.is_ok(), "Unexpected error during trimming");
        assert!(record.check().is_ok(), "Invalid record after trimming");
        assert_eq!(record.seq, "TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAG".to_string());
        assert_eq!(record.qual, "AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJ".to_string());
    }
}
