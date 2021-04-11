//! Functions for reading and writing FASTQ records.
//!
//! # Examples
//!
//! ## Reading FASTQ records from stdin
//!
//! Example of reading FASTQ records from stdin in which a new
//! [`FastqRecord`](struct.FastqRecord.html) struct is created for each record.
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
//! A single FastqRecord struct can be reused for each new record to avoid the
//! cost of allocating memory in each iteration.
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
//!
//! ## Reading and writing FASTQ records
//!
//! The following example reads FASTQ records from stdin, trims sequences to
//! 50 bases and writes the resulting records to stdout.
//!
//! ```
//! # use anyhow::Result;
//! use mga2::fastq::{FastqReader, FastqWriter, FastqRecord};
//! use std::io::{self, BufReader, BufWriter};
//!
//! # fn main() -> Result<()> {
//! let mut reader = FastqReader::new(BufReader::new(io::stdin()));
//! let mut record = FastqRecord::new();
//!
//! let mut writer = FastqWriter::new(BufWriter::new(io::stdout()));
//!
//! while reader.read_next_into(&mut record)? {
//!     record.trim_to_length(50)?;
//!     writer.write_fastq(&record)?;
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

/// A FASTQ record.
///
/// # Example
///
/// ```
/// use mga2::fastq::FastqRecord;
///
/// let record = FastqRecord {
///     id: "MDE123".to_string(),
///     desc: Some("An example FASTQ record".to_string()),
///     seq: "TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA".to_string(),
///     qual: "AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ".to_string(),
/// };
/// ```
#[derive(Clone, Default)]
pub struct FastqRecord {
    /// The sequence identifier - the first part of the first line of the FASTQ
    /// record immediately after the '@' character up to the first space
    /// character.
    pub id: String,

    /// An optional description - the remainder of the first line of the FASTQ
    /// record following the first space character.
    pub desc: Option<String>,

    /// The sequence string.
    pub seq: String,

    /// The quality score string - each character represents a Phred quality
    /// score representing the confidence of the base call.
    pub qual: String,
}

impl FastqRecord {
    /// Create an empty `FastqRecord`.
    ///
    /// # Example
    ///
    /// ```
    /// use mga2::fastq::FastqRecord;
    ///
    /// let mut record = FastqRecord::new();
    /// ```
    pub fn new() -> FastqRecord {
        FastqRecord {
            id: String::new(),
            desc: None,
            seq: String::new(),
            qual: String::new(),
        }
    }

    /// Create an empty `FastqRecord` allocating memory for the expected size of
    /// the record identifier and sequence length. The quality string is
    /// expected to be the same length as the sequence and allocated with the
    /// same capacity.
    ///
    /// # Example
    ///
    /// ```
    /// use mga2::fastq::FastqRecord;
    ///
    /// let mut record = FastqRecord::with_capacity(40, 150);
    /// ```
    pub fn with_capacity(id_capacity: usize, sequence_capacity: usize) -> FastqRecord {
        FastqRecord {
            id: String::with_capacity(id_capacity),
            desc: None,
            seq: String::with_capacity(sequence_capacity),
            qual: String::with_capacity(sequence_capacity),
        }
    }

    /// Clear the contents of the `FastqRecord`.
    ///
    /// # Example
    ///
    /// ```
    /// use mga2::fastq::FastqRecord;
    ///
    /// let mut record = FastqRecord {
    ///     id: "MDE123".to_string(),
    ///     desc: Some("An example FASTQ record".to_string()),
    ///     seq: "TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA".to_string(),
    ///     qual: "AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ".to_string(),
    /// };
    /// assert!(!record.is_empty());
    ///
    /// record.clear();
    /// assert!(record.is_empty());
    /// ```
    pub fn clear(&mut self) {
        self.id.clear();
        self.desc = None;
        self.seq.clear();
        self.qual.clear();
    }

    /// Returns `true` if the sequence is empty.
    ///
    /// # Example
    ///
    /// ```
    /// use mga2::fastq::FastqRecord;
    ///
    /// let mut record = FastqRecord::new();
    /// assert!(record.is_empty());
    ///
    /// record.id = "MDE123".to_string();
    /// record.seq = "TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA".to_string();
    /// record.qual = "AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ".to_string();
    /// assert!(!record.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Checks the validity of the `FastqRecord`.
    ///
    /// # Example
    ///
    /// ```
    /// use mga2::fastq::FastqRecord;
    ///
    /// let mut record = FastqRecord::new();
    ///
    /// record.id = "MDE123".to_string();
    /// record.seq = "AAATéTCCGGA".to_string();
    /// record.qual = "AAFFFJJJJJJ".to_string();
    ///
    /// let result = record.check();
    /// assert!(result.is_err());
    ///
    /// let error = result.unwrap_err();
    /// assert_eq!(error.to_string(), "Sequence string contains non-ASCII character(s)");
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The record identifier is empty
    /// - The sequence or quality score strings contain non-ASCII characters
    /// - The sequence and quality strings have differing lengths
    ///
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

    /// Trims the sequence and quality strings such that the retained portion
    /// starts at the given `start` position and, optionally, ends at the given
    /// `end` position.
    ///
    /// # Example
    ///
    /// ```
    /// use mga2::fastq::FastqRecord;
    ///
    /// let mut record = FastqRecord {
    ///     id: "MDE123".to_string(),
    ///     desc: Some("An example FASTQ record".to_string()),
    ///     seq: "TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA".to_string(),
    ///     qual: "AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ".to_string(),
    /// };
    ///
    /// record.trim(11, Some(30));
    /// assert_eq!(record.seq.len(), 20);
    /// assert_eq!(record.seq, "GAAGTTGTTAAAATTTCCGG");
    /// assert_eq!(record.qual.len(), 20);
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The sequence and quality strings have differing lengths
    /// - The given `start` position is zero - numbering starts from 1
    /// - `end` < `start`
    ///
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

    /// Trims the sequence and quality strings to the given `length`.
    ///
    /// Sequences that are shorter than the specified length are unchanged.
    ///
    /// # Example
    ///
    /// ```
    /// use mga2::fastq::FastqRecord;
    ///
    /// let mut record = FastqRecord {
    ///     id: "MDE123".to_string(),
    ///     desc: Some("An example FASTQ record".to_string()),
    ///     seq: "TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA".to_string(),
    ///     qual: "AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ".to_string(),
    /// };
    ///
    /// record.trim_to_length(30);
    /// assert_eq!(record.seq.len(), 30);
    /// assert_eq!(record.seq, "TGTGACCCAAGAAGTTGTTAAAATTTCCGG");
    /// assert_eq!(record.qual.len(), 30);
    /// assert_eq!(record.qual, "AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJ");
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The sequence and quality strings have differing lengths
    ///
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

/// A reader for FASTQ files.
///
/// # Examples
///
/// ## Reading from stdin
///
/// ```
/// # use anyhow::Result;
/// use mga2::fastq::FastqReader;
/// use std::io::{self, BufReader};
///
/// # fn main() -> Result<()> {
/// let mut reader = FastqReader::new(BufReader::new(io::stdin()));
///
/// let record = reader.read_next();
/// if let Some(record) = reader.read_next()? {
///     println!("{} {}", record.id, record.seq);
/// }
/// # Ok(())
/// # }
/// ```
///
/// ## Reading from a file
///
/// ```
/// # use anyhow::Result;
/// use mga2::fastq::FastqReader;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// # fn main() -> Result<()> {
/// let file = File::open("test.fq")?;
/// let mut reader = FastqReader::new(BufReader::new(file));
///
/// let record = reader.read_next();
/// if let Some(record) = reader.read_next()? {
///     println!("{} {}", record.id, record.seq);
/// }
/// # Ok(())
/// # }
/// ```
///
/// Note that FASTQ files are usually gzip compressed and the
/// [`create_fastq_reader`](fn.create_fastq_reader.html) function will wrap the
/// buffered reader within a decoder if the file name has a ".gz" extension.
/// The `create_fastq_reader` function is the preferred way of creating a
/// `FastqReader` for FASTQ files.
///
pub struct FastqReader<R: BufRead> {
    reader: R,
    name: String,
    id_capacity: usize,
    sequence_capacity: usize,
    test_separator_char: [u8; 1],
    line_count: u32,
}

impl<R: BufRead> FastqReader<R> {
    /// Create a new `FastqReader`.
    ///
    /// ```
    /// use mga2::fastq::FastqReader;
    /// use std::io::{self, BufReader};
    ///
    /// let mut reader = FastqReader::new(BufReader::new(io::stdin()));
    /// ```
    pub fn new(reader: R) -> Self {
        FastqReader::with_capacity(reader, "unnamed", 50, 160)
    }

    /// Create a new 'named' `FastqReader`
    ///
    /// The name for the reader is used only for error messages, e.g. to append
    /// the name of the file being read to the line number on which a problem
    /// occurred.
    ///
    /// ```
    /// use mga2::fastq::FastqReader;
    /// use std::fs::File;
    /// # use std::io;
    /// use std::io::BufReader;
    ///
    /// # fn main() -> io::Result<()> {
    /// let filename = "test.fq";
    /// let file = File::open(filename)?;
    /// let mut reader = FastqReader::with_name(BufReader::new(file), filename);
    /// # Ok(())
    /// # }
    /// ```
    pub fn with_name(reader: R, name: &str) -> Self {
        FastqReader::with_capacity(reader, name, 50, 160)
    }

    /// Create a new `FastqReader` that will allocate new records with the
    /// specified capacities for
    ///
    /// The name for the reader is used only for error messages, e.g. to append
    /// the name of the file being read to the line number on which a problem
    /// occurred.
    ///
    /// ```
    /// use mga2::fastq::FastqReader;
    /// use std::fs::File;
    /// # use std::io;
    /// use std::io::BufReader;
    ///
    /// # fn main() -> io::Result<()> {
    /// let filename = "test.fq";
    /// let file = File::open(filename)?;
    /// let mut reader = FastqReader::with_name(BufReader::new(file), filename);
    /// # Ok(())
    /// # }
    /// ```
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

    /// Reads the next line into a buffer and returns the number of bytes read.
    fn read_next_line(&mut self, buffer: &mut String) -> usize {
        let number_of_bytes = self.reader.read_line(buffer).unwrap();
        self.line_count += 1;
        number_of_bytes
    }

    /// Reads the next FASTQ record.
    ///
    /// This method creates a new `FastqRecord`. This may be desirable if
    /// records are being added to a collection for subsequent processing but if
    /// records are processed one at a time and are no longer required
    /// afterwards, the [`read_next_into`](struct.FastqReader.html#method.read_next_into)
    /// method may be preferable as it will avoid unnecessary and costly memory
    /// allocation.
    ///
    /// # Example
    ///
    /// ```
    /// # use anyhow::Result;
    /// use mga2::fastq::FastqReader;
    /// use std::io::{self, BufReader};
    ///
    /// # fn main() -> Result<()> {
    /// let mut reader = FastqReader::new(BufReader::new(io::stdin()));
    ///
    /// let record = reader.read_next();
    /// if let Some(record) = reader.read_next()? {
    ///     println!("{} {}", record.id, record.seq);
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn read_next(&mut self) -> Result<Option<FastqRecord>> {
        let mut record = FastqRecord::with_capacity(self.id_capacity, self.sequence_capacity);
        if self.read_next_into(&mut record)? {
            Ok(Some(record))
        } else {
            Ok(None)
        }
    }

    /// Reads the next FASTQ record into an existing `FastqRecord`.
    ///
    /// This function is preferred over [`read_next`](struct.FastqReader.html#method.read_next)
    /// when iterating over large numbers of FASTQ records where the record is
    /// only required for the duration of that iteration as it avoids the cost
    /// of allocating memory for a new `FastqRecord` for each record; instead
    /// the provided `FastqRecord` is reused with its contents overwritten.
    ///
    /// # Example
    ///
    /// ```
    /// # use anyhow::Result;
    /// use mga2::fastq::{FastqReader, FastqRecord};
    /// use std::io::{self, BufReader};
    ///
    /// # fn main() -> Result<()> {
    /// let mut reader = FastqReader::new(BufReader::new(io::stdin()));
    /// let mut record = FastqRecord::new();
    ///
    /// let mut number_of_records = 0;
    /// let mut number_of_bases = 0;
    ///
    /// while reader.read_next_into(&mut record)? {
    ///     number_of_records += 1;
    ///     number_of_bases += record.seq.len();
    /// }
    /// # Ok(())
    /// # }
    /// ```
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
                            "Expected '@' character at beginning of line {}, {}",
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
                            "Unexpected problem reading line {}, {}",
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
                "Error extracting sequence id for record at line {}, {}",
                self.line_count,
                self.name
            );
        }

        let mut length = 0;

        loop {
            let number_of_bytes = self.read_next_line(&mut record.seq);
            match number_of_bytes {
                0 => bail!(
                    "Incomplete record at line {}, {}",
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
                            "Sequence contains space(s) for record at line {}, {}",
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
                    "Incomplete record at line {}, {}",
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
                "Sequence and quality lengths differ for record ending at line {}, {}",
                self.line_count,
                self.name
            );
        }

        // check for spaces in quality scores
        if record.qual.find(' ').is_some() {
            bail!(
                "Qualities contains space(s) for record ending at line {}, {}",
                self.line_count,
                self.name
            );
        }

        Ok(true)
    }
}

/// A writer for FASTQ records.
///
/// A `FastqWriter` writes  `FastqRecord` structs to an underlying `BufWriter`,
/// usually to a file. The [`create_fastq_writer`](fn.create_fastq_writer.html)
/// function is the preferred way of creating a `FastqWriter` when writing
/// gzip compressed files.
///
/// # Examples
///
/// ## Writing to stdout
///
/// ```
/// # use anyhow::Result;
/// use mga2::fastq::{FastqWriter, FastqRecord};
/// use std::io::{self, BufWriter};
///
/// # fn main() -> Result<()> {
/// let mut writer = FastqWriter::new(BufWriter::new(io::stdout()));
///
/// let record = FastqRecord {
///     id: "MDE123".to_string(),
///     desc: Some("An example FASTQ record".to_string()),
///     seq: "TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA".to_string(),
///     qual: "AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ".to_string(),
/// };
///
/// writer.write_fastq(&record)?;
///
/// writer.flush()?;
/// # Ok(())
/// # }
/// ```
///
/// ## Writing to a file
///
/// ```
/// # use anyhow::Result;
/// use mga2::fastq::{FastqWriter, FastqRecord};
/// use std::io::BufWriter;
/// use std::fs::File;
///
/// # fn main() -> Result<()> {
/// let file = File::create("test.fq")?;
/// let mut writer = FastqWriter::new(BufWriter::new(file));
///
/// let record = FastqRecord {
///     id: "MDE123".to_string(),
///     desc: Some("An example FASTQ record".to_string()),
///     seq: "TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA".to_string(),
///     qual: "AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ".to_string(),
/// };
///
/// writer.write_fastq(&record)?;
///
/// writer.flush()?;
/// # Ok(())
/// # }
/// ```
///
/// ## Writing to a gzip compressed file
///
/// ```
/// # use anyhow::Result;
/// use mga2::fastq::{FastqRecord, create_fastq_writer};
/// use std::path::PathBuf;
///
/// # fn main() -> Result<()> {
/// let path = PathBuf::from("test.fq.gz");
/// let mut writer = create_fastq_writer(&Some(path))?;
///
/// let record = FastqRecord {
///     id: "MDE123".to_string(),
///     desc: Some("An example FASTQ record".to_string()),
///     seq: "TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA".to_string(),
///     qual: "AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ".to_string(),
/// };
///
/// writer.write_fastq(&record)?;
///
/// writer.flush()?;
/// # Ok(())
/// # }
/// ```
pub struct FastqWriter<W: Write> {
    writer: BufWriter<W>,
}

impl<W: Write> FastqWriter<W> {
    /// Create a new `FastqWriter`.
    ///
    /// ```
    /// use mga2::fastq::FastqWriter;
    /// use std::io::{self, BufWriter};
    ///
    /// let mut writer = FastqWriter::new(BufWriter::new(io::stdout()));
    /// ```
    pub fn new(writer: BufWriter<W>) -> Self {
        FastqWriter { writer }
    }

    /// Write a FASTQ record.
    ///
    /// # Example
    ///
    /// ```
    /// # use anyhow::Result;
    /// use mga2::fastq::{FastqWriter, FastqRecord};
    /// use std::io::{self, BufWriter};
    ///
    /// # fn main() -> Result<()> {
    /// let mut writer = FastqWriter::new(BufWriter::new(io::stdout()));
    ///
    /// let record = FastqRecord {
    ///     id: "MDE123".to_string(),
    ///     desc: Some("An example FASTQ record".to_string()),
    ///     seq: "TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA".to_string(),
    ///     qual: "AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ".to_string(),
    /// };
    ///
    /// writer.write_fastq(&record)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_fastq(&mut self, record: &FastqRecord) -> Result<()> {
        write_fastq(&mut self.writer, record)
            .with_context(|| format!("Error writing FASTQ format for record: {}", record.id))?;
        Ok(())
    }

    /// Write a `FastqRecord` in FASTA format.
    ///
    /// The FASTA format contains a header similar to that for a FASTQ record
    /// but with a '>' character in place of the '@' and with only the sequence,
    /// no quality score sting.
    ///
    /// The sequences will be split across multiple lines each up to 80
    /// characters long.
    ///
    /// # Example
    ///
    /// ```
    /// # use anyhow::Result;
    /// use mga2::fastq::{FastqWriter, FastqRecord};
    /// use std::io::{self, BufWriter};
    ///
    /// # fn main() -> Result<()> {
    /// let mut writer = FastqWriter::new(BufWriter::new(io::stdout()));
    ///
    /// let record = FastqRecord {
    ///     id: "MDE123".to_string(),
    ///     desc: Some("An example FASTQ record".to_string()),
    ///     seq: "TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA".to_string(),
    ///     qual: "AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ".to_string(),
    /// };
    ///
    /// writer.write_fasta(&record)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_fasta(&mut self, record: &FastqRecord) -> Result<()> {
        write_fasta(&mut self.writer, record)
            .with_context(|| format!("Error writing FASTA format for record: {}", record.id))?;
        Ok(())
    }

    /// Flush this output stream, ensuring that all intermediate, buffered
    /// contents are written.
    ///
    /// ```
    /// # use anyhow::Result;
    /// use mga2::fastq::{FastqWriter, FastqRecord};
    /// use std::io::{self, BufWriter};
    ///
    /// # fn main() -> Result<()> {
    /// let mut writer = FastqWriter::new(BufWriter::new(io::stdout()));
    ///
    /// let record = FastqRecord {
    ///     id: "MDE123".to_string(),
    ///     desc: Some("An example FASTQ record".to_string()),
    ///     seq: "TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA".to_string(),
    ///     qual: "AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ".to_string(),
    /// };
    ///
    /// writer.write_fastq(&record)?;
    ///
    /// writer.flush()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()?;
        Ok(())
    }
}

/// Create a `FastqReader` for reading `FastqRecord`s from a file.
///
/// This is the preferred way of creating a `FastqReader` for reading
/// gzip compressed files which are recognized by having a '.gz' extension.
///
/// # Example
///
/// ## Reading from a gzip compressed file
///
/// ```
/// # use anyhow::Result;
/// use mga2::fastq::{FastqRecord, create_fastq_reader};
/// use std::path::PathBuf;
///
/// # fn main() -> Result<()> {
/// let path = PathBuf::from("test.fq.gz");
/// let mut reader = create_fastq_reader(&path)?;
///
/// let record = reader.read_next();
/// if let Some(record) = reader.read_next()? {
///     println!("{} {}", record.id, record.seq);
/// }
/// # Ok(())
/// # }
/// ```
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

/// Create a `FastqWriter` for writing `FastqRecord`s to a file.
///
/// This is the preferred way of creating a `FastqWriter` when writing
/// gzip compressed files which are recognized by having a '.gz' extension.
///
/// # Examples
///
/// ## Writing to a gzip compressed file
///
/// ```
/// # use anyhow::Result;
/// use mga2::fastq::{FastqRecord, create_fastq_writer};
/// use std::path::PathBuf;
///
/// # fn main() -> Result<()> {
/// let path = PathBuf::from("test.fq.gz");
/// let mut writer = create_fastq_writer(&Some(path))?;
///
/// let record = FastqRecord {
///     id: "MDE123".to_string(),
///     desc: Some("An example FASTQ record".to_string()),
///     seq: "TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA".to_string(),
///     qual: "AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ".to_string(),
/// };
///
/// writer.write_fastq(&record)?;
///
/// writer.flush()?;
/// # Ok(())
/// # }
/// ```
///
/// ## Writing to stdout
///
/// ```
/// # use anyhow::Result;
/// use mga2::fastq::{FastqRecord, create_fastq_writer};
///
/// # fn main() -> Result<()> {
/// let mut writer = create_fastq_writer(&None)?;
///
/// let record = FastqRecord {
///     id: "MDE123".to_string(),
///     desc: Some("An example FASTQ record".to_string()),
///     seq: "TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA".to_string(),
///     qual: "AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ".to_string(),
/// };
///
/// writer.write_fastq(&record)?;
///
/// writer.flush()?;
/// # Ok(())
/// # }
/// ```
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
        assert!(error.to_string().starts_with("Incomplete record at line"));
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
            .starts_with("Error extracting sequence id for record at line"));
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
        assert_eq!(
            record.seq,
            "GAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA".to_string()
        );
        assert_eq!(
            record.qual,
            "JJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ".to_string()
        );
    }

    #[test]
    fn trim_beyond_end() {
        let mut record = create_record();
        let result = record.trim(11, Some(1000));
        assert!(result.is_ok(), "Unexpected error during trimming");
        assert!(record.check().is_ok(), "Invalid record after trimming");
        assert_eq!(
            record.seq,
            "GAAGTTGTTAAAATTTCCGGAGGTAGCCATTATATACCAA".to_string()
        );
        assert_eq!(
            record.qual,
            "JJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ".to_string()
        );
    }
    #[test]
    fn trim_to_length() {
        let mut record = create_record();
        let result = record.trim_to_length(36);
        assert!(result.is_ok(), "Unexpected error during trimming");
        assert!(record.check().is_ok(), "Invalid record after trimming");
        assert_eq!(
            record.seq,
            "TGTGACCCAAGAAGTTGTTAAAATTTCCGGAGGTAG".to_string()
        );
        assert_eq!(
            record.qual,
            "AAFFFJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJ".to_string()
        );
    }
}
