use anyhow::{bail, Context, Result};
use flate2::bufread::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{stdout, BufRead, BufReader, BufWriter, ErrorKind, Read, Write};
use std::path::PathBuf;

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

    pub fn trim(&mut self, start: usize, end: Option<usize>) -> Result<()> {
        if self.seq.len() != self.qual.len() {
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

        if start > self.seq.len() {
            self.seq = String::new();
            self.qual = String::new();
        } else {
            match end {
                Some(end) => {
                    self.seq = self.seq.as_str()[(start - 1)..end].to_string();
                    self.qual = self.qual.as_str()[(start - 1)..end].to_string();
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
    id_capacity: usize,
    sequence_capacity: usize,
    reader: R,
    test_separator_char: [u8; 1],
    line_count: u32,
}

impl<R: BufRead> FastqReader<R> {
    pub fn new(reader: R) -> Self {
        FastqReader::with_capacity(reader, 50, 160)
    }

    pub fn with_capacity(reader: R, id_capacity: usize, sequence_capacity: usize) -> Self {
        FastqReader {
            id_capacity,
            sequence_capacity,
            reader,
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

        match self.reader.read_exact(&mut self.test_separator_char) {
            Ok(_) => {
                if self.test_separator_char[0] != b'@' {
                    bail!(
                        "expected '@' character at beginning of line {}",
                        self.line_count + 1
                    );
                }
            }
            Err(error) => {
                if error.kind() == ErrorKind::UnexpectedEof {
                    return Ok(false);
                } else {
                    bail!("unexpected end of file at line {}", self.line_count + 1);
                }
            }
        }

        if self.read_next_line(&mut record.id) == 0 {
            return Ok(false);
        }

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
                "error extracting sequence id for record at line {}",
                self.line_count
            );
        }

        let mut length = 0;

        loop {
            let number_of_bytes = self.read_next_line(&mut record.seq);
            match number_of_bytes {
                0 => bail!("incomplete record at line {}", self.line_count - 1),
                1 => bail!("empty line at line number {}", self.line_count),
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
                            "sequence contains space(s) for record at line {}",
                            self.line_count
                        );
                    }
                    length = record.seq.len();
                }
            };
        }

        if record.seq.is_empty() {
            bail!("zero length sequence at line {}", self.line_count);
        }

        let mut qual_length = 0;
        while qual_length < length {
            let number_of_bytes = self.read_next_line(&mut record.qual);
            match number_of_bytes {
                0 => bail!("incomplete record at line {}", self.line_count - 1),
                1 => bail!("empty line at line number {}", self.line_count),
                _ => record.qual.pop(),
            };
            // check for spaces in quality scores
            if record.qual[qual_length..].find(' ').is_some() {
                bail!(
                    "qualities contains space(s) for record at line {}",
                    self.line_count
                );
            }
            qual_length = record.qual.len();
        }

        if record.qual.len() != length {
            bail!(
                "sequence and quality lengths differ for record ending at line {}",
                self.line_count
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

pub fn create_fastq_reader(fastq_file: &PathBuf) -> Result<FastqReader<BufReader<Box<dyn Read>>>> {
    let file =
        File::open(fastq_file).with_context(|| format!("Error opening file {:?}", fastq_file))?;

    let reader: Box<dyn Read> = if fastq_file.to_str().unwrap().ends_with(".gz") {
        Box::new(MultiGzDecoder::new(BufReader::with_capacity(
            64 * 1024,
            file,
        )))
    } else {
        Box::new(file)
    };

    let buffered_reader = BufReader::with_capacity(64 * 1024, reader);
    let fastq_reader = FastqReader::new(buffered_reader);

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
