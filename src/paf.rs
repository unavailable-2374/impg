//! PAF (Pairwise Alignment Format) parsing
//!
//! This module provides functions for parsing PAF format alignment files.
//! Supports both uncompressed and BGZF-compressed files with optional GZI indices.

use crate::alignment_record::{AlignmentRecord, Strand};
use crate::seqidx::SequenceIndex;
use log::debug;
use noodles::bgzf;
use std::cell::RefCell;
use std::fs::File;
use std::io::{BufRead, BufReader, Error as IoError, Read, Seek, SeekFrom};
use std::num::{NonZeroUsize, ParseIntError};

#[derive(Debug)]
pub enum ParseErr {
    NotEnoughFields,
    IoError(IoError),
    InvalidField(ParseIntError),
    InvalidStrand,
    InvalidCigarFormat,
    UnsupportedCigarOperation,
    InvalidFormat(String),
}

impl std::fmt::Display for ParseErr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseErr::NotEnoughFields => write!(f, "Not enough fields in PAF record"),
            ParseErr::IoError(e) => write!(f, "IO error: {}", e),
            ParseErr::InvalidField(e) => write!(f, "Invalid field: {}", e),
            ParseErr::InvalidStrand => write!(f, "Invalid strand"),
            ParseErr::InvalidCigarFormat => write!(f, "Invalid CIGAR format"),
            ParseErr::UnsupportedCigarOperation => write!(f, "Unsupported CIGAR operation"),
            ParseErr::InvalidFormat(msg) => write!(f, "{}", msg),
        }
    }
}

impl std::error::Error for ParseErr {}

enum PafHandle {
    Plain(File),
    Compressed(bgzf::io::Reader<File>),
}

const BGZF_HEADER_SIZE: usize = 18;

/// Check whether a file starts with a valid BGZF header.
/// Returns `Ok(false)` for regular gzip, too-small files, or plain text.
fn is_bgzf<R: Read + Seek>(reader: &mut R) -> std::io::Result<bool> {
    let mut header = [0u8; BGZF_HEADER_SIZE];
    let result = match reader.read_exact(&mut header) {
        Ok(()) => {
            Ok(header[0..2] == [0x1f, 0x8b]      // gzip magic
                && header[2] == 0x08              // DEFLATE
                && header[3] == 0x04              // FEXTRA
                && header[10..12] == [0x06, 0x00] // XLEN=6
                && header[12..14] == [b'B', b'C'] // BC subfield
                && header[14..16] == [0x02, 0x00]) // SLEN=2
        }
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(false),
        Err(e) => Err(e),
    };
    reader.seek(SeekFrom::Start(0))?;
    result
}

/// Open an alignment file for CIGAR random access, validating BGZF framing once.
///
/// This is the *only* place the open + `is_bgzf` header probe happens; callers
/// go through [`read_cigar_data`], which caches the returned handle per thread so
/// repeated reads from the same file reuse one open fd + decoder.
fn open_cigar_handle(alignment_file: &str) -> Result<PafHandle, String> {
    let is_compressed = [".gz", ".bgz"]
        .iter()
        .any(|extension| alignment_file.ends_with(extension));

    if is_compressed {
        let mut file = File::open(alignment_file)
            .map_err(|e| format!("Failed to open compressed file '{}': {}", alignment_file, e))?;
        if !is_bgzf(&mut file)
            .map_err(|e| format!("Failed to read header of '{}': {}", alignment_file, e))?
        {
            return Err(format!(
                "'{}' is regular gzip, not BGZF. Convert with: zcat '{}' | bgzip > output.paf.gz",
                alignment_file, alignment_file
            ));
        }
        Ok(PafHandle::Compressed(bgzf::io::Reader::new(file)))
    } else {
        let file = File::open(alignment_file)
            .map_err(|e| format!("Failed to open file '{}': {}", alignment_file, e))?;
        Ok(PafHandle::Plain(file))
    }
}

/// Seek to `offset` (a BGZF virtual position for compressed files, a plain byte
/// offset otherwise) and fill `buffer`. Every read seeks first, so reusing a
/// handle across unrelated reads is safe.
fn read_cigar_from_handle(
    handle: &mut PafHandle,
    alignment_file: &str,
    offset: u64,
    buffer: &mut [u8],
) -> Result<(), String> {
    match handle {
        PafHandle::Compressed(reader) => {
            let virtual_position = bgzf::VirtualPosition::from(offset);
            reader.seek(virtual_position).map_err(|e| {
                format!(
                    "Failed to seek in compressed file '{}': {}",
                    alignment_file, e
                )
            })?;
            reader.read_exact(buffer).map_err(|e| {
                format!(
                    "Failed to read data from compressed file '{}': {}",
                    alignment_file, e
                )
            })
        }
        PafHandle::Plain(file) => {
            file.seek(SeekFrom::Start(offset))
                .map_err(|e| format!("Failed to seek in file '{}': {}", alignment_file, e))?;
            file.read_exact(buffer)
                .map_err(|e| format!("Failed to read data from file '{}': {}", alignment_file, e))
        }
    }
}

/// Thread-local LRU of open CIGAR readers keyed by alignment-file path.
///
/// The depth `--cigar-precise` path reads one CIGAR per overlapping alignment via
/// [`read_cigar_data`]; the previous implementation did a fresh `File::open` +
/// `is_bgzf` probe + decoder construction on *every* read. At 10^5 per-file
/// indices on a parallel filesystem (Lustre) with 112 threads, that `open()`
/// metadata storm — not the O(1) virtual seek — dominated wall-clock. Reusing an
/// already-open handle collapses the K reopens of one file within a single query
/// (and repeat queries of the same file on the same worker) to at most one open.
struct CigarReaderCache {
    /// Front = most-recently-used. Small `cap`, so linear scan is fine.
    entries: Vec<(String, PafHandle)>,
    cap: usize,
}

impl CigarReaderCache {
    fn new() -> Self {
        // Per-thread open-fd budget. cap=1 already captures the dominant win
        // (all CIGAR reads within one sub-index query hit the same file); a
        // larger cap adds cross-query reuse. Bounded so cap × threads stays well
        // under typical `ulimit -n`. Override with IMPG_CIGAR_READER_CACHE.
        let cap = std::env::var("IMPG_CIGAR_READER_CACHE")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .filter(|&n| n > 0)
            .unwrap_or(16);
        Self {
            entries: Vec::with_capacity(cap),
            cap,
        }
    }

    fn read(&mut self, alignment_file: &str, offset: u64, buffer: &mut [u8]) -> Result<(), String> {
        if let Some(pos) = self.entries.iter().position(|(p, _)| p == alignment_file) {
            if pos != 0 {
                let entry = self.entries.remove(pos);
                self.entries.insert(0, entry);
            }
        } else {
            let handle = open_cigar_handle(alignment_file)?;
            if self.entries.len() >= self.cap {
                self.entries.pop(); // evict least-recently-used
            }
            self.entries.insert(0, (alignment_file.to_string(), handle));
        }

        let result = read_cigar_from_handle(&mut self.entries[0].1, alignment_file, offset, buffer);
        if result.is_err() {
            // Drop the (possibly desynced) reader so the next read re-opens fresh.
            self.entries.remove(0);
        }
        result
    }
}

thread_local! {
    static CIGAR_READER_CACHE: RefCell<CigarReaderCache> = RefCell::new(CigarReaderCache::new());
}

pub fn read_cigar_data(alignment_file: &str, offset: u64, buffer: &mut [u8]) -> Result<(), String> {
    CIGAR_READER_CACHE.with(|cache| cache.borrow_mut().read(alignment_file, offset, buffer))
}

/// Parse a single PAF line into an AlignmentRecord
/// This is PAF-format specific parsing logic
fn parse_paf_line(
    line: &str,
    file_pos: u64,
    seq_index: &mut SequenceIndex,
) -> Result<AlignmentRecord, ParseErr> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 12 {
        return Err(ParseErr::NotEnoughFields);
    }

    let query_name = fields[0].to_string();
    let query_length = fields[1].parse::<usize>().map_err(ParseErr::InvalidField)?;
    let query_start = fields[2].parse::<usize>().map_err(ParseErr::InvalidField)?;
    let query_end = fields[3].parse::<usize>().map_err(ParseErr::InvalidField)?;
    let target_name = fields[5].to_string();
    let target_length = fields[6].parse::<usize>().map_err(ParseErr::InvalidField)?;
    let target_start = fields[7].parse::<usize>().map_err(ParseErr::InvalidField)?;
    let target_end = fields[8].parse::<usize>().map_err(ParseErr::InvalidField)?;
    let strand_char = fields[4]
        .chars()
        .next()
        .ok_or_else(|| ParseErr::InvalidFormat("Expected '+' or '-' for strand".to_string()))?;
    let strand = match strand_char {
        '+' => Strand::Forward,
        '-' => Strand::Reverse,
        _ => return Err(ParseErr::InvalidStrand),
    };

    // Convert names to IDs using the SequenceIndex
    let query_id = seq_index.get_or_insert_id(&query_name, Some(query_length));
    let target_id = seq_index.get_or_insert_id(&target_name, Some(target_length));

    let mut cigar_offset: u64 = file_pos;
    let mut cigar_bytes: usize = 0;

    for tag_str in fields.iter() {
        if tag_str.starts_with("cg:Z:") {
            cigar_offset += 5;
            cigar_bytes = tag_str.len() - 5;
            break;
        } else {
            cigar_offset += (tag_str.len() + 1) as u64;
        }
    }

    // Create the record and set strand
    let mut record = AlignmentRecord {
        query_id,
        query_start,
        query_end,
        target_id,
        target_start,
        target_end,
        strand_and_data_offset: cigar_offset,
        data_bytes: cigar_bytes,
    };
    record.set_strand(strand);

    Ok(record)
}

pub fn parse_paf<R: BufRead>(
    reader: R,
    seq_index: &mut SequenceIndex,
) -> Result<Vec<AlignmentRecord>, ParseErr> {
    let mut bytes_read: u64 = 0;
    let mut records = Vec::new();
    for line_result in reader.lines() {
        let line = line_result.map_err(ParseErr::IoError)?;
        let record = parse_paf_line(&line, bytes_read, seq_index)?;
        records.push(record);

        // Size of line plus newline
        bytes_read += (line.len() + 1) as u64;
    }
    Ok(records)
}

/// Parse PAF from a BGZF-compressed file, storing virtual positions for seeking.
/// If a GZI index is provided, uses it for faster multithreaded decompression and converts
/// uncompressed offsets to virtual positions. Otherwise, reads with single-threaded BGZF reader.
pub fn parse_paf_bgzf<R: std::io::Read + std::io::Seek>(
    mut reader: noodles::bgzf::io::Reader<R>,
    seq_index: &mut SequenceIndex,
) -> Result<Vec<AlignmentRecord>, ParseErr> {
    use std::io::BufRead;

    let mut records = Vec::new();
    let mut line_bytes = Vec::new();

    loop {
        // Get virtual position BEFORE reading
        let line_start_vpos = reader.virtual_position();
        line_bytes.clear();

        let bytes_read = reader
            .read_until(b'\n', &mut line_bytes)
            .map_err(ParseErr::IoError)?;
        if bytes_read == 0 {
            break; // EOF
        }

        // Convert to string for parsing (excluding newline)
        let line_len = if line_bytes.ends_with(b"\n") {
            line_bytes.len() - 1
        } else {
            line_bytes.len()
        };
        let line = std::str::from_utf8(&line_bytes[..line_len])
            .map_err(|_| ParseErr::InvalidFormat("Invalid UTF-8".to_string()))?;

        if line.is_empty() {
            continue;
        }

        // Parse to get byte offset to CIGAR within the line
        let mut record = parse_paf_line(line, 0, seq_index)?;
        let cigar_byte_offset = record.strand_and_data_offset & !AlignmentRecord::STRAND_BIT;

        // Compute CIGAR virtual position by seeking back and advancing
        // This correctly handles BGZF block boundaries
        reader.seek(line_start_vpos).map_err(ParseErr::IoError)?;

        if cigar_byte_offset > 0 {
            std::io::copy(
                &mut reader.by_ref().take(cigar_byte_offset),
                &mut std::io::sink(),
            )
            .map_err(ParseErr::IoError)?;
        }

        let cigar_vpos = reader.virtual_position();

        // Update record with correct BGZF virtual position
        let strand_bit = record.strand_and_data_offset & AlignmentRecord::STRAND_BIT;
        record.strand_and_data_offset = u64::from(cigar_vpos) | strand_bit;

        records.push(record);

        // Skip remaining bytes to end of line instead of seeking
        // This eliminates one seek operation per line
        let remaining_bytes = line_bytes.len() as u64 - cigar_byte_offset;
        if remaining_bytes > 0 {
            std::io::copy(
                &mut reader.by_ref().take(remaining_bytes),
                &mut std::io::sink(),
            )
            .map_err(ParseErr::IoError)?;
        }
    }

    Ok(records)
}

/// Parse PAF from a BGZF-compressed file using a GZI index for faster multithreaded decompression.
/// After parsing with uncompressed offsets, converts them to virtual positions for seeking.
pub fn parse_paf_bgzf_with_gzi<R: std::io::Read>(
    reader: R,
    gzi_index: noodles::bgzf::gzi::Index,
    seq_index: &mut SequenceIndex,
) -> Result<Vec<AlignmentRecord>, ParseErr> {
    // First pass: parse with uncompressed byte offsets
    let reader = std::io::BufReader::new(reader);
    let mut records = parse_paf(reader, seq_index)?;

    // Second pass: convert uncompressed offsets to virtual positions using GZI
    for record in &mut records {
        // Extract the uncompressed offset (ignoring the strand bit)
        let uncompressed_offset = record.strand_and_data_offset & !AlignmentRecord::STRAND_BIT;

        // Convert to virtual position using GZI query
        let virtual_pos = gzi_index.query(uncompressed_offset).map_err(|e| {
            ParseErr::InvalidFormat(format!(
                "Failed to find virtual position for offset {}: {:?}",
                uncompressed_offset, e
            ))
        })?;

        // Update the record with virtual position, preserving strand bit
        let strand_bit = record.strand_and_data_offset & AlignmentRecord::STRAND_BIT;
        record.strand_and_data_offset = u64::from(virtual_pos) | strand_bit;
    }

    Ok(records)
}

/// Parse a PAF file with automatic format detection (compressed or uncompressed)
/// Handles .gz/.bgz compression with optional GZI index for multithreaded decompression
pub fn parse_paf_file(
    paf_file: &str,
    file: File,
    threads: NonZeroUsize,
    seq_index: &mut SequenceIndex,
) -> std::io::Result<Vec<AlignmentRecord>> {
    if [".gz", ".bgz"].iter().any(|e| paf_file.ends_with(e)) {
        let mut file = file;
        if !is_bgzf(&mut file)? {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "'{}' is regular gzip, not BGZF. Convert with: zcat '{}' | bgzip > output.paf.gz",
                    paf_file, paf_file
                ),
            ));
        }
        let gzi_path = format!("{}.gzi", paf_file);
        if std::path::Path::new(&gzi_path).exists() {
            debug!(
                "Found GZI index for {}, using multithreaded decompression",
                paf_file
            );
            let gzi_index = noodles::bgzf::gzi::fs::read(&gzi_path).map_err(|e| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("Failed to read GZI index {}: {}", gzi_path, e),
                )
            })?;
            let mt_reader =
                noodles::bgzf::io::MultithreadedReader::with_worker_count(threads, file);
            parse_paf_bgzf_with_gzi(mt_reader, gzi_index, seq_index).map_err(|e| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("Failed to parse PAF from {}: {}", paf_file, e),
                )
            })
        } else {
            debug!("No GZI index for {}, using BGZF reader", paf_file);
            let bgzf_reader = noodles::bgzf::io::Reader::new(file);
            parse_paf_bgzf(bgzf_reader, seq_index).map_err(|e| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("Failed to parse PAF from {}: {}", paf_file, e),
                )
            })
        }
    } else {
        let reader = BufReader::new(file);
        parse_paf(reader, seq_index).map_err(|e| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Failed to parse PAF from {}: {}", paf_file, e),
            )
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_paf_valid() {
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\t0\t100\t60\t100\t255";
        let mut seq_index = SequenceIndex::new();
        let record = parse_paf_line(line, 0, &mut seq_index).unwrap();

        // IDs should be 0 and 1 as they're the first entries in the SequenceIndex
        let query_id = seq_index.get_id("seq1").unwrap();
        let target_id = seq_index.get_id("seq2").unwrap();

        assert_eq!(
            record,
            AlignmentRecord {
                query_id,
                query_start: 0,
                query_end: 100,
                target_id,
                target_start: 0,
                target_end: 100,
                // If no cigar, offset is line length; data_bytes=0
                strand_and_data_offset: (line.len() + 1) as u64,
                data_bytes: 0,
            }
        );
    }

    #[test]
    fn test_parse_paf_valid_2() {
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\t0\t100\t60\t100\t255\tcg:Z:10=";
        let mut seq_index = SequenceIndex::new();
        assert!(parse_paf_line(line, 0, &mut seq_index).is_ok());
    }

    #[test]
    fn test_parse_paf_invalid() {
        // it's got a character 'z' in the length field
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\tz\t100\t60\t100\t255\tcg:Z:10M";
        let mut seq_index = SequenceIndex::new();
        assert!(parse_paf_line(line, 0, &mut seq_index).is_err());
    }

    #[test]
    fn test_parse_paf_cigar_invalid() {
        // it's got Q in the CIGAR string
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\tz\t100\t60\t100\t255\tcg:Z:10Q";
        let mut seq_index = SequenceIndex::new();
        assert!(parse_paf_line(line, 0, &mut seq_index).is_err());
    }
}
