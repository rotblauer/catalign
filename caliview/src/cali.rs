//! CALI file format reader
//!
//! Binary format for multi-scale alignment metrics.

use anyhow::{Context, Result};
use byteorder::{LittleEndian, ReadBytesExt};
use memmap2::Mmap;
use std::fs::File;
use std::io::{Cursor, Read};
use std::path::Path;

/// Magic number for CALI files
const CALI_MAGIC: &[u8; 4] = b"CALI";

/// Current format version
const CALI_VERSION: u32 = 1;

/// Default tile sizes
pub const DEFAULT_TILE_SIZES: &[u32] = &[1_000, 10_000, 100_000, 1_000_000];

/// CALI file header
#[derive(Debug, Clone)]
pub struct CaliHeader {
    pub version: u32,
    pub flags: u32,
    pub num_chromosomes: u32,
    pub num_samples: u32,
    pub tile_sizes: Vec<u32>,
    pub reference_name: String,
    pub sample_name: String,
}

impl CaliHeader {
    pub fn from_bytes(data: &[u8]) -> Result<Self> {
        let mut cursor = Cursor::new(data);

        // Check magic
        let mut magic = [0u8; 4];
        cursor.read_exact(&mut magic)?;
        if &magic != CALI_MAGIC {
            anyhow::bail!("Invalid CALI file: bad magic number");
        }

        let version = cursor.read_u32::<LittleEndian>()?;
        let flags = cursor.read_u32::<LittleEndian>()?;
        let num_chromosomes = cursor.read_u32::<LittleEndian>()?;
        let num_samples = cursor.read_u32::<LittleEndian>()?;

        // Read tile sizes
        let mut tile_sizes = Vec::new();
        for _ in 0..8 {
            let ts = cursor.read_u32::<LittleEndian>()?;
            if ts > 0 {
                tile_sizes.push(ts);
            }
        }

        // Read strings
        let mut ref_buf = [0u8; 64];
        cursor.read_exact(&mut ref_buf)?;
        let reference_name = String::from_utf8_lossy(&ref_buf)
            .trim_end_matches('\0')
            .to_string();

        let mut sample_buf = [0u8; 64];
        cursor.read_exact(&mut sample_buf)?;
        let sample_name = String::from_utf8_lossy(&sample_buf)
            .trim_end_matches('\0')
            .to_string();

        Ok(Self {
            version,
            flags,
            num_chromosomes,
            num_samples,
            tile_sizes,
            reference_name,
            sample_name,
        })
    }
}

/// Chromosome information
#[derive(Debug, Clone)]
pub struct ChromosomeInfo {
    pub name: String,
    pub length: u64,
    pub data_offset: u64,
    pub data_length: u64,
    pub tile_count: u32,
}

/// Metrics for a single tile
#[derive(Debug, Clone, Copy)]
pub struct TileMetrics {
    pub start: u32,
    pub end: u32,
    pub identity: f32,
    pub coverage: f32,
    pub energy: f32,
    pub gap_rate: f32,
    pub quality_score: f32,
    pub match_count: u32,
    pub mismatch_count: u32,
    pub insertion_count: u32,
    pub deletion_count: u32,
}

impl TileMetrics {
    pub fn from_bytes(data: &[u8]) -> Result<Self> {
        let mut cursor = Cursor::new(data);

        let start = cursor.read_u32::<LittleEndian>()?;
        let end = cursor.read_u32::<LittleEndian>()?;
        let identity = cursor.read_f32::<LittleEndian>()?;
        let coverage = cursor.read_f32::<LittleEndian>()?;
        let energy = cursor.read_f32::<LittleEndian>()?;
        let gap_rate = cursor.read_f32::<LittleEndian>()?;
        let quality_score = cursor.read_f32::<LittleEndian>()?;

        // Skip reserved fields
        let _ = cursor.read_f32::<LittleEndian>()?;
        let _ = cursor.read_f32::<LittleEndian>()?;
        let _ = cursor.read_f32::<LittleEndian>()?;

        let match_count = cursor.read_u32::<LittleEndian>()?;
        let mismatch_count = cursor.read_u32::<LittleEndian>()?;
        let insertion_count = cursor.read_u32::<LittleEndian>()?;
        let deletion_count = cursor.read_u32::<LittleEndian>()?;

        Ok(Self {
            start,
            end,
            identity,
            coverage,
            energy,
            gap_rate,
            quality_score,
            match_count,
            mismatch_count,
            insertion_count,
            deletion_count,
        })
    }

    pub const BYTE_SIZE: usize = 56;
}

/// CALI file reader with memory-mapped access
pub struct CaliReader {
    _file: File,
    mmap: Mmap,
    header: CaliHeader,
    chromosomes: Vec<ChromosomeInfo>,
}

impl CaliReader {
    /// Open a CALI file
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open {}", path.display()))?;

        let mmap = unsafe { Mmap::map(&file)? };

        // Read header (256 bytes)
        let header = CaliHeader::from_bytes(&mmap[..256])?;

        // Read chromosome table
        let mut chromosomes = Vec::new();
        let mut offset = 256;

        for _ in 0..header.num_chromosomes {
            let chrom = Self::read_chromosome_info(&mmap, offset)?;
            offset += 88; // Chromosome info size
            chromosomes.push(chrom);
        }

        Ok(Self {
            _file: file,
            mmap,
            header,
            chromosomes,
        })
    }

    fn read_chromosome_info(mmap: &Mmap, offset: usize) -> Result<ChromosomeInfo> {
        let data = &mmap[offset..offset + 88];
        let mut cursor = Cursor::new(data);

        let mut name_buf = [0u8; 32];
        cursor.read_exact(&mut name_buf)?;
        let name = String::from_utf8_lossy(&name_buf)
            .trim_end_matches('\0')
            .to_string();

        let length = cursor.read_u64::<LittleEndian>()?;
        let data_offset = cursor.read_u64::<LittleEndian>()?;
        let data_length = cursor.read_u64::<LittleEndian>()?;

        // Read first tile count (simplified)
        let tile_count = cursor.read_u32::<LittleEndian>()?;

        Ok(ChromosomeInfo {
            name,
            length,
            data_offset,
            data_length,
            tile_count,
        })
    }

    /// Get the file header
    pub fn header(&self) -> &CaliHeader {
        &self.header
    }

    /// Get chromosome information
    pub fn chromosomes(&self) -> &[ChromosomeInfo] {
        &self.chromosomes
    }

    /// Get chromosome by name
    pub fn get_chromosome(&self, name: &str) -> Option<&ChromosomeInfo> {
        self.chromosomes.iter().find(|c| c.name == name)
    }

    /// Get tiles for a region
    pub fn get_tiles(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        tile_size: u32,
    ) -> Result<Vec<TileMetrics>> {
        let chrom_info = self.get_chromosome(chrom)
            .ok_or_else(|| anyhow::anyhow!("Unknown chromosome: {}", chrom))?;

        let mut tiles = Vec::new();
        let data_start = chrom_info.data_offset as usize;
        let mut offset = data_start;

        // Find the right tile level
        for &ts in &self.header.tile_sizes {
            let count = self.read_tile_count(offset)?;
            offset += 4;

            if ts == tile_size {
                // Read tiles at this level
                for _ in 0..count {
                    let tile = TileMetrics::from_bytes(
                        &self.mmap[offset..offset + TileMetrics::BYTE_SIZE]
                    )?;

                    if tile.end as u64 >= start && (tile.start as u64) <= end {
                        tiles.push(tile);
                    }

                    offset += TileMetrics::BYTE_SIZE;
                }
                break;
            } else {
                // Skip this level
                offset += count as usize * TileMetrics::BYTE_SIZE;
            }
        }

        Ok(tiles)
    }

    fn read_tile_count(&self, offset: usize) -> Result<u32> {
        let mut cursor = Cursor::new(&self.mmap[offset..offset + 4]);
        Ok(cursor.read_u32::<LittleEndian>()?)
    }

    /// Get summary for a region (auto-select resolution)
    pub fn get_region_summary(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
    ) -> Result<RegionSummary> {
        let region_size = end - start;

        // Auto-select tile size (aim for ~100 tiles)
        let tile_size = self.header.tile_sizes
            .iter()
            .rev()
            .find(|&&ts| region_size / ts as u64 >= 10)
            .copied()
            .unwrap_or(self.header.tile_sizes[0]);

        let tiles = self.get_tiles(chrom, start, end, tile_size)?;

        if tiles.is_empty() {
            return Ok(RegionSummary::default());
        }

        let identities: Vec<f32> = tiles.iter()
            .filter(|t| t.coverage > 0.0)
            .map(|t| t.identity)
            .collect();

        let coverages: Vec<f32> = tiles.iter()
            .map(|t| t.coverage)
            .collect();

        Ok(RegionSummary {
            tile_size,
            num_tiles: tiles.len() as u32,
            mean_identity: if identities.is_empty() { 0.0 } else {
                identities.iter().sum::<f32>() / identities.len() as f32
            },
            mean_coverage: if coverages.is_empty() { 0.0 } else {
                coverages.iter().sum::<f32>() / coverages.len() as f32
            },
            min_identity: identities.iter().cloned().fold(f32::INFINITY, f32::min),
            max_identity: identities.iter().cloned().fold(f32::NEG_INFINITY, f32::max),
        })
    }
}

/// Summary statistics for a region
#[derive(Debug, Default)]
pub struct RegionSummary {
    pub tile_size: u32,
    pub num_tiles: u32,
    pub mean_identity: f32,
    pub mean_coverage: f32,
    pub min_identity: f32,
    pub max_identity: f32,
}
