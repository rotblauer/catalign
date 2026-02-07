//! Caliview - High-performance genome alignment viewer
//!
//! A fast, GPU-accelerated viewer for Catalign multi-scale alignment metrics.

mod app;
mod cali;
mod tracks;
mod ui;

use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "caliview")]
#[command(about = "High-performance genome alignment viewer", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// View a CALI file interactively
    View {
        /// Path to CALI file
        #[arg(value_name = "FILE")]
        file: PathBuf,

        /// Initial chromosome to display
        #[arg(short, long)]
        chrom: Option<String>,

        /// Initial position (e.g., "chr1:1000000-2000000")
        #[arg(short, long)]
        position: Option<String>,
    },

    /// Convert BAM/CRAM to CALI format
    Convert {
        /// Input BAM/CRAM file
        #[arg(value_name = "INPUT")]
        input: PathBuf,

        /// Output CALI file
        #[arg(short, long, value_name = "OUTPUT")]
        output: PathBuf,

        /// Reference FASTA (required for CRAM)
        #[arg(short, long)]
        reference: Option<PathBuf>,

        /// Minimum mapping quality
        #[arg(long, default_value = "0")]
        min_mapq: u8,
    },

    /// Export region to various formats
    Export {
        /// CALI file
        #[arg(value_name = "FILE")]
        file: PathBuf,

        /// Region (e.g., "chr1:1000000-2000000")
        #[arg(short, long)]
        region: String,

        /// Output file
        #[arg(short, long)]
        output: PathBuf,

        /// Output format (json, csv, bed)
        #[arg(short, long, default_value = "json")]
        format: String,
    },

    /// Show file statistics
    Stats {
        /// CALI file
        #[arg(value_name = "FILE")]
        file: PathBuf,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::View { file, chrom, position } => {
            run_viewer(file, chrom, position)
        }
        Commands::Convert { input, output, reference, min_mapq } => {
            convert_to_cali(input, output, reference, min_mapq)
        }
        Commands::Export { file, region, output, format } => {
            export_region(file, region, output, format)
        }
        Commands::Stats { file } => {
            show_stats(file)
        }
    }
}

fn run_viewer(file: PathBuf, chrom: Option<String>, position: Option<String>) -> Result<()> {
    println!("🧬 Caliview - Genome Alignment Viewer");
    println!("Loading: {}", file.display());

    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1400.0, 900.0])
            .with_min_inner_size([800.0, 600.0])
            .with_title("Caliview - Catalign Genome Viewer"),
        ..Default::default()
    };

    eframe::run_native(
        "Caliview",
        native_options,
        Box::new(|cc| Ok(Box::new(app::CaliviewApp::new(cc, file, chrom, position)))),
    ).map_err(|e| anyhow::anyhow!("Failed to start viewer: {}", e))
}

fn convert_to_cali(
    input: PathBuf,
    output: PathBuf,
    reference: Option<PathBuf>,
    min_mapq: u8,
) -> Result<()> {
    println!("Converting {} to CALI format...", input.display());

    // TODO: Implement conversion using noodles
    println!("Conversion complete: {}", output.display());
    Ok(())
}

fn export_region(
    file: PathBuf,
    region: String,
    output: PathBuf,
    format: String,
) -> Result<()> {
    println!("Exporting region {} from {}", region, file.display());

    // TODO: Implement export
    println!("Exported to: {}", output.display());
    Ok(())
}

fn show_stats(file: PathBuf) -> Result<()> {
    println!("📊 CALI File Statistics");
    println!("========================");
    println!("File: {}", file.display());

    let cali = cali::CaliReader::open(&file)?;

    println!("Version: {}", cali.header().version);
    println!("Reference: {}", cali.header().reference_name);
    println!("Sample: {}", cali.header().sample_name);
    println!("Chromosomes: {}", cali.chromosomes().len());
    println!();

    println!("{:<12} {:>15} {:>10}", "Chromosome", "Length", "Tiles");
    println!("{}", "-".repeat(40));

    for chrom in cali.chromosomes() {
        println!(
            "{:<12} {:>15} {:>10}",
            chrom.name,
            format_bp(chrom.length),
            chrom.tile_count
        );
    }

    Ok(())
}

fn format_bp(bp: u64) -> String {
    if bp >= 1_000_000_000 {
        format!("{:.2} Gb", bp as f64 / 1_000_000_000.0)
    } else if bp >= 1_000_000 {
        format!("{:.2} Mb", bp as f64 / 1_000_000.0)
    } else if bp >= 1_000 {
        format!("{:.2} kb", bp as f64 / 1_000.0)
    } else {
        format!("{} bp", bp)
    }
}
