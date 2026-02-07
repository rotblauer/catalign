//! Main application state and logic

use eframe::egui;
use std::path::PathBuf;
use std::sync::Arc;

use crate::cali::{CaliReader, TileMetrics};
use crate::tracks::TrackManager;
use crate::ui;

/// Main application state
pub struct CaliviewApp {
    /// CALI file reader
    cali: Option<CaliReader>,

    /// Track manager
    tracks: TrackManager,

    /// Current chromosome
    current_chrom: Option<String>,

    /// Current view range (start, end)
    view_range: (u64, u64),

    /// Cached tiles for current view
    cached_tiles: Vec<TileMetrics>,

    /// Current tile resolution
    current_resolution: u32,

    /// UI state
    show_settings: bool,
    show_stats: bool,

    /// Navigation
    goto_text: String,

    /// Error message
    error_message: Option<String>,
}

impl CaliviewApp {
    pub fn new(
        cc: &eframe::CreationContext<'_>,
        file: PathBuf,
        initial_chrom: Option<String>,
        initial_position: Option<String>,
    ) -> Self {
        // Load CALI file
        let cali = match CaliReader::open(&file) {
            Ok(c) => Some(c),
            Err(e) => {
                eprintln!("Error loading file: {}", e);
                None
            }
        };

        // Determine initial view
        let (current_chrom, view_range) = if let Some(ref c) = cali {
            let chrom = initial_chrom
                .or_else(|| c.chromosomes().first().map(|ch| ch.name.clone()));

            let range = if let Some(ref pos) = initial_position {
                parse_position(pos).unwrap_or((0, 1_000_000))
            } else if let Some(ref ch) = chrom {
                if let Some(info) = c.get_chromosome(ch) {
                    (0, info.length.min(1_000_000))
                } else {
                    (0, 1_000_000)
                }
            } else {
                (0, 1_000_000)
            };

            (chrom, range)
        } else {
            (None, (0, 1_000_000))
        };

        let mut app = Self {
            cali,
            tracks: TrackManager::new(),
            current_chrom,
            view_range,
            cached_tiles: Vec::new(),
            current_resolution: 10_000,
            show_settings: false,
            show_stats: false,
            goto_text: String::new(),
            error_message: None,
        };

        // Initial tile load
        app.refresh_tiles();

        app
    }

    fn refresh_tiles(&mut self) {
        if let (Some(cali), Some(chrom)) = (&self.cali, &self.current_chrom) {
            let region_size = self.view_range.1 - self.view_range.0;

            // Auto-select resolution
            let resolution = if region_size > 10_000_000 {
                1_000_000
            } else if region_size > 1_000_000 {
                100_000
            } else if region_size > 100_000 {
                10_000
            } else {
                1_000
            };

            self.current_resolution = resolution;

            match cali.get_tiles(chrom, self.view_range.0, self.view_range.1, resolution) {
                Ok(tiles) => {
                    self.cached_tiles = tiles;
                    self.error_message = None;
                }
                Err(e) => {
                    self.error_message = Some(format!("Error loading tiles: {}", e));
                }
            }
        }
    }

    fn zoom_in(&mut self) {
        let center = (self.view_range.0 + self.view_range.1) / 2;
        let half_width = (self.view_range.1 - self.view_range.0) / 4;
        self.view_range = (
            center.saturating_sub(half_width),
            center + half_width,
        );
        self.refresh_tiles();
    }

    fn zoom_out(&mut self) {
        let center = (self.view_range.0 + self.view_range.1) / 2;
        let half_width = self.view_range.1 - self.view_range.0;

        // Get chromosome length
        let max_len = self.cali.as_ref()
            .and_then(|c| self.current_chrom.as_ref().and_then(|ch| c.get_chromosome(ch)))
            .map(|info| info.length)
            .unwrap_or(u64::MAX);

        self.view_range = (
            center.saturating_sub(half_width),
            (center + half_width).min(max_len),
        );
        self.refresh_tiles();
    }

    fn pan_left(&mut self) {
        let width = self.view_range.1 - self.view_range.0;
        let delta = width / 4;
        self.view_range = (
            self.view_range.0.saturating_sub(delta),
            self.view_range.1.saturating_sub(delta),
        );
        self.refresh_tiles();
    }

    fn pan_right(&mut self) {
        let width = self.view_range.1 - self.view_range.0;
        let delta = width / 4;

        let max_len = self.cali.as_ref()
            .and_then(|c| self.current_chrom.as_ref().and_then(|ch| c.get_chromosome(ch)))
            .map(|info| info.length)
            .unwrap_or(u64::MAX);

        let new_end = (self.view_range.1 + delta).min(max_len);
        let new_start = new_end.saturating_sub(width);

        self.view_range = (new_start, new_end);
        self.refresh_tiles();
    }

    fn goto_position(&mut self) {
        if let Some((start, end)) = parse_position(&self.goto_text) {
            self.view_range = (start, end);
            self.refresh_tiles();
        }
    }
}

impl eframe::App for CaliviewApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Top panel - toolbar
        egui::TopBottomPanel::top("toolbar").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.heading("🧬 Caliview");
                ui.separator();

                // Chromosome selector - collect data first to avoid borrow issues
                let chrom_names: Vec<(String, u64)> = self.cali.as_ref()
                    .map(|c| c.chromosomes().iter().map(|ch| (ch.name.clone(), ch.length)).collect())
                    .unwrap_or_default();

                let mut selected_chrom: Option<(String, u64)> = None;

                if !chrom_names.is_empty() {
                    egui::ComboBox::from_label("Chromosome")
                        .selected_text(self.current_chrom.as_deref().unwrap_or("Select..."))
                        .show_ui(ui, |ui| {
                            for (name, length) in &chrom_names {
                                let is_selected = self.current_chrom.as_ref() == Some(name);
                                if ui.selectable_label(is_selected, name).clicked() {
                                    selected_chrom = Some((name.clone(), *length));
                                }
                            }
                        });
                }

                // Apply selection after combo box closes
                if let Some((name, length)) = selected_chrom {
                    self.current_chrom = Some(name);
                    self.view_range = (0, length.min(1_000_000));
                    self.refresh_tiles();
                }

                ui.separator();

                // Navigation
                if ui.button("⬅").clicked() { self.pan_left(); }
                if ui.button("➡").clicked() { self.pan_right(); }
                if ui.button("🔍+").clicked() { self.zoom_in(); }
                if ui.button("🔍-").clicked() { self.zoom_out(); }

                ui.separator();

                // Goto
                ui.label("Go to:");
                let response = ui.text_edit_singleline(&mut self.goto_text);
                if response.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                    self.goto_position();
                }

                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.button("⚙").clicked() {
                        self.show_settings = !self.show_settings;
                    }
                    if ui.button("📊").clicked() {
                        self.show_stats = !self.show_stats;
                    }
                });
            });
        });

        // Left panel - chromosome overview
        egui::SidePanel::left("overview").default_width(150.0).show(ctx, |ui| {
            ui.heading("Chromosomes");
            ui.separator();

            // Collect chromosome info first
            let chrom_list: Vec<(String, u64)> = self.cali.as_ref()
                .map(|c| c.chromosomes().iter().map(|ch| (ch.name.clone(), ch.length)).collect())
                .unwrap_or_default();

            let mut clicked_chrom: Option<(String, u64)> = None;

            egui::ScrollArea::vertical().show(ui, |ui| {
                for (name, length) in &chrom_list {
                    let selected = self.current_chrom.as_ref() == Some(name);
                    if ui.selectable_label(selected, format!("{} ({})", name, format_bp(*length))).clicked() {
                        clicked_chrom = Some((name.clone(), *length));
                    }
                }
            });

            // Apply click after scroll area
            if let Some((name, length)) = clicked_chrom {
                self.current_chrom = Some(name);
                self.view_range = (0, length.min(1_000_000));
                self.refresh_tiles();
            }
        });

        // Right panel - statistics (optional)
        if self.show_stats {
            egui::SidePanel::right("stats").default_width(250.0).show(ctx, |ui| {
                ui.heading("Region Statistics");
                ui.separator();

                if let (Some(cali), Some(chrom)) = (&self.cali, &self.current_chrom) {
                    ui.label(format!("Chromosome: {}", chrom));
                    ui.label(format!("Range: {}-{}",
                        format_bp(self.view_range.0),
                        format_bp(self.view_range.1)
                    ));
                    ui.label(format!("Resolution: {} tiles", self.current_resolution));
                    ui.label(format!("Tiles loaded: {}", self.cached_tiles.len()));

                    ui.separator();

                    if !self.cached_tiles.is_empty() {
                        let avg_identity: f32 = self.cached_tiles.iter()
                            .filter(|t| t.coverage > 0.0)
                            .map(|t| t.identity)
                            .sum::<f32>() / self.cached_tiles.len() as f32;

                        let avg_coverage: f32 = self.cached_tiles.iter()
                            .map(|t| t.coverage)
                            .sum::<f32>() / self.cached_tiles.len() as f32;

                        ui.label(format!("Avg Identity: {:.2}%", avg_identity * 100.0));
                        ui.label(format!("Avg Coverage: {:.2}%", avg_coverage * 100.0));
                    }
                }
            });
        }

        // Central panel - main tracks view
        egui::CentralPanel::default().show(ctx, |ui| {
            if let Some(ref error) = self.error_message {
                ui.colored_label(egui::Color32::RED, error);
            }

            // Position indicator
            ui.horizontal(|ui| {
                if let Some(chrom) = &self.current_chrom {
                    ui.label(format!(
                        "{}:{}-{} ({} bp)",
                        chrom,
                        format_bp(self.view_range.0),
                        format_bp(self.view_range.1),
                        format_bp(self.view_range.1 - self.view_range.0)
                    ));
                }
            });

            ui.separator();

            // Main tracks area
            let _available_height = ui.available_height();
            let track_height = 100.0;

            egui::ScrollArea::vertical().show(ui, |ui| {
                // Sequence track (shows ACTG when zoomed in, indel density when zoomed out)
                ui.group(|ui| {
                    ui.label("Sequence / Indels");
                    ui::draw_sequence_track(ui, &self.cached_tiles, self.view_range, track_height, None, None);
                });

                // Identity track
                ui.group(|ui| {
                    ui.label("Identity");
                    ui::draw_identity_track(ui, &self.cached_tiles, self.view_range, track_height);
                });

                // Coverage track
                ui.group(|ui| {
                    ui.label("Coverage");
                    ui::draw_coverage_track(ui, &self.cached_tiles, self.view_range, track_height);
                });

                // Quality track
                ui.group(|ui| {
                    ui.label("Quality Score");
                    ui::draw_quality_track(ui, &self.cached_tiles, self.view_range, track_height);
                });

                // Gap rate track
                ui.group(|ui| {
                    ui.label("Gap Rate");
                    ui::draw_gap_track(ui, &self.cached_tiles, self.view_range, track_height);
                });
            });
        });

        // Handle keyboard shortcuts
        ctx.input(|i| {
            if i.key_pressed(egui::Key::ArrowLeft) { self.pan_left(); }
            if i.key_pressed(egui::Key::ArrowRight) { self.pan_right(); }
            if i.key_pressed(egui::Key::ArrowUp) || i.key_pressed(egui::Key::Plus) { self.zoom_in(); }
            if i.key_pressed(egui::Key::ArrowDown) || i.key_pressed(egui::Key::Minus) { self.zoom_out(); }
        });
    }
}

fn parse_position(pos: &str) -> Option<(u64, u64)> {
    // Parse formats like "chr1:1000000-2000000" or "1000000-2000000"
    let pos = pos.trim();

    let range_part = if pos.contains(':') {
        pos.split(':').nth(1)?
    } else {
        pos
    };

    let parts: Vec<&str> = range_part.split('-').collect();
    if parts.len() == 2 {
        let start = parse_bp(parts[0])?;
        let end = parse_bp(parts[1])?;
        Some((start, end))
    } else if parts.len() == 1 {
        let center = parse_bp(parts[0])?;
        Some((center.saturating_sub(500_000), center + 500_000))
    } else {
        None
    }
}

fn parse_bp(s: &str) -> Option<u64> {
    let s = s.trim().to_lowercase();

    if s.ends_with("mb") || s.ends_with("m") {
        let num: f64 = s.trim_end_matches(|c| c == 'm' || c == 'b').parse().ok()?;
        Some((num * 1_000_000.0) as u64)
    } else if s.ends_with("kb") || s.ends_with("k") {
        let num: f64 = s.trim_end_matches(|c| c == 'k' || c == 'b').parse().ok()?;
        Some((num * 1_000.0) as u64)
    } else {
        s.replace(",", "").parse().ok()
    }
}

fn format_bp(bp: u64) -> String {
    if bp >= 1_000_000 {
        format!("{:.2}Mb", bp as f64 / 1_000_000.0)
    } else if bp >= 1_000 {
        format!("{:.2}kb", bp as f64 / 1_000.0)
    } else {
        format!("{}bp", bp)
    }
}
